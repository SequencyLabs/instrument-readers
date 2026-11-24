//! Minimal pure-Rust ZIP central-directory reader supporting Stored and Deflate entries.

use std::borrow::Cow;
use std::io::Read;

use flate2::read::DeflateDecoder;

use crate::binary::{read_u16_le, read_u32_le};
use crate::error::{InstrumentError, Result};
use serde::{Deserialize, Serialize};

/// Four-byte local file header magic (`PK\x03\x04`) or empty archive magic (`PK\x05\x06`).
pub const ZIP_LOCAL_MAGIC: [u8; 4] = [0x50, 0x4B, 0x03, 0x04];
pub const ZIP_EOCD_MAGIC: [u8; 4] = [0x50, 0x4B, 0x05, 0x06];
const ZIP_CDFH_MAGIC: u32 = 0x0201_4B50;
const ZIP_LFH_MAGIC: u32 = 0x0403_4B50;
const ZIP_EOCD_SIG: u32 = 0x0605_4B50;

const METHOD_STORED: u16 = 0;
const METHOD_DEFLATED: u16 = 8;
const MAX_UNCOMPRESSED_BYTES: usize = 256 * 1024 * 1024;

/// Metadata for a single member inside a ZIP container (`.dx` or `.rx`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ZipEntryInfo {
    /// Archive-relative member path, such as `"Base/InjectionACAML"` or `"12af0742.CH"`.
    pub name: String,
    /// Compressed payload size in bytes.
    pub compressed_size: usize,
    /// Uncompressed payload size in bytes.
    pub uncompressed_size: usize,
    /// `true` when deflated (`method == 8`), `false` when stored (`method == 0`).
    pub is_compressed: bool,
}

#[derive(Debug, Clone)]
struct RawZipEntry {
    info: ZipEntryInfo,
    method: u16,
    crc32: u32,
    local_header_offset: usize,
}

/// Borrowed reader over an in-memory ZIP archive.
#[derive(Debug)]
pub struct ZipArchive<'a> {
    data: &'a [u8],
    entries: Vec<RawZipEntry>,
}

impl<'a> ZipArchive<'a> {
    /// Parses the End of Central Directory and Central Directory records of `data`.
    pub fn parse(data: &'a [u8]) -> Result<Self> {
        let eocd_pos = find_eocd(data).ok_or(InstrumentError::InvalidZip(
            "missing End of Central Directory record",
        ))?;

        let total_entries = read_u16_le(data, eocd_pos + 10)
            .ok_or(InstrumentError::InvalidZip("truncated EOCD"))?
            as usize;
        let cd_size = read_u32_le(data, eocd_pos + 12)
            .ok_or(InstrumentError::InvalidZip("truncated EOCD"))? as usize;
        let cd_offset = read_u32_le(data, eocd_pos + 16)
            .ok_or(InstrumentError::InvalidZip("truncated EOCD"))? as usize;

        let cd_end = cd_offset
            .checked_add(cd_size)
            .filter(|&end| end <= eocd_pos)
            .ok_or(InstrumentError::InvalidZip(
                "central directory exceeds archive bounds",
            ))?;

        let mut entries = Vec::with_capacity(total_entries);
        let mut cursor = cd_offset;

        for _ in 0..total_entries {
            if cursor + 46 > cd_end {
                return Err(InstrumentError::InvalidZip(
                    "truncated central directory header",
                ));
            }
            if read_u32_le(data, cursor) != Some(ZIP_CDFH_MAGIC) {
                return Err(InstrumentError::InvalidZip(
                    "invalid central directory file signature",
                ));
            }
            let method = read_u16_le(data, cursor + 10).unwrap_or(0);
            let crc32 = read_u32_le(data, cursor + 16).unwrap_or(0);
            let compressed_size = read_u32_le(data, cursor + 20).unwrap_or(0) as usize;
            let uncompressed_size = read_u32_le(data, cursor + 24).unwrap_or(0) as usize;
            let name_len = read_u16_le(data, cursor + 28).unwrap_or(0) as usize;
            let extra_len = read_u16_le(data, cursor + 30).unwrap_or(0) as usize;
            let comment_len = read_u16_le(data, cursor + 32).unwrap_or(0) as usize;
            let local_header_offset = read_u32_le(data, cursor + 42).unwrap_or(0) as usize;

            let name_end = cursor
                .checked_add(46)
                .and_then(|n| n.checked_add(name_len))
                .filter(|&e| e <= cd_end)
                .ok_or(InstrumentError::InvalidZip(
                    "entry name exceeds central directory",
                ))?;
            let next_cursor = name_end
                .checked_add(extra_len)
                .and_then(|n| n.checked_add(comment_len))
                .filter(|&e| e <= cd_end)
                .ok_or(InstrumentError::InvalidZip(
                    "entry extra field exceeds central directory",
                ))?;

            let name = String::from_utf8_lossy(&data[cursor + 46..name_end]).into_owned();
            entries.push(RawZipEntry {
                info: ZipEntryInfo {
                    name,
                    compressed_size,
                    uncompressed_size,
                    is_compressed: method != METHOD_STORED,
                },
                method,
                crc32,
                local_header_offset,
            });
            cursor = next_cursor;
        }

        Ok(Self { data, entries })
    }

    /// Iterates over public metadata for every entry in central-directory order.
    pub fn entries(&self) -> impl Iterator<Item = &ZipEntryInfo> {
        self.entries.iter().map(|e| &e.info)
    }

    /// Collects owned [`ZipEntryInfo`] records for every entry in the archive.
    #[must_use]
    pub fn entry_infos(&self) -> Vec<ZipEntryInfo> {
        self.entries.iter().map(|e| e.info.clone()).collect()
    }

    /// Reads an entry by exact name, returning `Ok(None)` when no entry has that name.
    pub fn read_by_name(&self, target: &str) -> Result<Option<Cow<'a, [u8]>>> {
        let Some(entry) = self.entries.iter().find(|e| e.info.name == target) else {
            return Ok(None);
        };
        self.extract_entry(entry).map(Some)
    }

    /// Reads an entry by case-insensitive name match.
    pub fn read_by_name_ignore_ascii_case(&self, target: &str) -> Result<Option<Cow<'a, [u8]>>> {
        let Some(entry) = self
            .entries
            .iter()
            .find(|e| e.info.name.eq_ignore_ascii_case(target))
        else {
            return Ok(None);
        };
        self.extract_entry(entry).map(Some)
    }

    fn extract_entry(&self, entry: &RawZipEntry) -> Result<Cow<'a, [u8]>> {
        if entry.info.uncompressed_size > MAX_UNCOMPRESSED_BYTES {
            return Err(InstrumentError::InvalidZip(
                "entry uncompressed size exceeds safety cap",
            ));
        }

        let lfh = entry.local_header_offset;
        if lfh.checked_add(30).is_none_or(|end| end > self.data.len()) {
            return Err(InstrumentError::InvalidZip("truncated local file header"));
        }
        if read_u32_le(self.data, lfh) != Some(ZIP_LFH_MAGIC) {
            return Err(InstrumentError::InvalidZip(
                "invalid local file header signature",
            ));
        }

        let lfh_name_len = read_u16_le(self.data, lfh + 26).unwrap_or(0) as usize;
        let lfh_extra_len = read_u16_le(self.data, lfh + 28).unwrap_or(0) as usize;
        let payload_start = lfh
            .checked_add(30)
            .and_then(|n| n.checked_add(lfh_name_len))
            .and_then(|n| n.checked_add(lfh_extra_len))
            .ok_or(InstrumentError::InvalidZip("local header offset overflow"))?;
        let payload_end = payload_start
            .checked_add(entry.info.compressed_size)
            .filter(|&end| end <= self.data.len())
            .ok_or(InstrumentError::InvalidZip(
                "compressed payload exceeds archive bounds",
            ))?;
        let payload = &self.data[payload_start..payload_end];

        let out = match entry.method {
            METHOD_STORED => {
                if payload.len() != entry.info.uncompressed_size {
                    return Err(InstrumentError::InvalidZip("stored entry length mismatch"));
                }
                Cow::Borrowed(payload)
            }
            METHOD_DEFLATED => {
                let cap = entry.info.uncompressed_size;
                let mut decoded = Vec::with_capacity(cap.min(1024 * 1024));
                let limit = u64::try_from(cap).unwrap_or(u64::MAX).saturating_add(1);
                let mut decoder = DeflateDecoder::new(payload).take(limit);
                decoder
                    .read_to_end(&mut decoded)
                    .map_err(|_| InstrumentError::InvalidZip("deflate decompression failed"))?;
                if decoded.len() != cap {
                    return Err(InstrumentError::InvalidZip(
                        "deflated output length does not match header",
                    ));
                }
                Cow::Owned(decoded)
            }
            _ => {
                return Err(InstrumentError::InvalidZip(
                    "unsupported ZIP compression method",
                ));
            }
        };

        if crc32fast::hash(&out) != entry.crc32 {
            return Err(InstrumentError::InvalidZip("entry CRC-32 mismatch"));
        }
        Ok(out)
    }
}

fn find_eocd(data: &[u8]) -> Option<usize> {
    if data.len() < 22 {
        return None;
    }
    let min_pos = data.len().saturating_sub(22 + 65_535);
    let mut pos = data.len() - 22;
    loop {
        if read_u32_le(data, pos) == Some(ZIP_EOCD_SIG) {
            let comment_len = read_u16_le(data, pos + 20).unwrap_or(0) as usize;
            if pos + 22 + comment_len == data.len() {
                return Some(pos);
            }
        }
        if pos == min_pos {
            break;
        }
        pos -= 1;
    }
    None
}
