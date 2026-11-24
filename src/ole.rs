//! Pure-Rust OLE2 Compound File Binary Format (MS-CFB) reader and sector-layout mapper.

use std::borrow::Cow;
use std::collections::BTreeMap;

use crate::binary::{decode_utf16le, read_u16_le, read_u32_le, read_u64_le};
use crate::error::{InstrumentError, Result};

/// Eight-byte magic signature at byte 0 of every OLE2 Compound File (`D0 CF 11 E0 A1 B1 1A E1`).
pub const OLE_MAGIC: [u8; 8] = [0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1];

const MAX_REG_SECTOR: u32 = 0xFFFF_FFFA;
const DIFAT_SECTOR: u32 = 0xFFFF_FFFC;
const FAT_SECTOR: u32 = 0xFFFF_FFFD;
const NO_STREAM: u32 = 0xFFFF_FFFF;

const ENTRY_STORAGE: u8 = 1;
const ENTRY_STREAM: u8 = 2;
const ENTRY_ROOT: u8 = 5;

const HEADER_LEN: usize = 512;
const DIR_ENTRY_LEN: usize = 128;

/// Metadata and file byte extents for a single stream inside an [`OleFile`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OleStreamInfo {
    /// Slash-joined storage path, such as `"PDA 3D Raw Data/Max Plot"`.
    pub name: String,
    /// Logical byte length of the stream.
    pub size: usize,
    /// Merged `(start, end)` byte ranges in the outer compound file holding this stream's bytes.
    pub extents: Vec<(usize, usize)>,
}

#[derive(Debug, Clone)]
struct DirEntry {
    name: String,
    entry_type: u8,
    left_sibling: u32,
    right_sibling: u32,
    child: u32,
    start_sector: u32,
    size: usize,
}

/// Borrowed view over an OLE2 Compound File Binary Format buffer.
#[derive(Debug)]
pub struct OleFile<'a> {
    data: &'a [u8],
    streams: BTreeMap<String, OleStreamInfo>,
    ordered_streams: Vec<OleStreamInfo>,
    unallocated: Vec<(usize, usize)>,
}

impl<'a> OleFile<'a> {
    /// Parses an OLE2 compound file's header, DIFAT, FAT, MiniFAT, directory tree, and byte extents.
    pub fn parse(data: &'a [u8]) -> Result<Self> {
        if data.len() < HEADER_LEN || data[..8] != OLE_MAGIC {
            return Err(InstrumentError::InvalidOle("missing OLE2 header signature"));
        }
        if read_u16_le(data, 28) != Some(0xFFFE) {
            return Err(InstrumentError::InvalidOle("invalid OLE2 byte-order mark"));
        }

        let sector_shift =
            read_u16_le(data, 30).ok_or(InstrumentError::InvalidOle("short header"))?;
        if sector_shift != 9 && sector_shift != 12 {
            return Err(InstrumentError::InvalidOle("unsupported sector shift"));
        }
        let mini_shift =
            read_u16_le(data, 32).ok_or(InstrumentError::InvalidOle("short header"))?;
        if mini_shift != 6 {
            return Err(InstrumentError::InvalidOle("unsupported mini-sector shift"));
        }

        let sector_size = 1usize << sector_shift;
        let mini_size = 1usize << mini_shift;
        if data.len() < sector_size {
            return Err(InstrumentError::InvalidOle("file shorter than one sector"));
        }

        let num_fat_sectors = read_u32_le(data, 44).unwrap_or(0) as usize;
        let first_dir_sector = read_u32_le(data, 48).unwrap_or(NO_STREAM);
        let mini_cutoff = read_u32_le(data, 56).unwrap_or(4096) as usize;
        let first_minifat_sector = read_u32_le(data, 60).unwrap_or(NO_STREAM);
        let first_difat_sector = read_u32_le(data, 68).unwrap_or(NO_STREAM);
        let num_difat_sectors = read_u32_le(data, 72).unwrap_or(0) as usize;

        let fat_sector_ids = collect_fat_sectors(
            data,
            sector_size,
            num_fat_sectors,
            first_difat_sector,
            num_difat_sectors,
        )?;
        let fat = read_u32_table(data, sector_size, &fat_sector_ids)?;
        if fat.is_empty() {
            return Err(InstrumentError::InvalidOle("empty FAT"));
        }

        let dir_sectors = follow_chain(&fat, first_dir_sector)?;
        if dir_sectors.is_empty() {
            return Err(InstrumentError::InvalidOle("empty directory chain"));
        }
        let entries = parse_directory(data, sector_size, sector_shift, &dir_sectors)?;
        let root = entries
            .first()
            .filter(|e| e.entry_type == ENTRY_ROOT)
            .ok_or(InstrumentError::InvalidOle("missing root directory entry"))?;

        let (container_offsets, minifat_sectors, minifat) = if root.size > 0 {
            let root_chain = follow_chain(&fat, root.start_sector)?;
            let offsets = sector_offsets(data.len(), sector_size, &root_chain)?;
            let mf_sectors = if first_minifat_sector < MAX_REG_SECTOR {
                follow_chain(&fat, first_minifat_sector)?
            } else {
                Vec::new()
            };
            let mf = read_u32_table(data, sector_size, &mf_sectors)?;
            (offsets, mf_sectors, mf)
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };

        let ctx = ExtentContext {
            file_len: data.len(),
            sector_size,
            mini_size,
            mini_cutoff,
            fat: &fat,
            minifat: &minifat,
            container_offsets: &container_offsets,
        };

        let mut ordered_streams = Vec::new();
        let mut visited = vec![false; entries.len()];
        visited[0] = true;
        walk_storage(
            &entries,
            root.child,
            "",
            &ctx,
            &mut visited,
            &mut ordered_streams,
        )?;

        let unallocated = compute_unallocated(
            data.len(),
            sector_size,
            &fat,
            &dir_sectors,
            &minifat_sectors,
            &ordered_streams,
        )?;

        let mut streams = BTreeMap::new();
        for info in &ordered_streams {
            streams.insert(info.name.clone(), info.clone());
        }

        Ok(Self {
            data,
            streams,
            ordered_streams,
            unallocated,
        })
    }

    /// Returns `true` if a stream with the given slash-separated path exists.
    #[must_use]
    pub fn exists(&self, path: &str) -> bool {
        self.streams.contains_key(path)
    }

    /// Returns the logical byte size of `path` if present.
    #[must_use]
    pub fn stream_size(&self, path: &str) -> Option<usize> {
        self.streams.get(path).map(|s| s.size)
    }

    /// Iterates over all stream paths in directory-tree order.
    pub fn stream_names(&self) -> impl Iterator<Item = &str> {
        self.ordered_streams.iter().map(|s| s.name.as_str())
    }

    /// Returns metadata and file extents for every stream in directory-tree order.
    #[must_use]
    pub fn streams(&self) -> &[OleStreamInfo] {
        &self.ordered_streams
    }

    /// Returns the file byte extents `(start, end)` of `path`.
    #[must_use]
    pub fn extents(&self, path: &str) -> Option<&[(usize, usize)]> {
        self.streams.get(path).map(|s| s.extents.as_slice())
    }

    /// Returns byte ranges in the compound file not held by the header, tables, or any stream.
    #[must_use]
    pub fn unallocated(&self) -> &[(usize, usize)] {
        &self.unallocated
    }

    /// Reads a stream's bytes, borrowing directly from the input slice when contiguous.
    #[must_use]
    pub fn stream(&self, path: &str) -> Option<Cow<'a, [u8]>> {
        let info = self.streams.get(path)?;
        match info.extents.as_slice() {
            [] => Some(Cow::Borrowed(&[])),
            &[(start, end)] => Some(Cow::Borrowed(&self.data[start..end])),
            many => {
                let mut buf = Vec::with_capacity(info.size);
                for &(start, end) in many {
                    buf.extend_from_slice(&self.data[start..end]);
                }
                Some(Cow::Owned(buf))
            }
        }
    }

    /// Reads a required stream or returns [`InstrumentError::MissingStream`].
    pub fn require_stream(&self, path: &str) -> Result<Cow<'a, [u8]>> {
        self.stream(path)
            .ok_or_else(|| InstrumentError::MissingStream(path.to_owned()))
    }
}

/// Maps a `[start, end)` slice within a stream to its `(file_start, file_end)` byte runs.
#[must_use]
pub fn place(extents: &[(usize, usize)], start: usize, end: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut offset = 0usize;
    for &(a, b) in extents {
        let span = b.saturating_sub(a);
        let lo = start.max(offset);
        let hi = end.min(offset + span);
        if lo < hi {
            out.push((a + (lo - offset), a + (hi - offset)));
        }
        offset += span;
    }
    out
}

struct ExtentContext<'a> {
    file_len: usize,
    sector_size: usize,
    mini_size: usize,
    mini_cutoff: usize,
    fat: &'a [u32],
    minifat: &'a [u32],
    container_offsets: &'a [usize],
}

fn sector_offset(file_len: usize, sector_size: usize, sector: u32) -> Result<usize> {
    let pos = (sector as usize)
        .checked_add(1)
        .and_then(|n| n.checked_mul(sector_size))
        .ok_or(InstrumentError::CorruptChain(sector))?;
    if pos
        .checked_add(sector_size)
        .is_none_or(|end| end > file_len)
    {
        return Err(InstrumentError::CorruptChain(sector));
    }
    Ok(pos)
}

fn sector_offsets(file_len: usize, sector_size: usize, sectors: &[u32]) -> Result<Vec<usize>> {
    sectors
        .iter()
        .map(|&s| sector_offset(file_len, sector_size, s))
        .collect()
}

fn collect_fat_sectors(
    data: &[u8],
    sector_size: usize,
    num_fat_sectors: usize,
    first_difat_sector: u32,
    num_difat_sectors: usize,
) -> Result<Vec<u32>> {
    let mut out = Vec::with_capacity(num_fat_sectors);
    for idx in 0..109 {
        if out.len() >= num_fat_sectors {
            break;
        }
        if let Some(sec) = read_u32_le(data, 76 + idx * 4)
            && sec < MAX_REG_SECTOR
        {
            out.push(sec);
        }
    }

    let max_sectors = data.len() / sector_size;
    let per_difat = (sector_size / 4).saturating_sub(1);
    let mut difat = first_difat_sector;
    let mut seen = 0usize;

    while difat < MAX_REG_SECTOR && seen < num_difat_sectors && out.len() < num_fat_sectors {
        if seen > max_sectors {
            return Err(InstrumentError::CorruptChain(difat));
        }
        let base = sector_offset(data.len(), sector_size, difat)?;
        for i in 0..per_difat {
            if out.len() >= num_fat_sectors {
                break;
            }
            if let Some(sec) = read_u32_le(data, base + i * 4)
                && sec < MAX_REG_SECTOR
            {
                out.push(sec);
            }
        }
        difat = read_u32_le(data, base + per_difat * 4).unwrap_or(NO_STREAM);
        seen += 1;
    }

    Ok(out)
}

fn read_u32_table(data: &[u8], sector_size: usize, sectors: &[u32]) -> Result<Vec<u32>> {
    let entries_per_sector = sector_size / 4;
    let mut table = Vec::with_capacity(sectors.len() * entries_per_sector);
    for &sec in sectors {
        let base = sector_offset(data.len(), sector_size, sec)?;
        for i in 0..entries_per_sector {
            let val = read_u32_le(data, base + i * 4).ok_or(InstrumentError::CorruptChain(sec))?;
            table.push(val);
        }
    }
    Ok(table)
}

fn follow_chain(table: &[u32], start: u32) -> Result<Vec<u32>> {
    let mut chain = Vec::new();
    let mut cur = start;
    while cur < MAX_REG_SECTOR {
        let &next = table
            .get(cur as usize)
            .ok_or(InstrumentError::CorruptChain(cur))?;
        if chain.len() >= table.len() {
            return Err(InstrumentError::CorruptChain(cur));
        }
        chain.push(cur);
        cur = next;
    }
    Ok(chain)
}

fn parse_directory(
    data: &[u8],
    sector_size: usize,
    sector_shift: u16,
    dir_sectors: &[u32],
) -> Result<Vec<DirEntry>> {
    let per_sector = sector_size / DIR_ENTRY_LEN;
    let mut entries = Vec::with_capacity(dir_sectors.len() * per_sector);

    for &sec in dir_sectors {
        let base = sector_offset(data.len(), sector_size, sec)?;
        for slot_idx in 0..per_sector {
            let slot =
                &data[base + slot_idx * DIR_ENTRY_LEN..base + (slot_idx + 1) * DIR_ENTRY_LEN];
            let raw_name_len = read_u16_le(slot, 64).unwrap_or(0) as usize;
            let name_bytes = raw_name_len.min(64);
            let name = decode_utf16le(&slot[..name_bytes]);
            let entry_type = slot[66];
            let left_sibling = read_u32_le(slot, 68).unwrap_or(NO_STREAM);
            let right_sibling = read_u32_le(slot, 72).unwrap_or(NO_STREAM);
            let child = read_u32_le(slot, 76).unwrap_or(NO_STREAM);
            let start_sector = read_u32_le(slot, 116).unwrap_or(NO_STREAM);
            let size = if sector_shift == 9 {
                read_u32_le(slot, 120).unwrap_or(0) as usize
            } else {
                usize::try_from(read_u64_le(slot, 120).unwrap_or(0))
                    .map_err(|_| InstrumentError::InvalidOle("stream size overflows usize"))?
            };
            entries.push(DirEntry {
                name,
                entry_type,
                left_sibling,
                right_sibling,
                child,
                start_sector,
                size,
            });
        }
    }

    Ok(entries)
}

fn walk_storage(
    entries: &[DirEntry],
    first_child: u32,
    prefix: &str,
    ctx: &ExtentContext<'_>,
    visited: &mut [bool],
    out: &mut Vec<OleStreamInfo>,
) -> Result<()> {
    if first_child == NO_STREAM {
        return Ok(());
    }

    let mut siblings = Vec::new();
    collect_siblings(entries, first_child, visited, &mut siblings)?;

    for idx in siblings {
        let entry = &entries[idx];
        let full_name = if prefix.is_empty() {
            entry.name.clone()
        } else {
            format!("{prefix}/{}", entry.name)
        };

        match entry.entry_type {
            ENTRY_STREAM => {
                let extents = stream_extents(entry.size, entry.start_sector, ctx)?;
                out.push(OleStreamInfo {
                    name: full_name,
                    size: entry.size,
                    extents,
                });
            }
            ENTRY_STORAGE => {
                walk_storage(entries, entry.child, &full_name, ctx, visited, out)?;
            }
            _ => {}
        }
    }

    Ok(())
}

fn collect_siblings(
    entries: &[DirEntry],
    root_id: u32,
    visited: &mut [bool],
    order: &mut Vec<usize>,
) -> Result<()> {
    let mut stack: Vec<(u32, bool)> = vec![(root_id, false)];
    while let Some((id, expanded)) = stack.pop() {
        if id == NO_STREAM {
            continue;
        }
        let idx = id as usize;
        let entry = entries.get(idx).ok_or(InstrumentError::CorruptChain(id))?;

        if expanded {
            order.push(idx);
            continue;
        }
        if visited[idx] {
            return Err(InstrumentError::CorruptChain(id));
        }
        visited[idx] = true;

        stack.push((entry.right_sibling, false));
        stack.push((id, true));
        stack.push((entry.left_sibling, false));
    }
    Ok(())
}

fn stream_extents(
    size: usize,
    start_sector: u32,
    ctx: &ExtentContext<'_>,
) -> Result<Vec<(usize, usize)>> {
    if size == 0 {
        return Ok(Vec::new());
    }
    if size >= ctx.mini_cutoff {
        let chain = follow_chain(ctx.fat, start_sector)?;
        let offsets = sector_offsets(ctx.file_len, ctx.sector_size, &chain)?;
        merge_runs(&offsets, ctx.sector_size, size, ctx.file_len)
    } else {
        let chain = follow_chain(ctx.minifat, start_sector)?;
        let mut offsets = Vec::with_capacity(chain.len());
        for m in chain {
            let byte_pos = (m as usize)
                .checked_mul(ctx.mini_size)
                .ok_or(InstrumentError::CorruptChain(m))?;
            let c_idx = byte_pos / ctx.sector_size;
            let c_rem = byte_pos % ctx.sector_size;
            let base = *ctx
                .container_offsets
                .get(c_idx)
                .ok_or(InstrumentError::CorruptChain(m))?;
            offsets.push(base + c_rem);
        }
        merge_runs(&offsets, ctx.mini_size, size, ctx.file_len)
    }
}

fn merge_runs(
    offsets: &[usize],
    unit_size: usize,
    length: usize,
    file_len: usize,
) -> Result<Vec<(usize, usize)>> {
    let needed = length.div_ceil(unit_size);
    if offsets.len() < needed {
        return Err(InstrumentError::InvalidOle(
            "stream sector chain shorter than declared size",
        ));
    }
    let mut runs: Vec<(usize, usize)> = Vec::new();
    for (i, &at) in offsets[..needed].iter().enumerate() {
        let chunk = unit_size.min(length - i * unit_size);
        let end =
            at.checked_add(chunk)
                .filter(|&e| e <= file_len)
                .ok_or(InstrumentError::InvalidOle(
                    "stream extent exceeds file bounds",
                ))?;
        if let Some(last) = runs.last_mut()
            && last.1 == at
        {
            last.1 = end;
        } else {
            runs.push((at, end));
        }
    }
    Ok(runs)
}

fn compute_unallocated(
    file_len: usize,
    sector_size: usize,
    fat: &[u32],
    dir_sectors: &[u32],
    minifat_sectors: &[u32],
    streams: &[OleStreamInfo],
) -> Result<Vec<(usize, usize)>> {
    let mut held: Vec<(usize, usize)> = Vec::new();
    held.push((0, sector_size));

    for (sec_idx, &entry) in fat.iter().enumerate() {
        if entry == FAT_SECTOR || entry == DIFAT_SECTOR {
            let sec = u32::try_from(sec_idx)
                .map_err(|_| InstrumentError::InvalidOle("sector index overflow"))?;
            if let Ok(at) = sector_offset(file_len, sector_size, sec) {
                held.push((at, at + sector_size));
            }
        }
    }

    for &sec in dir_sectors.iter().chain(minifat_sectors.iter()) {
        let at = sector_offset(file_len, sector_size, sec)?;
        held.push((at, at + sector_size));
    }

    for stream in streams {
        held.extend_from_slice(&stream.extents);
    }

    held.sort_unstable();
    held.push((file_len, file_len));

    let mut free = Vec::new();
    let mut cursor = 0usize;
    for (a, b) in held {
        if a > cursor {
            free.push((cursor, a.min(file_len)));
        }
        cursor = cursor.max(b);
    }
    Ok(free)
}
