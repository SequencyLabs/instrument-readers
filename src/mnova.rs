//! Pure-Rust binary container reader for Mestrelab MestReNova (`.mnova`) documents.
//!
//! Parses Qt big-endian data streams, embedded zlib file tables (`flate2`), parameter table rows,
//! 1D NMR axis records, Bruker `acqus` experiment metadata, measurement byte spans (JPEG 2000
//! frames and FIDs), and framed 1D JPEG 2000 codestreams ready for image decoding.

use std::collections::{BTreeMap, BTreeSet};

use flate2::{Decompress, FlushDecompress, Status};

use crate::binary::{read_f32_be, read_f64_be, read_u16_be, read_u32_be};
use crate::bruker::{self, AcqusParams};
use crate::error::{InstrumentError, Result};

const NULL_STR: u32 = 0xFFFF_FFFF;
const ENTRY_META_LEN: usize = 39;
const MAX_ZLIB_STREAM: usize = 8 << 20;
const AXIS_VALUES_LEN: usize = 7 * 8;
const J2K_SOC_SIZ: [u8; 4] = [0xFF, 0x4F, 0xFF, 0x51];
const ROW_PREFIX: [u8; 12] = [0, 0, 0, 0x0D, 0, 0, 0, 0x05, 0, 0, 0, 0];

/// One file embedded in a MestReNova document's zlib file table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MnovaEmbedded {
    pub name: String,
    pub entry: (usize, usize),
    pub start: usize,
    pub end: usize,
    pub path: (usize, usize),
    pub after: (usize, usize),
    pub text: Vec<u8>,
}

/// One parameter row from a MestReNova parameter table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MnovaRow {
    pub span: (usize, usize),
    pub label: Option<String>,
    pub key: Option<String>,
    pub formula: Option<String>,
    pub value: (usize, usize),
}

/// One 1D NMR axis record in a MestReNova document.
#[derive(Debug, Clone, PartialEq)]
pub struct MnovaAxis {
    pub at: usize,
    pub low: f64,
    pub mhz: f64,
    pub width: f64,
    pub nucleus: String,
}

impl MnovaAxis {
    /// Returns `n` chemical shift values in ppm spread across the axis from high to low.
    #[must_use]
    pub fn ppm(&self, n: usize) -> Vec<f64> {
        if n == 0 || self.mhz == 0.0 {
            return Vec::new();
        }
        let step = self.width / (n as f64);
        (0..n)
            .map(|i| (self.low + self.width - (i as f64) * step) / self.mhz)
            .collect()
    }
}

/// One framed 1D JPEG 2000 processed spectrum record extracted from a MestReNova document.
#[derive(Debug, Clone, PartialEq)]
pub struct MnovaJ2kFrame {
    pub at: usize,
    pub position: Option<usize>,
    pub nucleus: String,
    pub low_hz: f64,
    pub mhz: f64,
    pub width_hz: f64,
    pub low_val: f32,
    pub high_val: f32,
    pub depth: u8,
    pub codestream: Vec<u8>,
}

/// One Bruker experiment embedded in a MestReNova document.
#[derive(Debug, Clone, PartialEq)]
pub struct MnovaSpectrum {
    pub dataset: String,
    pub expno: String,
    pub position: usize,
    pub acqus: AcqusParams,
}

struct EmbeddedExperiment {
    at: usize,
    position: usize,
    key: Option<(String, String)>,
    text: String,
}

enum InflateError {
    Corrupt,
    TooLarge,
}

/// Scans `content` for every complete zlib stream, returning `(start, end, decompressed)`.
pub fn zlib_streams(content: &[u8]) -> Vec<(usize, usize, Vec<u8>)> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i + 2 <= content.len() {
        if content[i] == 0x78 && matches!(content[i + 1], 0x01 | 0x5E | 0x9C | 0xDA) {
            let window_end = content.len().min(i.saturating_add(MAX_ZLIB_STREAM));
            if let Ok((consumed, bytes)) = inflate_zlib(&content[i..window_end], MAX_ZLIB_STREAM)
                && !bytes.is_empty()
            {
                out.push((i, i + consumed, bytes));
            }
        }
        i += 1;
    }
    out
}

/// Parses every embedded file table in a MestReNova document.
pub fn parse_embedded(content: &[u8]) -> Result<Vec<MnovaEmbedded>> {
    let scanned = zlib_streams(content);
    let mut cache: BTreeMap<usize, (usize, Vec<u8>)> = scanned
        .into_iter()
        .map(|(start, end, bytes)| (start, (end, bytes)))
        .collect();
    let starts: Vec<usize> = cache.keys().copied().collect();

    let mut found: Vec<MnovaEmbedded> = Vec::new();
    for start in starts {
        if found.iter().any(|e| e.start == start) {
            continue;
        }
        let name_end = start
            .checked_sub(8 + ENTRY_META_LEN)
            .ok_or_else(|| invalid(format!("the stream at byte {start} is in no file table")))?;
        let name = string_before(content, name_end)
            .ok_or_else(|| invalid(format!("the stream at byte {start} is in no file table")))?;
        let count_off = name
            .0
            .checked_sub(8)
            .ok_or_else(|| invalid(format!("the stream at byte {start} is in no file table")))?;
        let count = read_u32_be(content, count_off)
            .ok_or_else(|| invalid(format!("the stream at byte {start} is in no file table")))?
            as usize;
        if count == 0 || count >= 1024 {
            return Err(invalid(format!(
                "the stream at byte {start} is in no file table"
            )));
        }
        let mut at = name.0 - 4;
        for _ in 0..count {
            let entry = parse_entry(content, at, &mut cache)?;
            at = entry.entry.1;
            found.push(entry);
        }
        if !found.iter().any(|e| e.start == start) {
            return Err(invalid(format!(
                "the stream at byte {start} is in no file table"
            )));
        }
    }
    Ok(found)
}

/// Parses every parameter table row in a MestReNova document.
pub fn parse_rows(content: &[u8]) -> Result<Vec<MnovaRow>> {
    let mut found = Vec::new();
    let mut i = 0usize;
    while i + 13 <= content.len() {
        if content[i..i + 12] == ROW_PREFIX && matches!(content[i + 12], 0 | 1) {
            let at = i + 17;
            let shows = match (
                read_u32_be(content, at),
                read_u32_be(content, at + 4),
                read_u32_be(content, at + 8),
                read_u32_be(content, at + 12),
            ) {
                (Some(1), Some(0), _, _) => at + 8,
                (Some(1), Some(8), Some(NULL_STR), Some(0)) => at + 16,
                _ => {
                    i += 1;
                    continue;
                }
            };
            let source = read_qt_string(content, shows)?;
            let (key, formula, value_at) = if source.0 == source.1 {
                let f = read_qt_string(content, source.1)?;
                (None, Some(f), f.1)
            } else {
                let k = read_qt_string(content, source.1)?;
                (Some(k), None, k.1)
            };
            let value = read_qt_string(content, value_at + 1)?;
            let label = read_qt_string(content, value.1 + 1)?;
            let label_str = decode_utf16be(&content[label.0..label.1]);
            found.push(MnovaRow {
                span: (i, label.1),
                label: (!label_str.is_empty()).then_some(label_str),
                key: key.map(|(a, b)| decode_utf16be(&content[a..b])),
                formula: formula.map(|(a, b)| decode_utf16be(&content[a..b])),
                value,
            });
            i = label.1;
        } else {
            i += 1;
        }
    }
    Ok(found)
}

/// Finds every measurement byte span (framed JPEG 2000 spectra and raw FID blocks).
#[must_use]
pub fn parse_measurements(content: &[u8]) -> Vec<(usize, usize)> {
    let mut spans = Vec::new();
    for at in find_subslice(content, &J2K_SOC_SIZ) {
        if let Some(end) = framed_end(content, at) {
            spans.push((at - 12, end));
        }
    }
    let mut i = 0usize;
    while i + 19 <= content.len() {
        if content[i + 4..i + 11] == [0, 0, 0, 0, 1, 0, 0]
            && matches!(content[i + 11], 2 | 4 | 8 | 32)
            && content[i + 12..i + 15] == [0, 5, 0]
            && content[i + 15] == content[i + 11]
            && content[i + 16..i + 19] == [0, 0, 1]
            && let Some(size) = read_u32_be(content, i)
        {
            let size = size as usize;
            let width = content[i + 11] as usize;
            if size >= 15
                && (size - 15).is_multiple_of(width)
                && let Some(end) = i.checked_add(4).and_then(|s| s.checked_add(size))
                && end <= content.len()
            {
                spans.push((i, end));
            }
        }
        i += 1;
    }
    spans.sort_unstable();
    spans
}

/// Extracts every valid 1D NMR axis record in `content`.
#[must_use]
pub fn parse_axes(content: &[u8]) -> Vec<MnovaAxis> {
    let mut found = Vec::new();
    let mut i = AXIS_VALUES_LEN;
    while i + 5 <= content.len() {
        if content[i..i + 4] == [1, 0, 0, 0]
            && (4..=10).contains(&content[i + 4])
            && content[i + 4].is_multiple_of(2)
        {
            let nbytes = content[i + 4] as usize;
            let str_start = i + 5;
            let str_end = str_start + nbytes;
            if str_end <= content.len()
                && let Some(nucleus) = parse_nucleus_utf16be(&content[str_start..str_end])
            {
                let at = i - AXIS_VALUES_LEN;
                if let (Some(low), Some(mhz), Some(width)) = (
                    read_f64_be(content, at),
                    read_f64_be(content, at + 8),
                    read_f64_be(content, at + 16),
                ) && mhz > 1.0
                    && mhz < 2000.0
                    && width > 0.0
                    && width < 1e7
                    && low.abs() < 1e7
                {
                    found.push(MnovaAxis {
                        at,
                        low,
                        mhz,
                        width,
                        nucleus,
                    });
                }
            }
        }
        i += 1;
    }
    found
}

/// Extracts every embedded Bruker `acqus` experiment in document order (deduplicated by `(dataset, expno)`).
#[must_use]
pub fn parse_spectra(content: &[u8]) -> Vec<MnovaSpectrum> {
    let mut seen: BTreeSet<(String, String)> = BTreeSet::new();
    let mut out = Vec::new();
    for exp in find_experiments(content) {
        if let Some((dataset, expno)) = exp.key
            && seen.insert((dataset.clone(), expno.clone()))
        {
            out.push(MnovaSpectrum {
                dataset,
                expno,
                position: exp.position,
                acqus: bruker::parse_acqus_fields(&exp.text),
            });
        }
    }
    out
}

/// Extracts every framed 1D JPEG 2000 spectrum codestream paired with its axis and experiment position.
#[must_use]
pub fn parse_j2k_frames(content: &[u8]) -> Vec<MnovaJ2kFrame> {
    let axes = parse_axes(content);
    let experiments = find_experiments(content);
    let mut frames = Vec::new();

    for at in find_subslice(content, &J2K_SOC_SIZ) {
        let Some(axis) = axes.iter().rev().find(|a| a.at < at) else {
            continue;
        };
        if framed_end(content, at).is_none() {
            continue;
        }
        let Some(size) = read_u32_be(content, at - 4).map(|v| v as usize) else {
            continue;
        };
        let (Some(low_val), Some(high_val), Some(rows), Some(components), Some(&depth)) = (
            read_f32_be(content, at - 12),
            read_f32_be(content, at - 8),
            read_u32_be(content, at + 12),
            read_u16_be(content, at + 40),
            content.get(at + 42),
        ) else {
            continue;
        };
        if rows != 1
            || components != 1
            || (depth & 0x7F) == 0
            || !low_val.is_finite()
            || !high_val.is_finite()
        {
            continue;
        }
        let Some(codestream) = content.get(at..at + size) else {
            continue;
        };
        let position = experiments
            .iter()
            .rev()
            .find(|exp| exp.at < axis.at)
            .map(|exp| exp.position);
        frames.push(MnovaJ2kFrame {
            at,
            position,
            nucleus: axis.nucleus.clone(),
            low_hz: axis.low,
            mhz: axis.mhz,
            width_hz: axis.width,
            low_val,
            high_val,
            depth,
            codestream: codestream.to_vec(),
        });
    }

    frames
}

fn find_experiments(content: &[u8]) -> Vec<EmbeddedExperiment> {
    let mut out = Vec::new();
    let mut position = 0usize;
    for (at, _, stream) in zlib_streams(content) {
        if stream.starts_with(b"##TITLE= Parameter file")
            && stream.windows(8).any(|w| w == b"##$NUC1=")
        {
            let text: String = stream.iter().map(|&b| char::from(b)).collect();
            let key = extract_acqus_dataset_expno(&text);
            out.push(EmbeddedExperiment {
                at,
                position,
                key,
                text,
            });
            position += 1;
        }
    }
    out
}

fn extract_acqus_dataset_expno(text: &str) -> Option<(String, String)> {
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("$$ ") else {
            continue;
        };
        let mut parts = rest.trim_end().rsplit(['/', '\\']);
        let (Some(file), Some(expno), Some(dataset)) = (parts.next(), parts.next(), parts.next())
        else {
            continue;
        };
        if file == "acqus"
            && !dataset.is_empty()
            && !expno.is_empty()
            && expno.bytes().all(|b| b.is_ascii_digit())
        {
            return Some((dataset.to_owned(), expno.to_owned()));
        }
    }
    None
}

fn parse_entry(
    content: &[u8],
    at: usize,
    cache: &mut BTreeMap<usize, (usize, Vec<u8>)>,
) -> Result<MnovaEmbedded> {
    let name = read_qt_string(content, at)?;
    let head = name.1;
    let magic_ok = matches!(
        (
            read_u32_be(content, head),
            read_u32_be(content, head + 4),
            read_u32_be(content, head + ENTRY_META_LEN - 4),
        ),
        (Some(0x1F), Some(0), Some(0))
    );
    if !magic_ok {
        return Err(invalid(format!(
            "the file table's entry at byte {at} is of a kind no reader knows"
        )));
    }
    let size = read_u32_be(content, head + ENTRY_META_LEN)
        .ok_or_else(|| invalid("truncated embedded entry size"))? as usize;
    let length = read_u32_be(content, head + ENTRY_META_LEN + 4)
        .ok_or_else(|| invalid("truncated embedded entry length"))? as usize;
    let start = head + ENTRY_META_LEN + 8;
    let (end, text) = if size == 4 && length == 0 {
        (start, Vec::new())
    } else {
        if length > MAX_ZLIB_STREAM {
            return Err(invalid(format!(
                "the file at byte {start} exceeds the embedded size limit"
            )));
        }
        let comp_len = size
            .checked_sub(4)
            .ok_or_else(|| invalid(format!("the file at byte {start} is not its frame's")))?;
        let end = start
            .checked_add(comp_len)
            .filter(|&e| e <= content.len())
            .ok_or_else(|| invalid(format!("the file at byte {start} is not its frame's")))?;

        let decoded = if let Some((cached_end, cached_bytes)) = cache.remove(&start) {
            if cached_end != end || cached_bytes.len() != length {
                return Err(invalid(format!(
                    "the file at byte {start} is not its frame's"
                )));
            }
            cached_bytes
        } else {
            match inflate_zlib(&content[start..end], length.max(1)) {
                Ok((consumed, decoded)) if consumed == comp_len && decoded.len() == length => {
                    decoded
                }
                Ok(_) | Err(InflateError::TooLarge) => {
                    return Err(invalid(format!(
                        "the file at byte {start} is not its frame's"
                    )));
                }
                Err(InflateError::Corrupt) => {
                    return Err(invalid(format!("corrupt zlib stream at byte {start}")));
                }
            }
        };
        (end, decoded)
    };
    let path = read_qt_string(content, end)?;
    let after = read_qt_string(content, path.1)?;
    Ok(MnovaEmbedded {
        name: decode_utf16be(&content[name.0..name.1]),
        entry: (at, after.1),
        start,
        end,
        path,
        after,
        text,
    })
}

fn read_qt_string(content: &[u8], at: usize) -> Result<(usize, usize)> {
    let size = read_u32_be(content, at)
        .ok_or_else(|| invalid(format!("the document ends before byte {}", at + 4)))?;
    if size == NULL_STR {
        return Ok((at + 4, at + 4));
    }
    let size = size as usize;
    let end = at
        .checked_add(4)
        .and_then(|s| s.checked_add(size))
        .ok_or_else(|| invalid(format!("no string at byte {at}")))?;
    if !size.is_multiple_of(2) || end > content.len() {
        return Err(invalid(format!("no string at byte {at}")));
    }
    Ok((at + 4, end))
}

fn string_before(content: &[u8], end: usize) -> Option<(usize, usize)> {
    let mut at = end.checked_sub(6)?;
    let floor = end.saturating_sub(136);
    loop {
        if let Some(len) = read_u32_be(content, at)
            && (len as usize) == end - at - 4
        {
            return Some((at + 4, end));
        }
        if at < floor + 2 {
            break;
        }
        at -= 2;
    }
    None
}

fn framed_end(content: &[u8], at: usize) -> Option<usize> {
    if at < 12 || at + 43 > content.len() {
        return None;
    }
    let size = read_u32_be(content, at - 4)? as usize;
    let tail_start = at.checked_add(size)?;
    let tail_end = tail_start.checked_add(8)?;
    if tail_end > content.len() {
        return None;
    }
    (content[tail_start..tail_end] == content[at - 12..at - 4]).then_some(tail_end)
}

fn parse_nucleus_utf16be(bytes: &[u8]) -> Option<String> {
    let s = decode_utf16be(bytes);
    let b = s.as_bytes();
    let digits = b.iter().take_while(|c| c.is_ascii_digit()).count();
    if !(1..=3).contains(&digits) || digits >= b.len() {
        return None;
    }
    let sym = &b[digits..];
    let valid_symbol = (sym.len() == 1 && sym[0].is_ascii_uppercase())
        || (sym.len() == 2 && sym[0].is_ascii_uppercase() && sym[1].is_ascii_lowercase());
    valid_symbol.then_some(s)
}

fn decode_utf16be(bytes: &[u8]) -> String {
    let (cow, _) = encoding_rs::UTF_16BE.decode_without_bom_handling(bytes);
    cow.into_owned()
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Vec<usize> {
    if needle.is_empty() || haystack.len() < needle.len() {
        return Vec::new();
    }
    haystack
        .windows(needle.len())
        .enumerate()
        .filter_map(|(i, w)| (w == needle).then_some(i))
        .collect()
}

fn inflate_zlib(
    input: &[u8],
    max_out: usize,
) -> std::result::Result<(usize, Vec<u8>), InflateError> {
    let mut decomp = Decompress::new(true);
    let mut out = Vec::new();
    let mut buf = [0u8; 8192];
    loop {
        let before_in = decomp.total_in() as usize;
        let before_out = decomp.total_out() as usize;
        let status = decomp
            .decompress(&input[before_in..], &mut buf, FlushDecompress::None)
            .map_err(|_| InflateError::Corrupt)?;
        let produced = (decomp.total_out() as usize) - before_out;
        if produced > max_out.saturating_sub(out.len()) {
            return Err(InflateError::TooLarge);
        }
        out.extend_from_slice(&buf[..produced]);
        match status {
            Status::StreamEnd => return Ok((decomp.total_in() as usize, out)),
            Status::Ok => {
                if produced == 0 && (decomp.total_in() as usize) == before_in {
                    return Err(InflateError::Corrupt);
                }
            }
            Status::BufError => return Err(InflateError::Corrupt),
        }
    }
}

fn invalid(msg: impl Into<String>) -> InstrumentError {
    InstrumentError::InvalidMnova(msg.into())
}
