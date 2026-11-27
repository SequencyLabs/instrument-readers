//! Pure-Rust reader for Thermo Finnigan v66 `.raw` runs and `.sld` sequence files.
//!
//! Reads the 1,356-byte Finnigan header, audit timestamps and Windows logins, 31-string sample
//! rows, mass spectrometer run header (controller type 0), scan index entries, and centroided
//! `(m/z, intensity)` peak packets.

use std::collections::BTreeSet;

use crate::binary::{
    decode_utf16le, file_name, format_unix_timestamp, read_f32_le, read_f64_le, read_u32_le,
    read_u64_le,
};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Technique, Trace, VendorFormat};

pub const MAGIC: &[u8] = b"\x01\xa1F\0i\0n\0n\0i\0g\0a\0n\0";
const VERSION: u32 = 66;
const HEADER_LEN: usize = 1356;
const INJECTION_LEN: usize = 64;
const MIN_ROW_BYTES: usize = INJECTION_LEN + 31 * 4 + 4;
const LOGIN_OFFSETS: [usize; 4] = [0x30, 0x62, 0xA0, 0xD2];
const LOGIN_BYTES: usize = 50;
const SEQUENCE_LABELS_OFF: usize = 0x5B0;
const CONTROLLERS_OFF: usize = 0x1C;
const CONTROLLER_TABLE_OFF: usize = 0x330;
const INSTRUMENT_ID_OFF: usize = 0x1D98;
const SCAN_INDEX_OFF: usize = 0x1CF0;
const SCAN_ENTRY_LEN: usize = 88;
pub const PEAK_RECORD_LEN: usize = 56;
const FILETIME_EPOCH_DIFF_SECS: i64 = 11_644_473_600;

/// One sample row from a Thermo Finnigan `.raw` run or `.sld` sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermoRow {
    pub sample_name: String,
    pub sample_id: String,
    pub data_file: String,
    pub vial: String,
    pub method: String,
}

/// Parsed Thermo Finnigan `.raw` or `.sld` metadata header and sample rows.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoFile {
    pub first: String,
    pub second: String,
    pub logins: BTreeSet<String>,
    pub rows: Vec<ThermoRow>,
    pub instrument: Option<String>,
    pub scans: u32,
    pub minutes: Option<(f64, f64)>,
}

/// Scan index summary from a Thermo Finnigan `.raw` mass spectrometer controller.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoScans {
    pub minutes: Vec<f64>,
    pub tic: Vec<f64>,
    pub offsets: Vec<u64>,
    pub data_offset: u64,
}

/// Physical byte spans of user-written login slots and UTF-16LE strings in a Thermo Finnigan file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermoTextSpans {
    pub logins: Vec<(usize, usize, usize)>,
    pub strings: Vec<(usize, usize)>,
}

fn parse_run_with_offset(data: &[u8]) -> Result<(ThermoFile, Option<usize>)> {
    let mut file = parse_header(data)?;
    let (row, ms_header) = parse_ms_header(data)?;
    file.rows = vec![row];

    let Some(header_u64) = ms_header else {
        return Ok((file, None));
    };
    let Ok(header) = usize::try_from(header_u64) else {
        return Ok((file, None));
    };

    let Some(min_end) = header.checked_add(0x58) else {
        return Ok((file, Some(header)));
    };
    if min_end <= data.len() {
        let first = read_u32_le(data, header + 8).ok_or_else(truncated)?;
        let last = read_u32_le(data, header + 12).ok_or_else(truncated)?;
        let count_u64 = if last != 0 && last >= first {
            u64::from(last) - u64::from(first) + 1
        } else {
            0
        };
        file.scans = u32::try_from(count_u64).unwrap_or(u32::MAX);
        if let (Some(t0), Some(t1)) = (
            read_f64_le(data, header + 0x48),
            read_f64_le(data, header + 0x50),
        ) {
            file.minutes = Some((t0, t1));
        }
    }

    if let Some(inst_off) = header.checked_add(INSTRUMENT_ID_OFF + 12)
        && inst_off < data.len()
        && let Ok((name, next)) = read_utf16_field(data, inst_off)
        && let Ok((model, _)) = read_utf16_field(data, next)
    {
        file.instrument = if !model.is_empty() {
            Some(format!("Thermo {model}"))
        } else if !name.is_empty() {
            Some(name)
        } else {
            None
        };
    }

    Ok((file, Some(header)))
}

/// Parses the header, sample row, scan count, and instrument model from a Thermo `.raw` buffer.
pub fn parse_run(data: &[u8]) -> Result<ThermoFile> {
    Ok(parse_run_with_offset(data)?.0)
}

/// Parses the header and all sample rows from a Thermo `.sld` sequence file.
pub fn parse_sequence(data: &[u8]) -> Result<ThermoFile> {
    let mut file = parse_header(data)?;
    let (_, row_span_list) = sequence_spans(data)?;
    file.rows = row_span_list
        .iter()
        .map(|spans| decode_row(data, spans))
        .collect();
    Ok(file)
}

/// Returns physical byte spans for the four header login slots and every UTF-16LE string field.
pub fn text_spans(data: &[u8], is_sequence: bool) -> Result<ThermoTextSpans> {
    if !data.starts_with(MAGIC) {
        return Err(InstrumentError::InvalidThermo(
            "not a Thermo Finnigan file".to_owned(),
        ));
    }
    let mut logins = Vec::with_capacity(LOGIN_OFFSETS.len());
    for off in LOGIN_OFFSETS {
        let slice = data.get(off..off + LOGIN_BYTES).ok_or_else(truncated)?;
        let units: Vec<u16> = slice
            .as_chunks::<2>()
            .0
            .iter()
            .map(|&c| u16::from_le_bytes(c))
            .collect();
        String::from_utf16(&units).map_err(|_| {
            InstrumentError::InvalidThermo(format!("invalid UTF-16LE in login slot at byte {off}"))
        })?;
        let nul_units = units.iter().position(|&u| u == 0).unwrap_or(units.len());
        logins.push((off, off + nul_units * 2, off + LOGIN_BYTES));
    }
    let strings = if is_sequence {
        let (labels, rows) = sequence_spans(data)?;
        let mut all = labels;
        for row in rows {
            all.extend(row);
        }
        all
    } else {
        row_spans(data, HEADER_LEN)?.0
    };
    Ok(ThermoTextSpans { logins, strings })
}

/// Returns the `(offset, length)` byte ranges to read scan bounds and index pointers from an MS run header.
pub fn ms_header_ranges(header_offset: u64) -> Result<((u64, usize), (u64, usize))> {
    let bounds_off = header_offset.checked_add(8).ok_or_else(truncated)?;
    let index_off = header_offset
        .checked_add(SCAN_INDEX_OFF as u64)
        .ok_or_else(truncated)?;
    Ok(((bounds_off, 8), (index_off, 16)))
}

/// Decodes the 8-byte scan bounds and 16-byte index/data pointers into `(index_offset, index_bytes, data_offset)`.
pub fn decode_scan_range(bounds: &[u8], pointers: &[u8]) -> Result<(u64, usize, u64)> {
    let first = read_u32_le(bounds, 0).ok_or_else(truncated)?;
    let last = read_u32_le(bounds, 4).ok_or_else(truncated)?;
    let index_offset = read_u64_le(pointers, 0).ok_or_else(truncated)?;
    let data_offset = read_u64_le(pointers, 8).ok_or_else(truncated)?;
    let count = if last != 0 && last >= first {
        u64::from(last) - u64::from(first) + 1
    } else {
        0
    };
    let count_usize = usize::try_from(count).map_err(|_| truncated())?;
    let index_bytes = count_usize
        .checked_mul(SCAN_ENTRY_LEN)
        .ok_or_else(truncated)?;
    Ok((index_offset, index_bytes, data_offset))
}

type SequenceSpans = (Vec<(usize, usize)>, Vec<Vec<(usize, usize)>>);

fn sequence_spans(data: &[u8]) -> Result<SequenceSpans> {
    let mut labels = Vec::with_capacity(5);
    let mut at = SEQUENCE_LABELS_OFF;
    for _ in 0..5 {
        let span = string_span(data, at)?;
        at = span.1;
        labels.push(span);
    }
    let count_off = at.checked_add(64).ok_or_else(truncated)?;
    let count = read_u32_le(data, count_off).ok_or_else(truncated)? as usize;
    at = at.checked_add(68).ok_or_else(truncated)?;

    let max_rows = data.len().saturating_sub(at) / MIN_ROW_BYTES;
    if count > max_rows {
        return Err(truncated());
    }

    let mut rows = Vec::with_capacity(count);
    for _ in 0..count {
        let (spans, next_at) = row_spans(data, at)?;
        rows.push(spans);
        at = next_at;
    }
    Ok((labels, rows))
}

/// Returns the mass spectrometer run header file offset (controller type 0), if present.
pub fn ms_header_offset(data: &[u8]) -> Result<Option<u64>> {
    parse_header(data)?;
    let (_, header) = parse_ms_header(data)?;
    Ok(header)
}

/// Decodes raw 88-byte scan index entries into retention times (min), TIC values, and peak offsets.
pub fn decode_scan_index(entries: &[u8], data_offset: u64) -> Result<ThermoScans> {
    if !entries.len().is_multiple_of(SCAN_ENTRY_LEN) {
        return Err(InstrumentError::InvalidThermo(
            "scan index byte length is not a multiple of 88".to_owned(),
        ));
    }
    let count = entries.len() / SCAN_ENTRY_LEN;
    let mut minutes = Vec::with_capacity(count);
    let mut tic = Vec::with_capacity(count);
    let mut offsets = Vec::with_capacity(count);

    for i in 0..count {
        let base = i * SCAN_ENTRY_LEN;
        minutes.push(read_f64_le(entries, base + 24).ok_or_else(truncated)?);
        tic.push(read_f64_le(entries, base + 32).ok_or_else(truncated)?);
        offsets.push(read_u64_le(entries, base + 72).ok_or_else(truncated)?);
    }

    Ok(ThermoScans {
        minutes,
        tic,
        offsets,
        data_offset,
    })
}

/// Decodes a 56-byte scan peak record into `(relative_start, byte_len)`.
pub fn decode_peak_record(record: &[u8]) -> Result<(u32, u32)> {
    if record.len() < PEAK_RECORD_LEN {
        return Err(truncated());
    }
    let start = read_u32_le(record, 20).ok_or_else(truncated)?;
    let size = read_u32_le(record, 48).ok_or_else(truncated)?;
    Ok((start, size))
}

/// Decodes packed little-endian `f32` `(m/z, intensity)` peak pairs.
pub fn decode_peak_pairs(bytes: &[u8]) -> Result<(Vec<f64>, Vec<f64>)> {
    if !bytes.len().is_multiple_of(8) {
        return Err(InstrumentError::InvalidThermo(
            "peak pairs byte length is not a multiple of 8".to_owned(),
        ));
    }
    let n = bytes.len() / 8;
    let mut mz = Vec::with_capacity(n);
    let mut intensity = Vec::with_capacity(n);
    for i in 0..n {
        let base = i * 8;
        mz.push(f64::from(read_f32_le(bytes, base).ok_or_else(truncated)?));
        intensity.push(f64::from(
            read_f32_le(bytes, base + 4).ok_or_else(truncated)?,
        ));
    }
    Ok((mz, intensity))
}

/// Parses a complete Thermo Finnigan `.raw` file into a unified [`InstrumentReading`].
pub fn read_raw(data: &[u8]) -> Result<InstrumentReading> {
    let (run, header_opt) = parse_run_with_offset(data)?;
    let mut traces = Vec::new();
    if let Some(header) = header_opt
        && run.scans > 0
        && let Some(idx_field) = header.checked_add(SCAN_INDEX_OFF)
        && let Some(data_field) = idx_field.checked_add(8)
        && let (Some(idx_off), Some(data_off)) =
            (read_u64_le(data, idx_field), read_u64_le(data, data_field))
        && let Ok(start) = usize::try_from(idx_off)
        && let Some(byte_len) = (run.scans as usize).checked_mul(SCAN_ENTRY_LEN)
        && let Some(end) = start.checked_add(byte_len)
        && let Some(slice) = data.get(start..end)
        && let Ok(scans) = decode_scan_index(slice, data_off)
    {
        let points = scans.minutes.into_iter().zip(scans.tic).collect();
        traces.push(Trace {
            label: "TIC".to_owned(),
            x_unit: "min".to_owned(),
            y_unit: "counts".to_owned(),
            points,
        });
    }

    let first_row = run.rows.first();
    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::ThermoRaw,
        acquired_at: Some(run.second),
        sample_name: first_row
            .map(|r| r.sample_name.clone())
            .filter(|s| !s.is_empty()),
        instrument: run.instrument,
        method: InstrumentMethod {
            name: first_row
                .map(|r| r.method.clone())
                .filter(|s| !s.is_empty()),
            ..InstrumentMethod::default()
        },
        traces,
        ..InstrumentReading::default()
    })
}

fn parse_header(data: &[u8]) -> Result<ThermoFile> {
    if !data.starts_with(MAGIC) {
        return Err(InstrumentError::InvalidThermo(
            "not a Thermo Finnigan file".to_owned(),
        ));
    }
    let version = read_u32_le(data, 0x24).ok_or_else(truncated)?;
    if version != VERSION {
        return Err(InstrumentError::InvalidThermo(format!(
            "Finnigan file version {version}, not {VERSION}"
        )));
    }

    let mut logins = BTreeSet::new();
    for off in LOGIN_OFFSETS {
        let slice = data.get(off..off + LOGIN_BYTES).ok_or_else(truncated)?;
        let decoded = decode_utf16le(slice);
        let first_nul = decoded.split('\0').next().unwrap_or("");
        logins.insert(first_nul.to_owned());
    }

    Ok(ThermoFile {
        first: parse_filetime_plus_0800(data, 0x28)?,
        second: parse_filetime_plus_0800(data, 0x98)?,
        logins,
        rows: Vec::new(),
        instrument: None,
        scans: 0,
        minutes: None,
    })
}

fn parse_ms_header(data: &[u8]) -> Result<(ThermoRow, Option<u64>)> {
    let (spans, at) = row_spans(data, HEADER_LEN)?;
    let row = decode_row(data, &spans);
    let (_, at) = read_utf16_field(data, at.checked_add(24).ok_or_else(truncated)?)?;
    let count_off = at.checked_add(CONTROLLERS_OFF).ok_or_else(truncated)?;
    let n_controllers = read_u32_le(data, count_off).ok_or_else(truncated)? as usize;
    let table_start = at.checked_add(CONTROLLER_TABLE_OFF).ok_or_else(truncated)?;
    let max_entries = data.len().saturating_sub(table_start) / 16;
    if n_controllers > max_entries {
        return Err(truncated());
    }
    for i in 0..n_controllers {
        let entry = table_start + 16 * i;
        let kind = read_u32_le(data, entry).ok_or_else(truncated)?;
        let header = read_u64_le(data, entry + 8).ok_or_else(truncated)?;
        if kind == 0 {
            return Ok((row, Some(header)));
        }
    }
    Ok((row, None))
}

fn string_span(data: &[u8], at: usize) -> Result<(usize, usize)> {
    let chars = read_u32_le(data, at).ok_or_else(truncated)? as usize;
    let byte_len = chars.checked_mul(2).ok_or_else(truncated)?;
    let start = at.checked_add(4).ok_or_else(truncated)?;
    let end = start.checked_add(byte_len).ok_or_else(truncated)?;
    if end > data.len() {
        return Err(InstrumentError::InvalidThermo(format!(
            "string at byte {at} runs past the end of the file"
        )));
    }
    Ok((start, end))
}

fn read_utf16_field(data: &[u8], at: usize) -> Result<(String, usize)> {
    let (start, end) = string_span(data, at)?;
    Ok((decode_utf16le(&data[start..end]), end))
}

fn row_spans(data: &[u8], start: usize) -> Result<(Vec<(usize, usize)>, usize)> {
    let mut at = start.checked_add(INJECTION_LEN).ok_or_else(truncated)?;
    let mut spans = Vec::with_capacity(31);
    for i in 0..31 {
        let field_at = if i == 16 {
            at.checked_add(4).ok_or_else(truncated)?
        } else {
            at
        };
        let span = string_span(data, field_at)?;
        at = span.1;
        spans.push(span);
    }
    Ok((spans, at))
}

fn decode_row(data: &[u8], spans: &[(usize, usize)]) -> ThermoRow {
    let field = |idx: usize| -> String {
        spans
            .get(idx)
            .and_then(|&(s, e)| data.get(s..e))
            .map(decode_utf16le)
            .unwrap_or_default()
    };
    ThermoRow {
        sample_name: field(1),
        sample_id: field(2),
        method: file_name(&field(9)).to_owned(),
        data_file: file_name(&field(11)).to_owned(),
        vial: field(13),
    }
}

fn parse_filetime_plus_0800(data: &[u8], at: usize) -> Result<String> {
    let ticks = read_u64_le(data, at).ok_or_else(truncated)?;
    let secs = (ticks / 10_000_000) as i64 - FILETIME_EPOCH_DIFF_SECS;
    Ok(format_unix_timestamp(secs, 0, 'T', "+08:00"))
}

fn truncated() -> InstrumentError {
    InstrumentError::InvalidThermo("Finnigan file truncated".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn utf16_field(text: &str) -> Vec<u8> {
        let units: Vec<u16> = text.encode_utf16().collect();
        let mut out = Vec::with_capacity(4 + units.len() * 2);
        out.extend_from_slice(&(units.len() as u32).to_le_bytes());
        for u in units {
            out.extend_from_slice(&u.to_le_bytes());
        }
        out
    }

    fn synthetic_header(login: &str) -> Vec<u8> {
        let mut buf = vec![0u8; HEADER_LEN];
        buf[..MAGIC.len()].copy_from_slice(MAGIC);
        buf[0x24..0x28].copy_from_slice(&VERSION.to_le_bytes());
        let login_bytes: Vec<u8> = login.encode_utf16().flat_map(u16::to_le_bytes).collect();
        buf[0x30..0x30 + login_bytes.len()].copy_from_slice(&login_bytes);
        buf
    }

    fn synthetic_row(sample: &str) -> Vec<u8> {
        let mut out = vec![0u8; INJECTION_LEN];
        for i in 0..31 {
            if i == 16 {
                out.extend_from_slice(&0u32.to_le_bytes());
            }
            let text = if i == 1 { sample } else { "" };
            out.extend_from_slice(&utf16_field(text));
        }
        out
    }

    #[test]
    fn text_spans_finds_logins_and_row_strings_in_run_and_sequence() {
        let mut run = synthetic_header("LAB\\alice");
        run.extend_from_slice(&synthetic_row("TBP-1"));
        let spans = text_spans(&run, false).expect("run spans");
        assert_eq!(spans.logins.len(), 4);
        let (s, e, slot_end) = spans.logins[0];
        assert_eq!(decode_utf16le(&run[s..e]), "LAB\\alice");
        assert_eq!(slot_end - s, LOGIN_BYTES);
        assert_eq!(spans.strings.len(), 31);
        let (rs, re) = spans.strings[1];
        assert_eq!(decode_utf16le(&run[rs..re]), "TBP-1");

        let mut seq = synthetic_header("LAB\\bob");
        seq.resize(SEQUENCE_LABELS_OFF, 0);
        for label in ["L1", "L2", "L3", "L4", "L5"] {
            seq.extend_from_slice(&utf16_field(label));
        }
        seq.extend_from_slice(&[0u8; 64]);
        seq.extend_from_slice(&1u32.to_le_bytes());
        seq.extend_from_slice(&synthetic_row("SEQ-1"));
        let seq_spans = text_spans(&seq, true).expect("sequence spans");
        assert_eq!(seq_spans.strings.len(), 5 + 31);
        assert_eq!(
            decode_utf16le(&seq[seq_spans.strings[0].0..seq_spans.strings[0].1]),
            "L1"
        );
        assert_eq!(
            decode_utf16le(&seq[seq_spans.strings[6].0..seq_spans.strings[6].1]),
            "SEQ-1"
        );
    }

    #[test]
    fn ms_header_ranges_and_decode_scan_range_match_layout() {
        let ((b_off, b_len), (p_off, p_len)) = ms_header_ranges(1000).expect("ranges");
        assert_eq!((b_off, b_len), (1008, 8));
        assert_eq!((p_off, p_len), (1000 + SCAN_INDEX_OFF as u64, 16));

        let mut bounds = Vec::new();
        bounds.extend_from_slice(&1u32.to_le_bytes());
        bounds.extend_from_slice(&10u32.to_le_bytes());
        let mut pointers = Vec::new();
        pointers.extend_from_slice(&4096u64.to_le_bytes());
        pointers.extend_from_slice(&8192u64.to_le_bytes());
        assert_eq!(
            decode_scan_range(&bounds, &pointers).expect("scan range"),
            (4096, 10 * SCAN_ENTRY_LEN, 8192)
        );
    }
}
