//! Pure-Rust reader for Waters MassLynx `.raw` LC-MS directories.
//!
//! Decodes `_HEADER.TXT` metadata, `_extern.inf` MS polarities, `_FUNCTNS.INF` channel labels,
//! `_FUNC*.IDX` scan indexes, and `_FUNC*.DAT` 2-, 4-, 6-, and 8-byte intensity records into
//! both 1D TIC/UV traces and 2D `(rt, mz, intensity)` matrices.

use std::collections::BTreeMap;

use crate::binary::{file_name, read_f32_le, read_i16_le, read_u16_le, read_u32_le, read_u64_le};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentReading, Technique, Trace, VendorFormat};

const IDX_RECORD_LEN: usize = 22;
const FUNCTNS_RECORD_LEN: usize = 416;
const FUNCTNS_MZ_OFFSET: usize = 160;
const FUNCTNS_MZ_SLOTS: usize = 32;
const MAX_MATRIX_CELLS: usize = 50_000_000;

/// One UV channel extracted from a Waters MassLynx `.raw` run.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersUvChannel {
    /// Wavelength or channel label from `_FUNCTNS.INF`.
    pub wavelength: f64,
    /// Absorbance values aligned to `uv_rt`.
    pub values: Vec<f64>,
}

/// Full 2D MS + UV polarity slice parsed from a Waters MassLynx `.raw` directory.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersPolarityRun {
    /// Ion polarity (`"+"` or `"-"`).
    pub polarity: String,
    /// Acquisition date and time string from `_HEADER.TXT`.
    pub acquired: Option<String>,
    /// Instrument model from `_HEADER.TXT`.
    pub instrument: Option<String>,
    /// Inlet method filename from `_HEADER.TXT`.
    pub inlet_method: Option<String>,
    /// MS method filename from `_HEADER.TXT`.
    pub ms_method: Option<String>,
    /// Retention times in minutes for the MS scans.
    pub rt: Vec<f64>,
    /// Binned m/z axis values.
    pub mz: Vec<f64>,
    /// Row-major intensity matrix (`rt.len() * mz.len()`).
    pub intensity: Vec<f64>,
    /// Retention times in minutes for the UV channels, if present.
    pub uv_rt: Option<Vec<f64>>,
    /// UV channels aligned to `uv_rt`.
    pub uv: Vec<WatersUvChannel>,
}

#[derive(Clone, Copy)]
enum DatFormat {
    Dat2,
    Dat4,
    Dat6,
    Dat8,
}

impl DatFormat {
    const fn from_bytes_per_pair(n: usize) -> Option<Self> {
        match n {
            2 => Some(Self::Dat2),
            4 => Some(Self::Dat4),
            6 => Some(Self::Dat6),
            8 => Some(Self::Dat8),
            _ => None,
        }
    }

    const fn record_len(self) -> usize {
        match self {
            Self::Dat2 => 2,
            Self::Dat4 => 4,
            Self::Dat6 => 6,
            Self::Dat8 => 8,
        }
    }

    fn sum_scan(self, slice: &[u8]) -> f64 {
        let step = self.record_len();
        let mut total = 0.0;
        let mut pos = 0;
        while pos + step <= slice.len() {
            total += match self {
                Self::Dat2 => decode_dat2_value(read_u16_le(slice, pos).unwrap_or(0)),
                Self::Dat4 => decode_dat4_value(read_u32_le(slice, pos).unwrap_or(0)),
                Self::Dat6 => decode_dat6_pair(slice, pos).1,
                Self::Dat8 => decode_dat8_pair(slice, pos).1,
            };
            pos += step;
        }
        total
    }
}

struct WatersBundle<'a> {
    meta: HeaderMeta,
    polarities: Vec<char>,
    functns_bytes: Option<&'a [u8]>,
    func_dats: BTreeMap<u32, &'a [u8]>,
    func_idxs: BTreeMap<u32, &'a [u8]>,
}

struct HeaderMeta {
    sample_name: Option<String>,
    acquired_at: Option<String>,
    instrument: Option<String>,
    inlet_method: Option<String>,
    ms_method: Option<String>,
    calibrations: Vec<Vec<f64>>,
}

struct IdxScan {
    offset: usize,
    pair_count: usize,
    time_min: f64,
}

struct DecodedMatrix {
    rt: Vec<f64>,
    ylabels: Vec<f64>,
    matrix: Vec<f64>,
}

impl HeaderMeta {
    const fn empty() -> Self {
        Self {
            sample_name: None,
            acquired_at: None,
            instrument: None,
            inlet_method: None,
            ms_method: None,
            calibrations: Vec::new(),
        }
    }
}

fn collect_bundle<'a>(files: &[(&str, &'a [u8])]) -> Result<WatersBundle<'a>> {
    let mut header_bytes: Option<&[u8]> = None;
    let mut extern_bytes: Option<&[u8]> = None;
    let mut functns_bytes: Option<&[u8]> = None;
    let mut func_dats: BTreeMap<u32, &[u8]> = BTreeMap::new();
    let mut func_idxs: BTreeMap<u32, &[u8]> = BTreeMap::new();

    for &(path, bytes) in files {
        let upper = file_name(path).to_ascii_uppercase();
        if upper == "_HEADER.TXT" {
            header_bytes = Some(bytes);
        } else if upper == "_EXTERN.INF" {
            extern_bytes = Some(bytes);
        } else if upper == "_FUNCTNS.INF" {
            functns_bytes = Some(bytes);
        } else if let Some(num) = parse_func_number(&upper, ".DAT") {
            func_dats.insert(num, bytes);
        } else if let Some(num) = parse_func_number(&upper, ".IDX") {
            func_idxs.insert(num, bytes);
        }
    }

    let meta = header_bytes
        .map(parse_header_txt)
        .transpose()?
        .unwrap_or_else(HeaderMeta::empty);
    let polarities = extern_bytes
        .map(parse_extern_polarities)
        .unwrap_or_default();

    Ok(WatersBundle {
        meta,
        polarities,
        functns_bytes,
        func_dats,
        func_idxs,
    })
}

/// Parses a Waters MassLynx `_HEADER.TXT` file into an [`InstrumentReading`] with metadata only.
pub fn read_header(data: &[u8]) -> Result<InstrumentReading> {
    let meta = parse_header_txt(data)?;
    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::WatersRaw,
        acquired_at: meta.acquired_at,
        sample_name: meta.sample_name,
        instrument: meta.instrument,
        ..InstrumentReading::default()
    })
}

/// Parses a Waters MassLynx `.raw` directory bundle into an [`InstrumentReading`].
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let bundle = collect_bundle(files)?;
    let mut traces = Vec::new();

    for (func_idx, (&func_num, &dat_bytes)) in bundle.func_dats.iter().enumerate() {
        let Some(&idx_bytes) = bundle.func_idxs.get(&func_num) else {
            continue;
        };
        if dat_bytes.is_empty() || idx_bytes.len() < IDX_RECORD_LEN {
            continue;
        }
        let Some((scans, fmt)) = parse_idx(idx_bytes, dat_bytes.len()) else {
            continue;
        };
        let Some(points) = decode_function_trace(dat_bytes, &scans, fmt) else {
            continue;
        };
        if points.is_empty() {
            continue;
        }

        let polarity = bundle.polarities.get(func_idx).copied();
        let (label, y_unit) = match polarity {
            Some('+') => (format!("ESI+ TIC (Func {func_num})"), "counts"),
            Some('-') => (format!("ESI- TIC (Func {func_num})"), "counts"),
            _ => (format!("UV (Func {func_num})"), "AU"),
        };
        traces.push(Trace {
            label,
            x_unit: "min".to_owned(),
            y_unit: y_unit.to_owned(),
            points,
        });
    }

    if traces.is_empty() && bundle.meta.sample_name.is_none() && bundle.meta.instrument.is_none() {
        return Err(InstrumentError::InvalidWaters(
            "bundle contains no readable _HEADER.TXT or _FUNC*.DAT/IDX files".to_owned(),
        ));
    }

    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::WatersRaw,
        acquired_at: bundle.meta.acquired_at,
        sample_name: bundle.meta.sample_name,
        instrument: bundle.meta.instrument,
        traces,
        ..InstrumentReading::default()
    })
}

/// Parses every MS polarity (`"+"` and `"-"`) and UV channel from a Waters MassLynx `.raw` bundle.
pub fn read_polarities(
    files: &[(&str, &[u8])],
    bin_width: f64,
) -> Result<BTreeMap<String, WatersPolarityRun>> {
    if !bin_width.is_finite() || bin_width <= 0.0 {
        return Err(InstrumentError::InvalidWaters(
            "bin_width must be positive and finite".to_owned(),
        ));
    }

    let bundle = collect_bundle(files)?;
    if let Some(inf) = bundle.functns_bytes {
        let expected = FUNCTNS_RECORD_LEN
            .checked_mul(bundle.func_dats.len())
            .ok_or_else(|| InstrumentError::InvalidWaters("function count overflow".to_owned()))?;
        if inf.len() != expected {
            return Err(InstrumentError::InvalidWaters(format!(
                "_FUNCTNS.INF size {} does not match {} functions",
                inf.len(),
                bundle.func_dats.len()
            )));
        }
    }

    let func_mzs = bundle
        .functns_bytes
        .map(parse_functns_inf)
        .unwrap_or_default();

    let mut uv_rt: Option<Vec<f64>> = None;
    let mut uv_channels: Vec<WatersUvChannel> = Vec::new();
    let mut ms_by_polarity: BTreeMap<String, DecodedMatrix> = BTreeMap::new();

    for (func_idx, (&func_num, &dat_bytes)) in bundle.func_dats.iter().enumerate() {
        let Some(&idx_bytes) = bundle.func_idxs.get(&func_num) else {
            continue;
        };
        let Some((scans, fmt)) = parse_idx(idx_bytes, dat_bytes.len()) else {
            return Err(InstrumentError::InvalidWaters(format!(
                "unreadable index for _FUNC{func_num:03}.DAT"
            )));
        };

        let polarity = bundle.polarities.get(func_idx).copied();
        let calib = bundle
            .meta
            .calibrations
            .get(func_idx)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let inf_labels = func_num
            .checked_sub(1)
            .and_then(|idx| func_mzs.get(idx as usize))
            .map(Vec::as_slice)
            .unwrap_or(&[]);

        let decoded = decode_function_matrix(
            dat_bytes,
            &scans,
            fmt,
            inf_labels,
            calib,
            if polarity.is_some() { bin_width } else { 1.0 },
            polarity.is_none(),
        )?;

        if let Some(pol) = polarity {
            ms_by_polarity.entry(pol.to_string()).or_insert(decoded);
        } else if uv_rt.is_none() {
            let nrows = decoded.rt.len();
            let ncols = decoded.ylabels.len();
            let mut chans = Vec::with_capacity(ncols);
            for (col, &wl) in decoded.ylabels.iter().enumerate() {
                let mut vals = Vec::with_capacity(nrows);
                for row in 0..nrows {
                    vals.push(decoded.matrix[row * ncols + col]);
                }
                chans.push(WatersUvChannel {
                    wavelength: wl,
                    values: vals,
                });
            }
            uv_rt = Some(decoded.rt);
            uv_channels = chans;
        }
    }

    let mut out = BTreeMap::new();
    for (polarity, ms) in ms_by_polarity {
        out.insert(
            polarity.clone(),
            WatersPolarityRun {
                polarity,
                acquired: bundle.meta.acquired_at.clone(),
                instrument: bundle.meta.instrument.clone(),
                inlet_method: bundle.meta.inlet_method.clone(),
                ms_method: bundle.meta.ms_method.clone(),
                rt: ms.rt,
                mz: ms.ylabels,
                intensity: ms.matrix,
                uv_rt: uv_rt.clone(),
                uv: uv_channels.clone(),
            },
        );
    }

    Ok(out)
}

fn parse_func_number(upper_name: &str, suffix: &str) -> Option<u32> {
    let middle = upper_name.strip_prefix("_FUNC")?.strip_suffix(suffix)?;
    if middle.is_empty() || !middle.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    middle.parse::<u32>().ok()
}

fn parse_header_txt(data: &[u8]) -> Result<HeaderMeta> {
    let text = String::from_utf8_lossy(data);
    let mut sample_name = None;
    let mut acquired_date = None;
    let mut acquired_time = None;
    let mut instrument = None;
    let mut inlet_method = None;
    let mut ms_method = None;
    let mut calibrations = Vec::new();
    let mut saw_marker = false;

    for line in text.lines() {
        let trimmed = line.trim();
        let Some(rest) = trimmed.strip_prefix("$$ ") else {
            continue;
        };
        saw_marker = true;
        let Some((key, val)) = rest.split_once(':') else {
            continue;
        };
        let key = key.trim();
        let value = val.trim();
        if key.starts_with("Cal Function") {
            let coeffs: Vec<f64> = value
                .split(',')
                .filter_map(|p| p.trim().parse::<f64>().ok())
                .collect();
            calibrations.push(coeffs);
            continue;
        }
        if value.is_empty() {
            continue;
        }
        match key {
            "Acquired Name" => sample_name = Some(value.to_owned()),
            "Acquired Date" => acquired_date = Some(value.to_owned()),
            "Acquired Time" => acquired_time = Some(value.to_owned()),
            "Instrument" => instrument = Some(value.to_owned()),
            "Inlet Method" => {
                let name = file_name(value);
                if !name.is_empty() {
                    inlet_method = Some(name.to_owned());
                }
            }
            "MS Method" => {
                let name = file_name(value);
                if !name.is_empty() {
                    ms_method = Some(name.to_owned());
                }
            }
            _ => {}
        }
    }

    if !saw_marker {
        return Err(InstrumentError::InvalidWaters(
            "missing $$ header markers in _HEADER.TXT".to_owned(),
        ));
    }

    let acquired_at = match (acquired_date, acquired_time) {
        (Some(d), Some(t)) => Some(format!("{d} {t}")),
        (Some(d), None) => Some(d),
        _ => None,
    };

    Ok(HeaderMeta {
        sample_name,
        acquired_at,
        instrument,
        inlet_method,
        ms_method,
        calibrations,
    })
}

fn parse_functns_inf(data: &[u8]) -> Vec<Vec<f64>> {
    let num_funcs = data.len() / FUNCTNS_RECORD_LEN;
    let mut out = Vec::with_capacity(num_funcs);
    for f in 0..num_funcs {
        let base = f * FUNCTNS_RECORD_LEN + FUNCTNS_MZ_OFFSET;
        let mut labels = Vec::new();
        for slot in 0..FUNCTNS_MZ_SLOTS {
            if let Some(val) = read_f32_le(data, base + slot * 4)
                && val != 0.0
                && val.is_finite()
            {
                labels.push(f64::from(val));
            }
        }
        out.push(labels);
    }
    out
}

fn parse_extern_polarities(data: &[u8]) -> Vec<char> {
    let text = String::from_utf8_lossy(data);
    let lines: Vec<&str> = text.lines().collect();
    let mut out = Vec::new();

    for (i, line) in lines.iter().enumerate() {
        if line.starts_with("Instrument Parameters") {
            let candidate = lines
                .get(i + 1)
                .filter(|l| l.trim_start().starts_with("Polarity"))
                .or_else(|| {
                    lines
                        .get(i + 2)
                        .filter(|l| l.trim_start().starts_with("Polarity"))
                });
            if let Some(pol_line) = candidate
                && let Some(pol) = extract_polarity_char(pol_line)
            {
                out.push(pol);
            }
        } else if line.starts_with("Function ") {
            for next in lines.iter().skip(i + 1).take(4) {
                if next.trim_start().starts_with("Polarity")
                    && let Some(pol) = extract_polarity_char(next)
                {
                    out.push(pol);
                    break;
                }
            }
        }
    }

    out
}

fn extract_polarity_char(line: &str) -> Option<char> {
    let last = line.trim().chars().last()?;
    matches!(last, '+' | '-').then_some(last)
}

fn parse_idx(idx_bytes: &[u8], dat_len: usize) -> Option<(Vec<IdxScan>, DatFormat)> {
    let num_scans = idx_bytes.len() / IDX_RECORD_LEN;
    if num_scans == 0 {
        return None;
    }

    let mut scans = Vec::with_capacity(num_scans);
    let mut last_nonzero: Option<(usize, usize)> = None;

    for i in 0..num_scans {
        let base = i * IDX_RECORD_LEN;
        let offset = read_u32_le(idx_bytes, base)? as usize;
        let pair_count = (read_u32_le(idx_bytes, base + 4)? & 0x3F_FFFF) as usize;
        let time_min = f64::from(read_f32_le(idx_bytes, base + 12)?);
        if pair_count > 0 {
            last_nonzero = Some((offset, pair_count));
        }
        scans.push(IdxScan {
            offset,
            pair_count,
            time_min,
        });
    }

    let (final_offset, final_count) = last_nonzero?;
    let tail_bytes = dat_len.checked_sub(final_offset)?;
    if tail_bytes == 0 || !tail_bytes.is_multiple_of(final_count) {
        return None;
    }
    let fmt = DatFormat::from_bytes_per_pair(tail_bytes / final_count)?;
    Some((scans, fmt))
}

fn decode_function_trace(
    dat_bytes: &[u8],
    scans: &[IdxScan],
    fmt: DatFormat,
) -> Option<Vec<(f64, f64)>> {
    let mut points = Vec::with_capacity(scans.len());
    let record_len = fmt.record_len();

    for scan in scans {
        if !scan.time_min.is_finite() {
            continue;
        }
        if scan.pair_count == 0 {
            points.push((scan.time_min, 0.0));
            continue;
        }
        let byte_len = scan.pair_count.checked_mul(record_len)?;
        let slice = dat_bytes.get(scan.offset..scan.offset.checked_add(byte_len)?)?;
        points.push((scan.time_min, fmt.sum_scan(slice)));
    }

    Some(points)
}

fn alloc_matrix(nrows: usize, ncols: usize) -> Result<Vec<f64>> {
    let cells = nrows
        .checked_mul(ncols)
        .filter(|&c| c <= MAX_MATRIX_CELLS)
        .ok_or_else(|| {
            InstrumentError::InvalidWaters("decoded matrix exceeds size limit".to_owned())
        })?;
    Ok(vec![0.0; cells])
}

fn decode_function_matrix(
    dat_bytes: &[u8],
    scans: &[IdxScan],
    fmt: DatFormat,
    inf_labels: &[f64],
    calib: &[f64],
    bin_width: f64,
    sparse_keys: bool,
) -> Result<DecodedMatrix> {
    let rt: Vec<f64> = scans.iter().map(|s| s.time_min).collect();
    let record_len = fmt.record_len();

    if matches!(fmt, DatFormat::Dat2 | DatFormat::Dat4) {
        let ncols = inf_labels.len().max(1);
        let ylabels = if inf_labels.is_empty() {
            vec![0.0]
        } else {
            inf_labels.to_vec()
        };
        let mut matrix = alloc_matrix(scans.len(), ncols)?;
        for (row, scan) in scans.iter().enumerate() {
            let byte_len = scan.pair_count.checked_mul(record_len).ok_or_else(|| {
                InstrumentError::InvalidWaters("scan byte length overflow".to_owned())
            })?;
            let end = scan
                .offset
                .checked_add(byte_len)
                .ok_or_else(|| InstrumentError::InvalidWaters("scan offset overflow".to_owned()))?;
            let slice = dat_bytes.get(scan.offset..end).ok_or_else(|| {
                InstrumentError::InvalidWaters("truncated _FUNC*.DAT scan".to_owned())
            })?;
            for col in 0..scan.pair_count.min(ncols) {
                let val = match fmt {
                    DatFormat::Dat2 => decode_dat2_value(read_u16_le(slice, col * 2).unwrap_or(0)),
                    DatFormat::Dat4 => decode_dat4_value(read_u32_le(slice, col * 4).unwrap_or(0)),
                    _ => 0.0,
                };
                matrix[row * ncols + col] = val;
            }
        }
        return Ok(DecodedMatrix {
            rt,
            ylabels,
            matrix,
        });
    }

    let mut scan_bins: Vec<Vec<(i64, f64)>> = Vec::with_capacity(scans.len());
    let mut min_bin: Option<i64> = None;
    let mut max_bin: Option<i64> = None;

    for scan in scans {
        let byte_len = scan.pair_count.checked_mul(record_len).ok_or_else(|| {
            InstrumentError::InvalidWaters("scan byte length overflow".to_owned())
        })?;
        let end = scan
            .offset
            .checked_add(byte_len)
            .ok_or_else(|| InstrumentError::InvalidWaters("scan offset overflow".to_owned()))?;
        let slice = dat_bytes.get(scan.offset..end).ok_or_else(|| {
            InstrumentError::InvalidWaters("truncated _FUNC*.DAT scan".to_owned())
        })?;
        let mut pairs = Vec::with_capacity(scan.pair_count);
        for i in 0..scan.pair_count {
            let (raw_key, val) = match fmt {
                DatFormat::Dat6 => decode_dat6_pair(slice, i * 6),
                DatFormat::Dat8 => decode_dat8_pair(slice, i * 8),
                _ => (0.0, 0.0),
            };
            let key = apply_calibration(raw_key, calib);
            if !key.is_finite() || !val.is_finite() {
                continue;
            }
            let scaled = (key / bin_width).round();
            if !scaled.is_finite() || scaled < (i64::MIN as f64) || scaled > (i64::MAX as f64) {
                continue;
            }
            let bin = scaled as i64;
            min_bin = Some(min_bin.map_or(bin, |m| m.min(bin)));
            max_bin = Some(max_bin.map_or(bin, |m| m.max(bin)));
            pairs.push((bin, val));
        }
        scan_bins.push(pairs);
    }

    let (Some(lo), Some(hi)) = (min_bin, max_bin) else {
        return Ok(DecodedMatrix {
            rt,
            ylabels: Vec::new(),
            matrix: Vec::new(),
        });
    };

    let (ylabels, sparse_lookup): (Vec<f64>, Option<Vec<i64>>) = if sparse_keys {
        let mut unique_bins: Vec<i64> = scan_bins
            .iter()
            .flat_map(|pairs| pairs.iter().map(|&(b, _)| b))
            .collect();
        unique_bins.sort_unstable();
        unique_bins.dedup();
        let labels = unique_bins
            .iter()
            .map(|&b| (b as f64) * bin_width)
            .collect();
        (labels, Some(unique_bins))
    } else {
        let span = hi
            .checked_sub(lo)
            .and_then(|d| d.checked_add(1))
            .and_then(|n| usize::try_from(n).ok())
            .filter(|&n| n <= MAX_MATRIX_CELLS)
            .ok_or_else(|| InstrumentError::InvalidWaters("m/z bin range overflow".to_owned()))?;
        let labels = (0..span)
            .map(|offset| ((lo + offset as i64) as f64) * bin_width)
            .collect();
        (labels, None)
    };

    let ncols = ylabels.len();
    let mut matrix = alloc_matrix(scans.len(), ncols)?;
    for (row, pairs) in scan_bins.iter().enumerate() {
        let row_base = row * ncols;
        for &(bin, val) in pairs {
            let col = if let Some(ref unique) = sparse_lookup {
                unique.binary_search(&bin).ok()
            } else {
                bin.checked_sub(lo)
                    .and_then(|d| usize::try_from(d).ok())
                    .filter(|&c| c < ncols)
            };
            if let Some(c) = col {
                matrix[row_base + c] += val;
            }
        }
    }

    Ok(DecodedMatrix {
        rt,
        ylabels,
        matrix,
    })
}

fn apply_calibration(key: f64, coeffs: &[f64]) -> f64 {
    if coeffs.is_empty() {
        return key;
    }
    let mut acc = 0.0;
    let mut power = 1.0;
    for &c in coeffs {
        acc += c * power;
        power *= key;
    }
    acc
}

fn decode_dat2_value(raw: u16) -> f64 {
    let base = f64::from(raw >> 3);
    let pow = i32::from(raw & 0x7);
    base * 4.0_f64.powi(pow)
}

fn decode_dat4_value(raw: u32) -> f64 {
    let powers = (raw >> 22) as i32;
    let base = f64::from(raw & 0x1F_FFFF) / 1024.0;
    base * 2.0_f64.powi(powers - 10)
}

fn decode_dat6_pair(slice: &[u8], pos: usize) -> (f64, f64) {
    let val_base = read_i16_le(slice, pos).unwrap_or(0);
    let raw = read_u32_le(slice, pos + 2).unwrap_or(0);
    let key_base = f64::from(raw >> 9);
    let key_exp = ((raw & 0x1F0) >> 4) as i32 - 23;
    let key = key_base * 2.0_f64.powi(key_exp);
    let val_pow = (raw & 0xF) as i32;
    let val = f64::from(val_base) * 4.0_f64.powi(val_pow);
    (key, val)
}

fn decode_dat8_pair(slice: &[u8], pos: usize) -> (f64, f64) {
    let raw = read_u64_le(slice, pos).unwrap_or(0);
    let key_bits = raw >> 28;
    let val_bits = raw & 0x0FFF_FFFF;

    let num_key_int = ((key_bits >> 31) as u32).min(31);
    let num_key_frac = 31 - num_key_int;
    let key_int_mask = (1u64 << num_key_int) - 1;
    let key_frac_mask = (1u64 << num_key_frac) - 1;
    let key_int = (key_bits >> num_key_frac) & key_int_mask;
    let key_frac = if num_key_frac == 0 {
        0.0
    } else {
        ((key_bits & key_frac_mask) as f64) / ((1u64 << num_key_frac) as f64)
    };
    let key = (key_int as f64) + key_frac;

    let raw_int_bits = (val_bits >> 22) as u32;
    let (num_int_bits, shift) = if raw_int_bits > 21 {
        (21, raw_int_bits - 21)
    } else {
        (raw_int_bits, 0)
    };
    let num_frac_bits = 21 - num_int_bits;
    let int_mask = (1u64 << num_int_bits) - 1;
    let frac_mask = (1u64 << num_frac_bits) - 1;
    let int_part = ((val_bits >> num_frac_bits) & int_mask) << shift;
    let frac_part = if num_frac_bits == 0 {
        0.0
    } else {
        ((val_bits & frac_mask) as f64) / ((1u64 << num_frac_bits) as f64)
    };
    let val = (int_part as f64) + frac_part;
    (key, val)
}
