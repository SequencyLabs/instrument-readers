//! Binary `_FUNC*.IDX` and `_FUNC*.DAT` decoders for Waters MassLynx `.raw` acquisitions.

use crate::binary::{read_f32_le, read_i16_le, read_u16_le, read_u32_le, read_u64_le};
use crate::error::{InstrumentError, Result};

use super::{IDX_RECORD_LEN, MAX_MATRIX_CELLS};

#[derive(Clone, Copy)]
pub(super) enum DatFormat {
    Dat2,
    Dat4,
    Dat6,
    Dat8,
}

impl DatFormat {
    pub(super) const fn from_bytes_per_pair(n: usize) -> Option<Self> {
        match n {
            2 => Some(Self::Dat2),
            4 => Some(Self::Dat4),
            6 => Some(Self::Dat6),
            8 => Some(Self::Dat8),
            _ => None,
        }
    }

    pub(super) const fn record_len(self) -> usize {
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

pub(super) struct IdxScan {
    pub(super) offset: usize,
    pub(super) pair_count: usize,
    pub(super) time_min: f64,
}

pub(super) struct DecodedMatrix {
    pub(super) rt: Vec<f64>,
    pub(super) ylabels: Vec<f64>,
    pub(super) matrix: Vec<f64>,
}

pub(super) fn parse_idx(idx_bytes: &[u8], dat_len: usize) -> Option<(Vec<IdxScan>, DatFormat)> {
    let num_scans = idx_bytes.len() / IDX_RECORD_LEN;
    if num_scans == 0 {
        return None;
    }

    let mut scans = Vec::with_capacity(num_scans);
    let mut nonzero: Vec<(usize, usize)> = Vec::with_capacity(num_scans);

    for i in 0..num_scans {
        let base = i * IDX_RECORD_LEN;
        let offset = read_u32_le(idx_bytes, base)? as usize;
        let pair_count = (read_u32_le(idx_bytes, base + 4)? & 0x3F_FFFF) as usize;
        let time_min = f64::from(read_f32_le(idx_bytes, base + 12)?);
        if pair_count > 0 {
            nonzero.push((offset, pair_count));
        }
        scans.push(IdxScan {
            offset,
            pair_count,
            time_min,
        });
    }

    let &(final_offset, final_count) = nonzero.last()?;
    let tail_bytes = dat_len.checked_sub(final_offset)?;
    if tail_bytes == 0 || !tail_bytes.is_multiple_of(final_count) {
        return None;
    }
    let tail_width = tail_bytes / final_count;
    let fmt = DatFormat::from_bytes_per_pair(tail_width)?;

    let mut inferred_width: Option<usize> = None;
    for pair in nonzero.windows(2) {
        let (off, count) = pair[0];
        let (next_off, _) = pair[1];
        let span = next_off.checked_sub(off)?;
        if span == 0 || !span.is_multiple_of(count) {
            return None;
        }
        let width = span / count;
        DatFormat::from_bytes_per_pair(width)?;
        if let Some(prev) = inferred_width
            && prev != width
        {
            return None;
        }
        inferred_width = Some(width);
    }
    if let Some(width) = inferred_width
        && width != tail_width
    {
        return None;
    }
    if !nonzero.iter().all(|&(off, count)| {
        count
            .checked_mul(tail_width)
            .and_then(|len| off.checked_add(len))
            .is_some_and(|end| end <= dat_len)
    }) {
        return None;
    }

    Some((scans, fmt))
}

pub(super) fn decode_function_trace(
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

pub(super) fn decode_function_matrix(
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
            let scaled = (key / bin_width).round_ties_even();
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
