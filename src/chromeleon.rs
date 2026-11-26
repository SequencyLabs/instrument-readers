//! Pure-Rust reader for Thermo Chromeleon 7 `SignHdr` `.raw` signal streams.
//!
//! Decodes `PtsLL2Df` second-difference varint blocks used by Chromeleon UV-Vis and detector channels.

use crate::binary::{read_f64_le, read_u64_le};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentReading, Technique, Trace, VendorFormat};

const SIGNAL_TAG: &[u8; 8] = b"SignHdr\0";
const POINTS_TAG: &[u8; 8] = b"PtsLL2Df";

/// Extracts the channel name (such as `"UV_VIS_1"`) from a Chromeleon `SignHdr\0` buffer.
#[must_use]
pub fn signal_name(data: &[u8]) -> Option<String> {
    if !data.starts_with(SIGNAL_TAG) {
        return None;
    }
    let body = data.get(16..)?;
    let name_len = *body.get(14)? as usize;
    let raw = body.get(15..15usize.checked_add(name_len)?)?;
    let name: String = raw.iter().map(|&b| char::from(b)).collect();
    let trimmed = name.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_owned())
}

/// Decodes a single Chromeleon `SignHdr\0` signal file into `(minutes, values)`.
pub fn decode_signal(data: &[u8]) -> Result<(Vec<f64>, Vec<f64>)> {
    if !data.starts_with(SIGNAL_TAG) {
        return Err(InstrumentError::InvalidChromeleon(
            "missing SignHdr magic".to_owned(),
        ));
    }

    let mut minutes = Vec::new();
    let mut values = Vec::new();
    let mut at = 0usize;

    while at < data.len() {
        let tag = data
            .get(at..at.checked_add(8).ok_or_else(truncated)?)
            .ok_or_else(truncated)?;
        let size_u64 = read_u64_le(data, at + 8).ok_or_else(truncated)?;
        let size = usize::try_from(size_u64).map_err(|_| truncated())?;
        let step = size.max(16);
        let body_end = at.checked_add(size).ok_or_else(truncated)?;
        let next_at = at.checked_add(step).ok_or_else(truncated)?;
        if body_end > data.len() || body_end < at + 16 {
            return Err(truncated());
        }
        let body = &data[at + 16..body_end];
        at = next_at;

        if tag == SIGNAL_TAG {
            continue;
        }
        if tag != POINTS_TAG {
            let tag_str = String::from_utf8_lossy(tag);
            return Err(InstrumentError::InvalidChromeleon(format!(
                "unsupported block type {tag_str}"
            )));
        }

        let per_minute = read_f64_le(body, 0).ok_or_else(truncated)?;
        let per_unit = read_f64_le(body, 8).ok_or_else(truncated)?;
        if per_minute == 0.0 || !per_minute.is_finite() || per_unit == 0.0 || !per_unit.is_finite()
        {
            return Err(InstrumentError::InvalidChromeleon(
                "invalid PtsLL2Df scale factor".to_owned(),
            ));
        }

        let pairs = decode_varints(&body[16..])?;
        if pairs.is_empty() || !pairs.len().is_multiple_of(2) {
            return Err(InstrumentError::InvalidChromeleon(
                "PtsLL2Df block has unpaired times and values".to_owned(),
            ));
        }

        let (pair_chunks, _) = pairs.as_chunks::<2>();
        if let Some((&[t0, v0], rest)) = pair_chunks.split_first() {
            let mut t_acc = t0;
            let mut v_acc = v0;
            let mut t_step = 0i64;
            let mut v_step = 0i64;
            minutes.push((t_acc as f64) / per_minute);
            values.push((v_acc as f64) / per_unit);
            for &[dt, dv] in rest {
                t_step = t_step.wrapping_add(dt);
                t_acc = t_acc.wrapping_add(t_step);
                v_step = v_step.wrapping_add(dv);
                v_acc = v_acc.wrapping_add(v_step);
                minutes.push((t_acc as f64) / per_minute);
                values.push((v_acc as f64) / per_unit);
            }
        }
    }

    if minutes.is_empty() {
        return Err(InstrumentError::InvalidChromeleon(
            "signal contains no points".to_owned(),
        ));
    }

    Ok((minutes, values))
}

/// Parses a single Chromeleon `SignHdr\0` `.raw` file into an [`InstrumentReading`].
pub fn read_signal(data: &[u8]) -> Result<InstrumentReading> {
    let label = signal_name(data).unwrap_or_else(|| "Signal".to_owned());
    let (minutes, values) = decode_signal(data)?;
    let y_unit = unit_for_channel(&label).to_owned();
    let points = minutes.into_iter().zip(values).collect();

    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::ChromeleonRaw,
        traces: vec![Trace {
            label,
            x_unit: "min".to_owned(),
            y_unit,
            points,
        }],
        ..InstrumentReading::default()
    })
}

/// Parses a bundle of Chromeleon `.raw` channel files, skipping unsupported auxiliary blocks
/// (such as `PtsLDiff` pump pressure or oven temperature streams) when readable UV channels exist.
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let mut traces = Vec::new();
    let mut last_err: Option<InstrumentError> = None;

    for &(_, bytes) in files {
        if !bytes.starts_with(SIGNAL_TAG) {
            continue;
        }
        let label = signal_name(bytes).unwrap_or_else(|| "Signal".to_owned());
        match decode_signal(bytes) {
            Ok((minutes, values)) => {
                traces.push(Trace {
                    y_unit: unit_for_channel(&label).to_owned(),
                    label,
                    x_unit: "min".to_owned(),
                    points: minutes.into_iter().zip(values).collect(),
                });
            }
            Err(err) => {
                last_err = Some(err);
            }
        }
    }

    if traces.is_empty() {
        return Err(last_err.unwrap_or_else(|| {
            InstrumentError::InvalidChromeleon("bundle contains no Chromeleon signals".to_owned())
        }));
    }

    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::ChromeleonRaw,
        traces,
        ..InstrumentReading::default()
    })
}

fn unit_for_channel(label: &str) -> &'static str {
    let upper = label.to_ascii_uppercase();
    if upper.starts_with("UV") || upper.starts_with("DAD") || upper.starts_with("VWD") {
        "mAU"
    } else if upper.contains("TEMP") {
        "C"
    } else if upper.contains("PRESSURE") {
        "bar"
    } else {
        "intensity"
    }
}

fn decode_varints(body: &[u8]) -> Result<Vec<i64>> {
    let mut out = Vec::new();
    let mut i = 0usize;

    while i < body.len() {
        let mut byte = body[i];
        let negative = (byte & 0x40) != 0;
        let mut value = u64::from(byte & 0x3F);
        let mut shift = 6u32;
        i += 1;

        while (byte & 0x80) != 0 {
            if i >= body.len() || shift >= 64 {
                return Err(truncated());
            }
            byte = body[i];
            if shift == 62 && (byte & 0x7F) > 1 {
                return Err(InstrumentError::InvalidChromeleon(
                    "varint overflow in PtsLL2Df block".to_owned(),
                ));
            }
            value |= u64::from(byte & 0x7F) << shift;
            shift += 7;
            i += 1;
        }

        if negative && value == 0 {
            return Ok(out);
        }

        let signed = i64::try_from(value).map_err(|_| {
            InstrumentError::InvalidChromeleon("varint overflow in PtsLL2Df block".to_owned())
        })?;
        out.push(if negative { -signed } else { signed });
    }

    Err(InstrumentError::InvalidChromeleon(
        "PtsLL2Df block missing terminating negative zero".to_owned(),
    ))
}

fn truncated() -> InstrumentError {
    InstrumentError::InvalidChromeleon("signal cut short".to_owned())
}
