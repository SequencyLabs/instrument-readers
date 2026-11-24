//! Binary ChemStation `.CH` chromatogram decoder used inside Agilent OpenLab `.dx` archives.

use crate::binary::{
    decode_utf16le, file_name, read_f32_be, read_f64_be, read_i16_be, read_i32_be,
};
use crate::types::Trace;

pub(crate) struct ParsedCh {
    pub sample_name: Option<String>,
    pub method_name: Option<String>,
    pub instrument: Option<String>,
    pub wavelength_nm: Option<f64>,
    pub trace: Trace,
}

enum ChEncoding {
    Float64Le,
    DoubleDelta,
    SingleDelta,
}

struct ChLayout {
    encoding: ChEncoding,
    float_times: bool,
    utf16: bool,
    data_start: usize,
    scale_off: usize,
    sample_off: usize,
    method_off: usize,
    inst_off: usize,
    unit_off: usize,
    sig_off: Option<usize>,
}

const CH_179: ChLayout = ChLayout {
    encoding: ChEncoding::Float64Le,
    float_times: true,
    utf16: true,
    data_start: 0x1800,
    scale_off: 0x127C,
    sample_off: 0x35A,
    method_off: 0xA0E,
    inst_off: 0xC11,
    unit_off: 0x104C,
    sig_off: Some(0x1075),
};

const CH_181: ChLayout = ChLayout {
    encoding: ChEncoding::DoubleDelta,
    sig_off: None,
    ..CH_179
};

const CH_130: ChLayout = ChLayout {
    encoding: ChEncoding::SingleDelta,
    float_times: false,
    ..CH_179
};

const CH_30: ChLayout = ChLayout {
    encoding: ChEncoding::SingleDelta,
    float_times: false,
    utf16: false,
    data_start: 0x0400,
    scale_off: 0x0284,
    sample_off: 0x018,
    method_off: 0x0E4,
    inst_off: 0x0DA,
    unit_off: 0x0244,
    sig_off: Some(0x0254),
};

pub(crate) fn parse_chemstation_ch(data: &[u8]) -> Option<ParsedCh> {
    let layout = match read_pascal_ascii(data, 0)?.as_str() {
        "179" => &CH_179,
        "181" => &CH_181,
        "130" => &CH_130,
        "30" => &CH_30,
        _ => return None,
    };
    if data.len() <= layout.data_start {
        return None;
    }

    let payload = &data[layout.data_start..];
    let raw_values: Vec<f64> = match layout.encoding {
        ChEncoding::Float64Le => {
            let (chunks, _) = payload.as_chunks::<8>();
            chunks.iter().map(|&b| f64::from_le_bytes(b)).collect()
        }
        ChEncoding::DoubleDelta => decode_ch_double_delta(payload),
        ChEncoding::SingleDelta => decode_ch_delta(payload),
    };
    if raw_values.is_empty() {
        return None;
    }

    let (start_ms, end_ms) = if layout.float_times {
        (
            f64::from(read_f32_be(data, 0x11A)?),
            f64::from(read_f32_be(data, 0x11E)?),
        )
    } else {
        (
            f64::from(read_i32_be(data, 0x11A)?),
            f64::from(read_i32_be(data, 0x11E)?),
        )
    };

    let raw_scale = read_f64_be(data, layout.scale_off)?;
    let scale = if raw_scale == 0.0 || !raw_scale.is_finite() {
        1.0
    } else {
        raw_scale
    };

    let read_str = |off: usize| {
        if layout.utf16 {
            read_pascal_utf16le(data, off)
        } else {
            read_pascal_ascii(data, off)
        }
    };

    let sample_name = read_str(layout.sample_off);
    let method_name = read_str(layout.method_off)
        .map(|s| file_name(&s).to_owned())
        .filter(|s| !s.is_empty());
    let instrument = read_str(layout.inst_off);
    let unit = read_str(layout.unit_off).unwrap_or_default();
    let signal = layout.sig_off.and_then(read_str);
    let wavelength_nm = signal.as_deref().and_then(extract_signal_wavelength);
    let label = signal.unwrap_or_else(|| "Signal".to_owned());
    let points = build_ch_points(&raw_values, start_ms, end_ms, scale);

    Some(ParsedCh {
        sample_name,
        method_name,
        instrument,
        wavelength_nm,
        trace: Trace {
            label,
            x_unit: "min".to_owned(),
            y_unit: if unit.is_empty() {
                "intensity".to_owned()
            } else {
                unit
            },
            points,
        },
    })
}

fn build_ch_points(raw: &[f64], start_ms: f64, end_ms: f64, scale: f64) -> Vec<(f64, f64)> {
    let n = raw.len();
    let start_min = start_ms / 60_000.0;
    let end_min = end_ms / 60_000.0;
    let has_span = n > 1 && end_min > start_min;
    let step_min = if has_span {
        (end_min - start_min) / ((n - 1) as f64)
    } else {
        1.0 / 60.0
    };
    let origin = if has_span { start_min } else { 0.0 };

    raw.iter()
        .enumerate()
        .map(|(i, &y)| (origin + (i as f64) * step_min, y * scale))
        .collect()
}

fn decode_ch_delta(buf: &[u8]) -> Vec<f64> {
    let mut out = Vec::new();
    let mut accum = 0i64;
    let mut pos = 0usize;
    while pos + 2 <= buf.len() {
        if buf[pos] != 0x10 {
            break;
        }
        let seg_count = buf[pos + 1] as usize;
        pos += 2;
        for _ in 0..seg_count {
            let Some(delta) = read_i16_be(buf, pos) else {
                return out;
            };
            pos += 2;
            if delta == i16::MIN {
                let Some(abs_val) = read_i32_be(buf, pos) else {
                    return out;
                };
                pos += 4;
                accum = i64::from(abs_val);
            } else {
                accum += i64::from(delta);
            }
            out.push(accum as f64);
        }
    }
    out
}

fn decode_ch_double_delta(buf: &[u8]) -> Vec<f64> {
    let mut out = Vec::new();
    let mut acc0 = 0i64;
    let mut acc1 = 0i64;
    let mut pos = 0usize;
    while let Some(d2) = read_i16_be(buf, pos) {
        pos += 2;
        if d2 == i16::MAX {
            let Some(hi) = read_i16_be(buf, pos) else {
                break;
            };
            let Some(lo) = read_i32_be(buf, pos + 2) else {
                break;
            };
            pos += 6;
            acc0 = (i64::from(hi) << 32) | i64::from(lo as u32);
            acc1 = 0;
        } else {
            acc1 += i64::from(d2);
            acc0 += acc1;
        }
        out.push(acc0 as f64);
    }
    out
}

fn parse_leading_pos_f64(s: &str) -> Option<f64> {
    let end = s
        .find(|c: char| !c.is_ascii_digit() && c != '.')
        .unwrap_or(s.len());
    s[..end].parse::<f64>().ok().filter(|&v| v > 0.0)
}

fn extract_signal_wavelength(desc: &str) -> Option<f64> {
    if let Some((_, after)) = desc.split_once("Sig=")
        && let Some(wl) = parse_leading_pos_f64(after)
    {
        return Some(wl);
    }
    let lower = desc.to_ascii_lowercase();
    let idx = lower.find("wavelength")?;
    let rest = desc[idx + "wavelength".len()..]
        .trim_start()
        .strip_prefix('=')?
        .trim_start();
    parse_leading_pos_f64(rest)
}

fn read_pascal_ascii(data: &[u8], offset: usize) -> Option<String> {
    let len = *data.get(offset)? as usize;
    let raw = data.get(offset + 1..offset + 1 + len)?;
    let text = String::from_utf8_lossy(raw).trim().to_owned();
    (!text.is_empty()).then_some(text)
}

fn read_pascal_utf16le(data: &[u8], offset: usize) -> Option<String> {
    let chars = *data.get(offset)? as usize;
    let byte_len = chars.checked_mul(2)?;
    let raw = data.get(offset + 1..offset + 1 + byte_len)?;
    let text = decode_utf16le(raw).trim().to_owned();
    (!text.is_empty()).then_some(text)
}
