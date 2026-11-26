//! Pure-Rust Agilent ChemStation `.CH` chromatogram and `.MS` spectral file decoders.
//!
//! Used both for standalone ChemStation `.D` folders (such as Agilent 1200 LC-MS) and for
//! `.CH` channels embedded inside Agilent OpenLab `.dx` archives.

use crate::binary::{
    decode_text_auto, decode_utf16le, file_name, read_f32_be, read_f64_be, read_i16_be,
    read_i32_be, read_u16_be, read_u16_le, read_u32_be,
};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Technique, Trace, VendorFormat};
use crate::xml::{elements, parse_xml, trimmed_text};

const MS_MAGIC: [u8; 4] = [0x01, 0x32, 0x00, 0x00];
const MS_INT_POW8: [u32; 4] = [1, 8, 64, 512];

/// Returns `true` if `data` starts with a recognized ChemStation `.CH` magic header (`30`, `130`, `179`, `181`).
#[must_use]
pub fn is_chemstation_ch(data: &[u8]) -> bool {
    matches!(
        read_pascal_ascii(data, 0).as_deref(),
        Some("30" | "130" | "179" | "181")
    )
}

/// Returns `true` if `data` starts with the ChemStation `.MS` magic (`\x012\x00\x00`).
#[must_use]
pub fn is_chemstation_ms(data: &[u8]) -> bool {
    data.starts_with(&MS_MAGIC)
}

/// Parses a standalone ChemStation `.CH` file into an [`InstrumentReading`].
pub fn read_ch(data: &[u8]) -> Result<InstrumentReading> {
    let parsed = parse_chemstation_ch(data).ok_or_else(|| {
        InstrumentError::InvalidChemStation("unrecognized or truncated .CH chromatogram".to_owned())
    })?;
    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::ChemStation,
        sample_name: parsed.sample_name,
        instrument: parsed.instrument,
        method: InstrumentMethod {
            name: parsed.method_name,
            wavelength_nm: parsed.wavelength_nm,
            ..InstrumentMethod::default()
        },
        traces: vec![parsed.trace],
        ..InstrumentReading::default()
    })
}

/// Parses a ChemStation `.MS` (`MSD1.MS`) binary file into a Total Ion Current (TIC) [`InstrumentReading`].
pub fn read_ms(data: &[u8]) -> Result<InstrumentReading> {
    if !is_chemstation_ms(data) {
        return Err(InstrumentError::InvalidChemStation(
            "missing ChemStation .MS magic header".to_owned(),
        ));
    }

    let ms_type = read_pascal_ascii(data, 0x04).unwrap_or_default();
    let num_times: usize = if ms_type == "MSD Spectral File" {
        read_u32_be(data, 0x116).ok_or_else(|| {
            InstrumentError::InvalidChemStation("truncated LC-MS scan count".to_owned())
        })? as usize
    } else {
        usize::from(read_u16_le(data, 0x142).ok_or_else(|| {
            InstrumentError::InvalidChemStation("truncated GC-MS scan count".to_owned())
        })?)
    };

    let raw_start = usize::from(read_u16_be(data, 0x10A).ok_or_else(|| {
        InstrumentError::InvalidChemStation("truncated .MS data offset".to_owned())
    })?);
    let mut pos = raw_start.saturating_mul(2).saturating_sub(2);
    if pos >= data.len() {
        return Err(InstrumentError::InvalidChemStation(
            "invalid .MS data start offset".to_owned(),
        ));
    }

    let mut points = Vec::with_capacity(num_times.min(data.len().saturating_sub(pos) / 28));
    for _ in 0..num_times {
        if pos + 18 > data.len() {
            break;
        }
        let Some(time_ms) = read_u32_be(data, pos + 2) else {
            break;
        };
        let Some(pair_count) = read_u16_be(data, pos + 12) else {
            break;
        };
        let pairs_start = pos + 18;
        let pairs_bytes = usize::from(pair_count) * 4;
        if pairs_start + pairs_bytes > data.len() {
            break;
        }

        let mut tic = 0u64;
        let mut p = pairs_start;
        for _ in 0..pair_count {
            let int_enc = u16::from_be_bytes([data[p + 2], data[p + 3]]);
            let head = usize::from(int_enc >> 14);
            let tail = u32::from(int_enc & 0x3FFF);
            tic += u64::from(MS_INT_POW8[head]) * u64::from(tail);
            p += 4;
        }

        points.push((f64::from(time_ms) / 60_000.0, tic as f64));
        pos = pairs_start + pairs_bytes + 10;
    }

    if points.is_empty() {
        return Err(InstrumentError::InvalidChemStation(
            "no scans decoded from .MS file".to_owned(),
        ));
    }

    let acquired_at = read_pascal_ascii(data, 0xB2);
    let method_name = read_pascal_ascii(data, 0xE4)
        .map(|s| file_name(&s).to_owned())
        .filter(|s| !s.is_empty());

    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::ChemStation,
        acquired_at,
        method: InstrumentMethod {
            name: method_name,
            ..InstrumentMethod::default()
        },
        traces: vec![Trace {
            label: "MSD1 TIC".to_owned(),
            x_unit: "min".to_owned(),
            y_unit: "counts".to_owned(),
            points,
        }],
        ..InstrumentReading::default()
    })
}

/// Parses an Agilent ChemStation `.D` directory bundle (`.CH`, `.MS`, and optional `SAMPLE.XML`).
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let mut reading = InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::ChemStation,
        ..InstrumentReading::default()
    };

    for &(path, bytes) in files {
        let base_lower = file_name(path).to_ascii_lowercase();
        if base_lower == "sample.xml" {
            if let Some(sample_name) = parse_sample_xml(bytes) {
                reading.sample_name = Some(sample_name);
            }
        } else if base_lower.ends_with(".ch") && is_chemstation_ch(bytes) {
            if let Some(parsed) = parse_chemstation_ch(bytes) {
                if reading.sample_name.is_none() {
                    reading.sample_name = parsed.sample_name;
                }
                if reading.instrument.is_none() {
                    reading.instrument = parsed.instrument;
                }
                if reading.method.name.is_none() {
                    reading.method.name = parsed.method_name;
                }
                if reading.method.wavelength_nm.is_none() {
                    reading.method.wavelength_nm = parsed.wavelength_nm;
                }
                reading.traces.push(parsed.trace);
            }
        } else if base_lower.ends_with(".ms")
            && is_chemstation_ms(bytes)
            && let Ok(ms_reading) = read_ms(bytes)
        {
            reading.technique = Technique::LcMs;
            if reading.acquired_at.is_none() {
                reading.acquired_at = ms_reading.acquired_at;
            }
            if reading.method.name.is_none() {
                reading.method.name = ms_reading.method.name;
            }
            reading.traces.extend(ms_reading.traces);
        }
    }

    if reading.traces.is_empty() {
        return Err(InstrumentError::InvalidChemStation(
            "bundle contains no readable .CH or .MS files".to_owned(),
        ));
    }

    Ok(reading)
}

fn parse_sample_xml(bytes: &[u8]) -> Option<String> {
    let text = decode_text_auto(bytes);
    let doc = parse_xml(&text).ok()?;
    elements(doc.root(), "Name").find_map(trimmed_text)
}

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
                accum = accum.wrapping_add(i64::from(delta));
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
            acc1 = acc1.wrapping_add(i64::from(d2));
            acc0 = acc0.wrapping_add(acc1);
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
