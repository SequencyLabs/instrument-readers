//! Parsing for Waters MassLynx text and binary metadata files (`_HEADER.TXT`, `_FUNCTNS.INF`,
//! `_extern.inf`).

use std::collections::BTreeMap;

use crate::binary::{file_name, read_f32_le};
use crate::error::{InstrumentError, Result};

use super::{
    FUNC_TYPE_DIODE_ARRAY, FUNC_TYPE_MASK, FUNCTNS_MZ_OFFSET, FUNCTNS_MZ_SLOTS, FUNCTNS_RECORD_LEN,
    IonPolarity, WatersExternFunction, WatersHeader,
};

pub(super) struct WatersBundle<'a> {
    pub(super) meta: WatersHeader,
    pub(super) extern_funcs: BTreeMap<u32, WatersExternFunction>,
    pub(super) functns_bytes: Option<&'a [u8]>,
    pub(super) func_dats: BTreeMap<u32, &'a [u8]>,
    pub(super) func_idxs: BTreeMap<u32, &'a [u8]>,
}

pub(super) struct WatersFunctionInfo {
    pub(super) func_type: u8,
    pub(super) labels: Vec<f64>,
}

pub(super) fn collect_bundle<'a>(files: &[(&str, &'a [u8])]) -> Result<WatersBundle<'a>> {
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
        .unwrap_or_default();
    let extern_funcs = extern_bytes.map(parse_extern_functions).unwrap_or_default();

    Ok(WatersBundle {
        meta,
        extern_funcs,
        functns_bytes,
        func_dats,
        func_idxs,
    })
}

pub(super) fn is_uv_function(
    has_functns: bool,
    func_info: Option<&WatersFunctionInfo>,
    polarity: Option<char>,
) -> bool {
    if has_functns {
        func_info.is_some_and(|info| info.func_type == FUNC_TYPE_DIODE_ARRAY)
    } else {
        polarity.is_none()
    }
}

pub(super) fn parse_func_number(upper_name: &str, suffix: &str) -> Option<u32> {
    let middle = upper_name.strip_prefix("_FUNC")?.strip_suffix(suffix)?;
    if middle.is_empty() || !middle.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    middle.parse::<u32>().ok()
}

pub(super) fn parse_header_txt(data: &[u8]) -> Result<WatersHeader> {
    let text = String::from_utf8_lossy(data);
    let mut sample_name = None;
    let mut acquired_date = None;
    let mut acquired_time = None;
    let mut instrument = None;
    let mut inlet_method = None;
    let mut ms_method = None;
    let mut calibrations = BTreeMap::new();
    let mut unnamed_cal_num: u32 = 0;
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
        if let Some(suffix) = key.strip_prefix("Cal Function") {
            let coeffs: Vec<f64> = value
                .split(',')
                .filter_map(|p| p.trim().parse::<f64>().ok())
                .collect();
            let func_num = match suffix
                .split_whitespace()
                .next()
                .and_then(|tok| tok.parse::<u32>().ok())
            {
                Some(n) => {
                    unnamed_cal_num = unnamed_cal_num.max(n);
                    n
                }
                None => {
                    unnamed_cal_num = unnamed_cal_num.saturating_add(1);
                    unnamed_cal_num
                }
            };
            calibrations.insert(func_num, coeffs);
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

    Ok(WatersHeader {
        sample_name,
        acquired_at,
        instrument,
        inlet_method,
        ms_method,
        calibrations,
    })
}

pub(super) fn parse_functns_inf(data: &[u8]) -> Vec<WatersFunctionInfo> {
    let num_funcs = data.len() / FUNCTNS_RECORD_LEN;
    let mut out = Vec::with_capacity(num_funcs);
    for f in 0..num_funcs {
        let record_start = f * FUNCTNS_RECORD_LEN;
        let func_type = data[record_start] & FUNC_TYPE_MASK;
        let base = record_start + FUNCTNS_MZ_OFFSET;
        let mut labels = Vec::new();
        for slot in 0..FUNCTNS_MZ_SLOTS {
            if let Some(val) = read_f32_le(data, base + slot * 4)
                && val != 0.0
                && val.is_finite()
            {
                labels.push(f64::from(val));
            }
        }
        out.push(WatersFunctionInfo { func_type, labels });
    }
    out
}

fn parse_section_function_number(heading: &str) -> Option<u32> {
    let idx = heading.find("Function ")?;
    heading[idx + "Function ".len()..]
        .split(|ch: char| !ch.is_ascii_digit())
        .next()
        .filter(|s| !s.is_empty())
        .and_then(|digits| digits.parse::<u32>().ok())
}

fn parse_extern_functions(data: &[u8]) -> BTreeMap<u32, WatersExternFunction> {
    let text = String::from_utf8_lossy(data);
    let mut out: BTreeMap<u32, WatersExternFunction> = BTreeMap::new();
    let mut current_func: Option<u32> = None;
    let mut unnamed_section: u32 = 0;

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("Instrument Parameters") || trimmed.starts_with("Function ") {
            if let Some(func_num) = parse_section_function_number(trimmed) {
                unnamed_section = unnamed_section.max(func_num);
                current_func = Some(func_num);
            } else if trimmed.starts_with("Instrument Parameters") {
                unnamed_section = unnamed_section.saturating_add(1);
                current_func = Some(unnamed_section);
            } else {
                current_func = None;
            }
            if let Some(func_num) = current_func {
                out.entry(func_num).or_default();
            }
            continue;
        }
        if let Some(func_num) = current_func {
            let entry = out.entry(func_num).or_default();
            if trimmed.starts_with("Polarity")
                && let Some(pol) = extract_polarity(trimmed)
            {
                entry.polarity.get_or_insert(pol);
            }
            if let Some((k, v)) = trimmed.split_once('\t') {
                let key = k.trim();
                let val = v.trim();
                if !key.is_empty() && !val.is_empty() {
                    entry.params.push((key.to_owned(), val.to_owned()));
                }
            }
        }
    }

    out
}

fn extract_polarity(line: &str) -> Option<IonPolarity> {
    match line.trim().chars().last()? {
        '+' => Some(IonPolarity::Positive),
        '-' => Some(IonPolarity::Negative),
        _ => None,
    }
}
