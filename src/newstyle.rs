//! Pure-Rust reader for Newstyle (Qingbohua) NU3000 prep-HPLC XML datasets.
//!
//! Each NU3000 run stores paired XML files under `Data/`:
//! - `<uuid>.Data` and `<uuid>.FilterData`: whitespace-separated intensity samples (`<Data unit="mAu" count="...">`)
//! - `<uuid>.icm`: integration peak table and `samplingRate` in samples per millisecond (typically `0.01` = 10 Hz)
//! - `<uuid>.sys`: connected pump, UV-Vis detector, and fraction collector device metadata.

use std::collections::BTreeMap;

use crate::binary::{decode_text_auto, file_name, round_digits};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat};
use crate::xml::{elements, parse_xml, text_content};

const DEFAULT_SAMPLING_RATE_PER_MS: f64 = 0.01;

struct ParsedIcm {
    sampling_rate_per_ms: f64,
    peaks: Vec<Peak>,
}

struct ParsedDataFile {
    y_unit: String,
    values: Vec<f64>,
}

/// Parses a single Newstyle `.Data` or `.FilterData` XML file into an [`InstrumentReading`].
pub fn read_data(data: &[u8]) -> Result<InstrumentReading> {
    let parsed = parse_data_xml(data)?;
    let trace = build_trace(
        "UV",
        &parsed.y_unit,
        &parsed.values,
        DEFAULT_SAMPLING_RATE_PER_MS,
    );
    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::NewstyleNu3000,
        traces: vec![trace],
        ..InstrumentReading::default()
    })
}

/// Parses a single Newstyle `.icm` peak-table XML file into an [`InstrumentReading`].
pub fn read_icm(data: &[u8]) -> Result<InstrumentReading> {
    let parsed = parse_icm_xml(data)?;
    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::NewstyleNu3000,
        peaks: parsed.peaks,
        ..InstrumentReading::default()
    })
}

/// Parses a multi-file Newstyle NU3000 dataset bundle (`.Data`, `.FilterData`, `.icm`, `.sys`).
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let mut raw_data: BTreeMap<String, &[u8]> = BTreeMap::new();
    let mut filter_data: BTreeMap<String, &[u8]> = BTreeMap::new();
    let mut icm_files: BTreeMap<String, &[u8]> = BTreeMap::new();
    let mut sys_files: Vec<&[u8]> = Vec::new();

    for &(path, bytes) in files {
        let base = file_name(path);
        let lower = base.to_ascii_lowercase();
        if let Some(stem) = strip_suffix_ci(base, &lower, ".data") {
            raw_data.insert(stem.to_owned(), bytes);
        } else if let Some(stem) = strip_suffix_ci(base, &lower, ".filterdata") {
            filter_data.insert(stem.to_owned(), bytes);
        } else if let Some(stem) = strip_suffix_ci(base, &lower, ".icm") {
            icm_files.insert(stem.to_owned(), bytes);
        } else if lower.ends_with(".sys") {
            sys_files.push(bytes);
        }
    }

    let instrument = sys_files
        .iter()
        .find_map(|bytes| parse_sys_instrument(bytes));

    let mut stems: Vec<String> = raw_data
        .keys()
        .chain(filter_data.keys())
        .chain(icm_files.keys())
        .cloned()
        .collect();
    stems.sort();
    stems.dedup();

    let mut traces = Vec::new();
    let mut peaks = Vec::new();

    for (idx, stem) in stems.iter().enumerate() {
        let icm = icm_files
            .get(stem)
            .map(|bytes| parse_icm_xml(bytes))
            .transpose()?;
        let sampling_rate = icm
            .as_ref()
            .map_or(DEFAULT_SAMPLING_RATE_PER_MS, |p| p.sampling_rate_per_ms);

        if let Some(bytes) = raw_data.get(stem).or_else(|| filter_data.get(stem)) {
            let parsed = parse_data_xml(bytes)?;
            let label = format!("UV {}", idx + 1);
            traces.push(build_trace(
                &label,
                &parsed.y_unit,
                &parsed.values,
                sampling_rate,
            ));
        }

        if let Some(parsed_icm) = icm
            && peaks.is_empty()
            && !parsed_icm.peaks.is_empty()
        {
            peaks = parsed_icm.peaks;
        }
    }

    if traces.is_empty() && peaks.is_empty() {
        return Err(InstrumentError::InvalidNewstyle(
            "bundle contains no readable .Data, .FilterData, or .icm files".to_owned(),
        ));
    }

    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::NewstyleNu3000,
        instrument,
        method: InstrumentMethod::default(),
        peaks,
        traces,
        ..InstrumentReading::default()
    })
}

fn strip_suffix_ci<'a>(base: &'a str, lower: &str, suffix: &str) -> Option<&'a str> {
    lower
        .ends_with(suffix)
        .then(|| &base[..base.len() - suffix.len()])
}

fn parse_data_xml(data: &[u8]) -> Result<ParsedDataFile> {
    let text = decode_text_auto(data);
    let doc = parse_xml(&text).map_err(|err| InstrumentError::Xml(err.to_string()))?;
    let root = doc.root_element();
    let tag = root.tag_name().name();
    if tag != "Data" && tag != "FilterData" {
        return Err(InstrumentError::InvalidNewstyle(format!(
            "unexpected root tag <{tag}>, expected <Data> or <FilterData>"
        )));
    }
    let y_unit = root
        .attribute("unit")
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or("mAU")
        .to_owned();
    let body = text_content(root);
    let values: Vec<f64> = body
        .split_ascii_whitespace()
        .filter_map(|tok| tok.parse::<f64>().ok())
        .filter(|v| v.is_finite())
        .collect();
    if values.is_empty() {
        return Err(InstrumentError::InvalidNewstyle(
            "no numeric points in <Data> element".to_owned(),
        ));
    }
    Ok(ParsedDataFile { y_unit, values })
}

fn parse_icm_xml(data: &[u8]) -> Result<ParsedIcm> {
    let text = decode_text_auto(data);
    let doc = parse_xml(&text).map_err(|err| InstrumentError::Xml(err.to_string()))?;
    let chrom = elements(doc.root(), "Chrom").next().ok_or_else(|| {
        InstrumentError::InvalidNewstyle("missing <Chrom> element in .icm file".to_owned())
    })?;

    let sampling_rate_per_ms = chrom
        .attribute("samplingRate")
        .and_then(|s| s.trim().parse::<f64>().ok())
        .filter(|&v| v > 0.0 && v.is_finite())
        .unwrap_or(DEFAULT_SAMPLING_RATE_PER_MS);

    let mut peaks = Vec::new();
    for node in elements(chrom, "Peak") {
        let rt_min = elements(node, "Feature")
            .next()
            .and_then(|f| f.attribute("SummitIndex"))
            .and_then(|s| s.trim().parse::<f64>().ok())
            .map(|idx| round_digits(idx / (sampling_rate_per_ms * 60_000.0), 4));
        let area = parse_opt_f64(node.attribute("AREA"));
        let area_pct = parse_opt_f64(node.attribute("AREA_PERCENT"));
        let height = parse_opt_f64(node.attribute("HEIGHT"));
        let compound = node
            .attribute("name")
            .map(str::trim)
            .filter(|s| !s.is_empty() && *s != "N.A.")
            .map(ToOwned::to_owned);
        peaks.push(Peak {
            rt_min,
            area,
            area_pct,
            height,
            compound,
        });
    }

    Ok(ParsedIcm {
        sampling_rate_per_ms,
        peaks,
    })
}

fn parse_sys_instrument(data: &[u8]) -> Option<String> {
    let text = decode_text_auto(data);
    let doc = parse_xml(&text).ok()?;
    let mut vendor: Option<String> = None;
    let mut devices = Vec::new();
    for dev in elements(doc.root(), "Device") {
        let cat = dev.attribute("Category").unwrap_or("");
        if matches!(cat, "Pump" | "Detector" | "Collector") {
            if vendor.is_none() {
                vendor = dev
                    .attribute("Manufactory")
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(ToOwned::to_owned);
            }
            if let Some(name) = dev
                .attribute("Name")
                .map(str::trim)
                .filter(|s| !s.is_empty())
                && !devices.iter().any(|d: &String| d == name)
            {
                devices.push(name.to_owned());
            }
        }
    }
    match (vendor, devices.is_empty()) {
        (Some(v), false) => Some(format!("{v} ({})", devices.join(", "))),
        (Some(v), true) => Some(v),
        (None, false) => Some(devices.join(", ")),
        (None, true) => None,
    }
}

fn parse_opt_f64(attr: Option<&str>) -> Option<f64> {
    attr.and_then(|s| s.trim().parse::<f64>().ok())
        .filter(|v| v.is_finite())
}

fn build_trace(label: &str, y_unit: &str, values: &[f64], sampling_rate_per_ms: f64) -> Trace {
    let step_min = 1.0 / (sampling_rate_per_ms * 60_000.0);
    let points = values
        .iter()
        .enumerate()
        .map(|(i, &y)| ((i as f64) * step_min, y))
        .collect();
    Trace {
        label: label.to_owned(),
        x_unit: "min".to_owned(),
        y_unit: y_unit.to_owned(),
        points,
    }
}
