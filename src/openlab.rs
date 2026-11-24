//! Agilent OpenLab CDS reader: `.acaml` XML metadata and zipped `.dx`/`.rx` container inspection.

use std::collections::{BTreeMap, BTreeSet};

use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use serde::{Deserialize, Serialize};

use crate::binary::{
    append_xml_ref, decode_text_auto, file_name, round_digits, xml_attr, xml_local_name,
};
use crate::chemstation::parse_chemstation_ch;
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat};
use crate::zip_archive::{ZipArchive, ZipEntryInfo};

const NO_TIME_PREFIX: &str = "1601-01-01";
const INJECTION_ACAML_PATH: &str = "Base/InjectionACAML";

/// Metadata for one injection row inside an OpenLab `.acaml` result set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OpenLabInjection {
    /// Injection identifier (`InjectionId`).
    pub id: Option<String>,
    /// Data file stem without `.dx` (`RawDataFileName`).
    pub data_file: String,
    /// Sample name (`SampleName`).
    pub sample_name: String,
    /// Sample role (`SampleType`, such as `"Sample"` or `"Blank"`).
    pub sample_type: String,
    /// Autosampler vial number (`VialNumber`).
    pub vial: String,
    /// Acquisition instrument (`InstrumentName`).
    pub instrument: Option<String>,
    /// UTC acquisition timestamp (`InjectionAcqDateTime`), excluding sentinel `1601-01-01`.
    pub acquired_utc: Option<String>,
    /// Acquisition method name (`AcqMethodName`).
    pub method: String,
    /// Processing method name (`DaMethodName`).
    pub processing_method: Option<String>,
    /// Completion status (`Completed`, `Aborted`, or `None`).
    pub status: Option<String>,
    /// Relative `.rx` result path (`ExternalResultPath`).
    pub result_file: Option<String>,
    /// Injection volume in microliters, when recorded.
    pub injection_volume_ul: Option<f64>,
}

/// Parsed OpenLab `.acaml` result set metadata and operator logins.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AcamlResultSet {
    /// Injection metadata rows in document order.
    pub injections: Vec<OpenLabInjection>,
    /// Distinct operator/system user names referenced in `<Username>` or `<CreatedByUser>`.
    pub logins: BTreeSet<String>,
}

/// Container contents and integrated peak table extracted from an OpenLab `.rx` archive.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RxInspection {
    /// Members inside the `.rx` ZIP archive.
    pub entries: Vec<ZipEntryInfo>,
    /// Integrated peaks parsed from `Base/InjectionACAML`, sorted by retention time.
    pub peaks: Vec<Peak>,
}

/// Container contents, injection metadata, and detector signals extracted from an OpenLab `.dx` archive.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DxInspection {
    /// Members inside the `.dx` ZIP archive.
    pub entries: Vec<ZipEntryInfo>,
    /// Sample name from `injection.acmd` or the ChemStation `.CH` header.
    pub sample_name: Option<String>,
    /// Acquisition method file name from `injection.acmd` or the `.CH` header.
    pub method_name: Option<String>,
    /// Instrument or workstation label from the `.CH` header, when available.
    pub instrument: Option<String>,
    /// Optical detection wavelength in nm parsed from a DAD/VWD signal description, when present.
    pub wavelength_nm: Option<f64>,
    /// Decoded detector chromatogram traces from `.CH` members.
    pub traces: Vec<Trace>,
}

#[derive(Default)]
struct MeasDataRow {
    status: Option<String>,
    result_file: Option<String>,
}

/// Parses an Agilent OpenLab `.acaml` XML document, matching `sequency_lab_data.instruments.openlab.result_set`.
pub fn parse_acaml(data: &[u8]) -> Result<AcamlResultSet> {
    let text = decode_text_auto(data);
    let mut reader = Reader::from_str(&text);
    reader.config_mut().trim_text(false);

    let mut saw_acaml_root = false;
    let mut stack: Vec<String> = Vec::new();
    let mut current_meas_id: Option<String> = None;
    let mut measured: BTreeMap<String, MeasDataRow> = BTreeMap::new();
    let mut raw_rows: Vec<OpenLabInjection> = Vec::new();
    let mut logins: BTreeSet<String> = BTreeSet::new();
    let mut text_buf = String::new();

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => {
                let name = xml_local_name(e.local_name());
                if name == "ACAML" {
                    saw_acaml_root = true;
                } else if name == "MeasData" {
                    current_meas_id = xml_attr(&e, "id");
                    if let Some(id) = current_meas_id.as_ref() {
                        measured.entry(id.clone()).or_default();
                    }
                } else if name == "InjectionMetaData" {
                    raw_rows.push(parse_injection_meta_attrs(&e));
                }
                text_buf.clear();
                stack.push(name);
            }
            Ok(Event::Empty(e)) => {
                let name = xml_local_name(e.local_name());
                if name == "ACAML" {
                    saw_acaml_root = true;
                } else if name == "InjectionMetaData" {
                    raw_rows.push(parse_injection_meta_attrs(&e));
                }
            }
            Ok(Event::Text(t)) => {
                text_buf.push_str(&t);
            }
            Ok(Event::GeneralRef(r)) => append_xml_ref(&mut text_buf, &r),
            Ok(Event::End(e)) => {
                let name = xml_local_name(e.local_name());
                let trimmed = text_buf.trim();
                if matches!(name.as_str(), "Username" | "CreatedByUser") {
                    logins.insert(trimmed.to_owned());
                }
                if let Some(id) = current_meas_id.as_ref()
                    && let Some(row) = measured.get_mut(id)
                {
                    if name == "Status" && parent_is_acq_status(&stack) && !trimmed.is_empty() {
                        row.status = Some(trimmed.to_owned());
                    } else if name == "ExternalResultPath" && !trimmed.is_empty() {
                        row.result_file = Some(trimmed.to_owned());
                    }
                }
                if name == "MeasData" {
                    current_meas_id = None;
                }
                text_buf.clear();
                stack.pop();
            }
            Ok(Event::Eof) => break,
            Ok(_) => {}
            Err(err) => return Err(InstrumentError::Xml(err.to_string())),
        }
    }

    if !saw_acaml_root {
        return Err(InstrumentError::InvalidOpenLab(
            "XML document has no ACAML root element".to_owned(),
        ));
    }

    for inj in &mut raw_rows {
        if let Some(id) = inj.id.as_ref()
            && let Some(m) = measured.get(id)
        {
            inj.status.clone_from(&m.status);
            inj.result_file.clone_from(&m.result_file);
        }
    }

    Ok(AcamlResultSet {
        injections: raw_rows,
        logins,
    })
}

/// Inspects an OpenLab `.rx` ZIP container and extracts peaks from `Base/InjectionACAML`.
pub fn inspect_rx(data: &[u8]) -> Result<RxInspection> {
    let archive = ZipArchive::parse(data)?;
    let entries = archive.entry_infos();
    let peaks = match archive.read_by_name(INJECTION_ACAML_PATH)? {
        Some(xml_bytes) => parse_injection_acaml_peaks(&xml_bytes)?,
        None => Vec::new(),
    };
    Ok(RxInspection { entries, peaks })
}

/// Parses peaks directly from an OpenLab `.rx` archive, matching `sequency_lab_data.instruments.openlab.peaks`.
pub fn parse_rx_peaks(data: &[u8]) -> Result<Vec<Peak>> {
    let archive = ZipArchive::parse(data)?;
    let xml_bytes = archive.read_by_name(INJECTION_ACAML_PATH)?.ok_or_else(|| {
        InstrumentError::InvalidOpenLab("missing Base/InjectionACAML in .rx archive".to_owned())
    })?;
    parse_injection_acaml_peaks(&xml_bytes)
}

/// Inspects an OpenLab `.dx` ZIP container, reading `injection.acmd` and `.CH` detector signals.
pub fn inspect_dx(data: &[u8]) -> Result<DxInspection> {
    let archive = ZipArchive::parse(data)?;
    let entries = archive.entry_infos();

    let (mut sample_name, mut method_name) =
        match archive.read_by_name_ignore_ascii_case("injection.acmd")? {
            Some(acmd) => parse_acmd_metadata(&acmd),
            None => (None, None),
        };

    let mut ch_names: Vec<String> = entries
        .iter()
        .filter(|e| has_ext_ignore_case(&e.name, ".ch"))
        .map(|e| e.name.clone())
        .collect();
    ch_names.sort();

    let mut traces = Vec::new();
    let mut instrument = None;
    let mut wavelength_nm = None;

    for member in ch_names {
        if let Some(ch_bytes) = archive.read_by_name(&member)?
            && let Some(parsed) = parse_chemstation_ch(&ch_bytes)
        {
            if sample_name.is_none() {
                sample_name = parsed.sample_name;
            }
            if method_name.is_none() {
                method_name = parsed.method_name;
            }
            if instrument.is_none() {
                instrument = parsed.instrument;
            }
            if wavelength_nm.is_none() {
                wavelength_nm = parsed.wavelength_nm;
            }
            traces.push(parsed.trace);
        }
    }

    Ok(DxInspection {
        entries,
        sample_name,
        method_name,
        instrument,
        wavelength_nm,
        traces,
    })
}

/// Converts an OpenLab `.acaml` file into a unified [`InstrumentReading`].
pub fn read_acaml(data: &[u8]) -> Result<InstrumentReading> {
    let rs = parse_acaml(data)?;
    Ok(reading_from_acaml(&rs))
}

fn reading_from_acaml(rs: &AcamlResultSet) -> InstrumentReading {
    let primary = choose_primary_injection(&rs.injections);
    let technique = primary
        .and_then(|i| infer_openlab_technique(i.instrument.as_deref(), Some(&i.method)))
        .unwrap_or(Technique::Hplc);

    InstrumentReading {
        technique,
        vendor_format: VendorFormat::OpenLabAcaml,
        acquired_at: primary.and_then(|i| i.acquired_utc.clone()),
        sample_name: primary
            .map(|i| i.sample_name.as_str())
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned),
        instrument: primary.and_then(|i| i.instrument.clone()),
        method: InstrumentMethod {
            name: primary
                .map(|i| i.method.as_str())
                .filter(|s| !s.is_empty())
                .map(ToOwned::to_owned),
            injection_volume_ul: primary.and_then(|i| i.injection_volume_ul),
            ..InstrumentMethod::default()
        },
        peaks: Vec::new(),
        traces: Vec::new(),
    }
}

/// Converts an OpenLab `.rx` archive into a unified [`InstrumentReading`].
pub fn read_rx(data: &[u8]) -> Result<InstrumentReading> {
    let inspection = inspect_rx(data)?;
    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::OpenLabRx,
        peaks: inspection.peaks,
        ..InstrumentReading::default()
    })
}

/// Converts an OpenLab `.dx` archive into a unified [`InstrumentReading`].
pub fn read_dx(data: &[u8]) -> Result<InstrumentReading> {
    let dx = inspect_dx(data)?;
    let technique = infer_openlab_technique(dx.instrument.as_deref(), dx.method_name.as_deref())
        .unwrap_or(Technique::Hplc);

    Ok(InstrumentReading {
        technique,
        vendor_format: VendorFormat::OpenLabDx,
        acquired_at: None,
        sample_name: dx.sample_name,
        instrument: dx.instrument,
        method: InstrumentMethod {
            name: dx.method_name,
            wavelength_nm: dx.wavelength_nm,
            ..InstrumentMethod::default()
        },
        peaks: Vec::new(),
        traces: dx.traces,
    })
}

/// Combines an OpenLab result-set bundle (`.acaml` plus optional `.rx` and `.dx`) into one [`InstrumentReading`].
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let acaml_bytes = files
        .iter()
        .find(|(name, _)| has_ext_ignore_case(name, ".acaml"))
        .map(|(_, b)| *b);

    let parsed_acaml = match acaml_bytes {
        Some(bytes) => Some(parse_acaml(bytes)?),
        None => None,
    };

    let mut reading = match parsed_acaml.as_ref() {
        Some(rs) => reading_from_acaml(rs),
        None => InstrumentReading {
            technique: Technique::Hplc,
            vendor_format: VendorFormat::OpenLabDx,
            ..InstrumentReading::default()
        },
    };

    let target_rx = parsed_acaml.as_ref().and_then(|rs| {
        choose_primary_injection(&rs.injections).and_then(|i| i.result_file.clone())
    });

    let rx_file = files.iter().find(|(path, _)| {
        if !has_ext_ignore_case(path, ".rx") {
            return false;
        }
        target_rx
            .as_deref()
            .is_none_or(|want| file_name(path).eq_ignore_ascii_case(file_name(want)))
    });
    if let Some((_, rx_bytes)) = rx_file {
        let rx = inspect_rx(rx_bytes)?;
        reading.peaks = rx.peaks;
    }

    let dx_candidates: Vec<&(&str, &[u8])> = files
        .iter()
        .filter(|(path, _)| {
            let base = file_name(path);
            has_ext_ignore_case(base, ".dx") && !starts_with_ignore_case(base, "snapshot-")
        })
        .collect();
    if dx_candidates.len() == 1 {
        let dx = inspect_dx(dx_candidates[0].1)?;
        if reading.sample_name.is_none() {
            reading.sample_name = dx.sample_name;
        }
        if reading.instrument.is_none() {
            reading.instrument = dx.instrument;
        }
        if reading.method.name.is_none() {
            reading.method.name = dx.method_name;
        }
        if reading.method.wavelength_nm.is_none() {
            reading.method.wavelength_nm = dx.wavelength_nm;
        }
        reading.traces = dx.traces;
    }

    Ok(reading)
}

fn choose_primary_injection(injections: &[OpenLabInjection]) -> Option<&OpenLabInjection> {
    injections
        .iter()
        .find(|i| i.sample_type.eq_ignore_ascii_case("Sample"))
        .or_else(|| injections.first())
}

fn infer_openlab_technique(instrument: Option<&str>, method: Option<&str>) -> Option<Technique> {
    for text in [instrument, method].into_iter().flatten() {
        let upper = text.to_ascii_uppercase();
        if upper.contains("LCMS")
            || upper.contains("GCMS")
            || upper.contains("MSD")
            || upper.contains("QTOF")
        {
            return Some(Technique::LcMs);
        }
        if upper.contains("HPLC")
            || upper.contains("LC")
            || upper.contains("DAD")
            || upper.contains("VWD")
            || upper.contains("GC")
        {
            return Some(Technique::Hplc);
        }
    }
    None
}

fn parse_injection_meta_attrs(e: &BytesStart<'_>) -> OpenLabInjection {
    let raw_data_file = xml_attr(e, "RawDataFileName").unwrap_or_default();
    let data_file = strip_suffix_ignore_case(&raw_data_file, ".dx").to_owned();
    let sample_name = xml_attr(e, "SampleName")
        .map(|s| s.trim().to_owned())
        .unwrap_or_default();
    let sample_type = xml_attr(e, "SampleType").unwrap_or_default();
    let vial = xml_attr(e, "VialNumber").unwrap_or_default();
    let instrument = xml_attr(e, "InstrumentName").filter(|s| !s.is_empty());
    let acquired_utc = xml_attr(e, "InjectionAcqDateTime")
        .filter(|at| !at.is_empty() && !at.starts_with(NO_TIME_PREFIX));
    let method = xml_attr(e, "AcqMethodName").unwrap_or_default();
    let processing_method = xml_attr(e, "DaMethodName").filter(|s| !s.is_empty());
    let injection_volume_ul = xml_attr(e, "InjectionVolume")
        .or_else(|| xml_attr(e, "InjVolume"))
        .and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v >= 0.0);

    OpenLabInjection {
        id: xml_attr(e, "InjectionId"),
        data_file,
        sample_name,
        sample_type,
        vial,
        instrument,
        acquired_utc,
        method,
        processing_method,
        status: None,
        result_file: None,
        injection_volume_ul,
    }
}

fn parse_injection_acaml_peaks(data: &[u8]) -> Result<Vec<Peak>> {
    let text = decode_text_auto(data);
    let mut reader = Reader::from_str(&text);
    reader.config_mut().trim_text(false);

    let mut in_compound = false;
    let mut compound_name = String::new();
    let mut compound_peak_ids: Vec<String> = Vec::new();
    let mut named_peaks: BTreeMap<String, String> = BTreeMap::new();

    let mut current_peak_id: Option<String> = None;
    let mut current_peak: Option<Peak> = None;
    let mut peaks: Vec<(Option<String>, Peak)> = Vec::new();
    let mut text_buf = String::new();

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => {
                let name = xml_local_name(e.local_name());
                match name.as_str() {
                    "InjectionCompound" => {
                        in_compound = true;
                        compound_name.clear();
                        compound_peak_ids.clear();
                    }
                    "Peak_ID" if in_compound => {
                        if let Some(id) = xml_attr(&e, "id") {
                            compound_peak_ids.push(id);
                        }
                    }
                    "Peak" => {
                        current_peak_id = xml_attr(&e, "id");
                        current_peak = Some(Peak {
                            rt_min: None,
                            area: None,
                            area_pct: None,
                            height: None,
                            compound: None,
                        });
                    }
                    _ => {
                        if let Some(peak) = current_peak.as_mut() {
                            record_peak_val_attr(&name, &e, peak);
                        }
                    }
                }
                text_buf.clear();
            }
            Ok(Event::Empty(e)) => {
                let name = xml_local_name(e.local_name());
                if name == "Peak_ID" && in_compound {
                    if let Some(id) = xml_attr(&e, "id") {
                        compound_peak_ids.push(id);
                    }
                } else if let Some(peak) = current_peak.as_mut() {
                    record_peak_val_attr(&name, &e, peak);
                }
            }
            Ok(Event::Text(t)) => {
                text_buf.push_str(&t);
            }
            Ok(Event::GeneralRef(r)) => append_xml_ref(&mut text_buf, &r),
            Ok(Event::End(e)) => {
                let name = xml_local_name(e.local_name());
                match name.as_str() {
                    "CompoundName" if in_compound => {
                        compound_name = text_buf.trim().to_owned();
                    }
                    "InjectionCompound" => {
                        in_compound = false;
                        if !compound_name.is_empty() {
                            for pid in compound_peak_ids.drain(..) {
                                named_peaks.insert(pid, compound_name.clone());
                            }
                        }
                    }
                    "Peak" => {
                        if let Some(peak) = current_peak.take() {
                            peaks.push((current_peak_id.take(), peak));
                        }
                    }
                    _ => {}
                }
                text_buf.clear();
            }
            Ok(Event::Eof) => break,
            Ok(_) => {}
            Err(err) => return Err(InstrumentError::Xml(err.to_string())),
        }
    }

    let mut resolved: Vec<Peak> = peaks
        .into_iter()
        .map(|(pid, mut peak)| {
            if let Some(id) = pid {
                peak.compound = named_peaks.get(&id).cloned();
            }
            peak
        })
        .collect();

    resolved.sort_by(|a, b| match (a.rt_min, b.rt_min) {
        (Some(ra), Some(rb)) => ra.total_cmp(&rb),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    Ok(resolved)
}

fn record_peak_val_attr(tag: &str, e: &BytesStart<'_>, peak: &mut Peak) {
    let Some(raw) = xml_attr(e, "val") else {
        return;
    };
    let Ok(val) = raw.parse::<f64>() else {
        return;
    };
    if !val.is_finite() {
        return;
    }
    match tag {
        "RetentionTime" => peak.rt_min = Some(round_digits(val, 3)),
        "Area" => peak.area = Some(round_digits(val, 3)),
        "AreaPercent" => peak.area_pct = Some(round_digits(val, 2)),
        "Height" => peak.height = Some(round_digits(val, 3)),
        _ => {}
    }
}

fn parse_acmd_metadata(data: &[u8]) -> (Option<String>, Option<String>) {
    let text = decode_text_auto(data);
    let mut reader = Reader::from_str(&text);
    let mut sample = None;
    let mut method = None;
    let mut current: Option<String> = None;
    let mut buf = String::new();

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => {
                current = Some(xml_local_name(e.local_name()));
                buf.clear();
            }
            Ok(Event::Text(t)) => {
                buf.push_str(&t);
            }
            Ok(Event::GeneralRef(r)) => append_xml_ref(&mut buf, &r),
            Ok(Event::End(e)) => {
                let name = xml_local_name(e.local_name());
                let trimmed = buf.trim();
                if current.as_deref() == Some(name.as_str()) && !trimmed.is_empty() {
                    if name == "SampleName" && sample.is_none() {
                        sample = Some(trimmed.to_owned());
                    } else if name == "AcquisitionMethod" && method.is_none() {
                        let m = file_name(trimmed);
                        if !m.is_empty() {
                            method = Some(m.to_owned());
                        }
                    }
                }
                current = None;
                buf.clear();
            }
            Ok(Event::Eof) | Err(_) => break,
            _ => {}
        }
    }

    (sample, method)
}

fn parent_is_acq_status(stack: &[String]) -> bool {
    stack
        .len()
        .checked_sub(2)
        .and_then(|idx| stack.get(idx))
        .is_some_and(|p| p == "AcquitionStatus" || p == "AcquisitionStatus")
}

fn strip_suffix_ignore_case<'a>(s: &'a str, suffix: &str) -> &'a str {
    let start = s.len().saturating_sub(suffix.len());
    if s.len() >= suffix.len()
        && s.is_char_boundary(start)
        && s[start..].eq_ignore_ascii_case(suffix)
    {
        &s[..start]
    } else {
        s
    }
}

fn has_ext_ignore_case(path: &str, ext: &str) -> bool {
    let start = path.len().saturating_sub(ext.len());
    path.len() >= ext.len()
        && path.is_char_boundary(start)
        && path[start..].eq_ignore_ascii_case(ext)
}

fn starts_with_ignore_case(s: &str, prefix: &str) -> bool {
    s.len() >= prefix.len()
        && s.is_char_boundary(prefix.len())
        && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}
