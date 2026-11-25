//! Agilent OpenLab CDS reader: `.acaml` XML metadata and zipped `.dx`/`.rx` container inspection.

use std::collections::{BTreeMap, BTreeSet};

use roxmltree::Node;
use serde::{Deserialize, Serialize};

use crate::binary::{decode_text_auto, file_name, round_digits};
use crate::chemstation::parse_chemstation_ch;
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat};
use crate::xml::{elements, parse_xml, text_content, trimmed_text};
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

/// Parses an Agilent OpenLab `.acaml` XML document, matching `sequency_lab_data.instruments.openlab.result_set`.
pub fn parse_acaml(data: &[u8]) -> Result<AcamlResultSet> {
    let text = decode_text_auto(data);
    let doc = parse_xml(&text).map_err(|err| InstrumentError::Xml(err.to_string()))?;
    if doc.root_element().tag_name().name() != "ACAML" {
        return Err(InstrumentError::InvalidOpenLab(
            "XML document has no ACAML root element".to_owned(),
        ));
    }

    let mut measured: BTreeMap<&str, (Option<String>, Option<String>)> = BTreeMap::new();
    for meas in elements(doc.root(), "MeasData") {
        let Some(id) = meas.attribute("id") else {
            continue;
        };
        let entry = measured.entry(id).or_default();
        if let Some(status) = elements(meas, "Status")
            .filter(|n| {
                n.parent().is_some_and(|p| {
                    matches!(p.tag_name().name(), "AcquitionStatus" | "AcquisitionStatus")
                })
            })
            .filter_map(trimmed_text)
            .next_back()
        {
            entry.0 = Some(status);
        }
        if let Some(result_file) = elements(meas, "ExternalResultPath")
            .filter_map(trimmed_text)
            .next_back()
        {
            entry.1 = Some(result_file);
        }
    }

    let injections = elements(doc.root(), "InjectionMetaData")
        .map(|n| parse_injection_meta(n, &measured))
        .collect();

    let logins = doc
        .descendants()
        .filter(|n| n.is_element() && matches!(n.tag_name().name(), "Username" | "CreatedByUser"))
        .map(|n| text_content(n).trim().to_owned())
        .collect();

    Ok(AcamlResultSet { injections, logins })
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
            sample_name = sample_name.or(parsed.sample_name);
            method_name = method_name.or(parsed.method_name);
            instrument = instrument.or(parsed.instrument);
            wavelength_nm = wavelength_nm.or(parsed.wavelength_nm);
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
        reading.sample_name = reading.sample_name.or(dx.sample_name);
        reading.instrument = reading.instrument.or(dx.instrument);
        reading.method.name = reading.method.name.or(dx.method_name);
        reading.method.wavelength_nm = reading.method.wavelength_nm.or(dx.wavelength_nm);
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

fn parse_injection_meta(
    node: Node<'_, '_>,
    measured: &BTreeMap<&str, (Option<String>, Option<String>)>,
) -> OpenLabInjection {
    let id = node.attribute("InjectionId");
    let (status, result_file) = id
        .and_then(|id| measured.get(id))
        .cloned()
        .unwrap_or((None, None));
    let raw_data_file = node.attribute("RawDataFileName").unwrap_or("");
    let injection_volume_ul = node
        .attribute("InjectionVolume")
        .or_else(|| node.attribute("InjVolume"))
        .and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite() && *v >= 0.0);

    OpenLabInjection {
        id: id.map(ToOwned::to_owned),
        data_file: strip_suffix_ignore_case(raw_data_file, ".dx").to_owned(),
        sample_name: node
            .attribute("SampleName")
            .map(str::trim)
            .unwrap_or("")
            .to_owned(),
        sample_type: node.attribute("SampleType").unwrap_or("").to_owned(),
        vial: node.attribute("VialNumber").unwrap_or("").to_owned(),
        instrument: node
            .attribute("InstrumentName")
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned),
        acquired_utc: node
            .attribute("InjectionAcqDateTime")
            .filter(|at| !at.is_empty() && !at.starts_with(NO_TIME_PREFIX))
            .map(ToOwned::to_owned),
        method: node.attribute("AcqMethodName").unwrap_or("").to_owned(),
        processing_method: node
            .attribute("DaMethodName")
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned),
        status,
        result_file,
        injection_volume_ul,
    }
}

fn parse_injection_acaml_peaks(data: &[u8]) -> Result<Vec<Peak>> {
    let text = decode_text_auto(data);
    let doc = parse_xml(&text).map_err(|err| InstrumentError::Xml(err.to_string()))?;

    let mut named_peaks: BTreeMap<&str, String> = BTreeMap::new();
    for compound in elements(doc.root(), "InjectionCompound") {
        let Some(name) = elements(compound, "CompoundName")
            .filter_map(trimmed_text)
            .next_back()
        else {
            continue;
        };
        for pid in elements(compound, "Peak_ID") {
            if let Some(id) = pid.attribute("id") {
                named_peaks.insert(id, name.clone());
            }
        }
    }

    let mut peaks: Vec<Peak> = elements(doc.root(), "Peak")
        .filter_map(|node| {
            let rt_min = peak_val(node, "RetentionTime", 3);
            let area = peak_val(node, "Area", 3);
            let area_pct = peak_val(node, "AreaPercent", 2);
            let height = peak_val(node, "Height", 3);
            if rt_min.is_none() && area.is_none() && area_pct.is_none() && height.is_none() {
                return None;
            }
            Some(Peak {
                rt_min,
                area,
                area_pct,
                height,
                compound: node
                    .attribute("id")
                    .and_then(|id| named_peaks.get(id))
                    .cloned(),
            })
        })
        .collect();

    peaks.sort_by(|a, b| match (a.rt_min, b.rt_min) {
        (Some(ra), Some(rb)) => ra.total_cmp(&rb),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    Ok(peaks)
}

fn peak_val(peak: Node<'_, '_>, tag: &str, decimals: i32) -> Option<f64> {
    elements(peak, tag)
        .filter_map(|n| n.attribute("val")?.parse::<f64>().ok())
        .rfind(|v| v.is_finite())
        .map(|v| round_digits(v, decimals))
}

fn parse_acmd_metadata(data: &[u8]) -> (Option<String>, Option<String>) {
    let text = decode_text_auto(data);
    let Ok(doc) = parse_xml(&text) else {
        return (None, None);
    };
    let sample = elements(doc.root(), "SampleName").find_map(trimmed_text);
    let method = elements(doc.root(), "AcquisitionMethod")
        .find_map(trimmed_text)
        .map(|s| file_name(&s).to_owned())
        .filter(|m| !m.is_empty());
    (sample, method)
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
    strip_suffix_ignore_case(path, ext).len() != path.len()
}

fn starts_with_ignore_case(s: &str, prefix: &str) -> bool {
    s.len() >= prefix.len()
        && s.is_char_boundary(prefix.len())
        && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_acaml_merges_duplicate_meas_data_and_joins_split_text() {
        let xml = b"  \n<?xml version=\"1.0\" encoding=\"utf-8\"?>
<ACAML>
  <Doc>
    <Username>ab<!--comment-->cd</Username>
    <MeasData id=\"inj-1\">
      <AcquisitionStatus><Status>Completed</Status></AcquisitionStatus>
    </MeasData>
    <MeasData id=\"inj-1\">
      <ExternalResultPath>results/run1.rx</ExternalResultPath>
    </MeasData>
    <InjectionMetaData InjectionId=\"inj-1\" RawDataFileName=\"run1.dx\" SampleName=\"S1\" SampleType=\"Sample\" VialNumber=\"1\" AcqMethodName=\"M1\"/>
  </Doc>
</ACAML>\0\0";

        let rs = parse_acaml(xml).unwrap();
        assert!(rs.logins.contains("abcd"));
        assert_eq!(rs.injections.len(), 1);
        assert_eq!(rs.injections[0].status.as_deref(), Some("Completed"));
        assert_eq!(
            rs.injections[0].result_file.as_deref(),
            Some("results/run1.rx")
        );
    }

    #[test]
    fn parse_injection_acaml_peaks_skips_empty_peak() {
        let xml = b"<ACAML>
  <InjectionCompound>
    <CompoundName>Caffeine</CompoundName>
    <Peak_ID id=\"p1\"/>
    <Peak_ID id=\"p2\"/>
  </InjectionCompound>
  <Peak id=\"p1\"/>
  <Peak id=\"p2\">
    <RetentionTime val=\"2.3456\"/>
    <Area val=\"100.5\"/>
  </Peak>
</ACAML>";

        let peaks = parse_injection_acaml_peaks(xml).unwrap();
        assert_eq!(peaks.len(), 1);
        assert_eq!(peaks[0].rt_min, Some(2.346));
        assert_eq!(peaks[0].area, Some(100.5));
        assert_eq!(peaks[0].compound.as_deref(), Some("Caffeine"));
    }
}
