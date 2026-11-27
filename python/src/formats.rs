//! Format-specific PyO3 wrappers for ANDI, Chromeleon, OMNIC, UV text, Thermo Finnigan,
//! Waters MassLynx, Agilent MassHunter, and MestReNova.

use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::PyBytes;
use sequency_instruments::{andi, chromeleon, masshunter, mnova, omnic, thermo, uvtext, waters};

use crate::{AcqusParams, Peak, Trace, py_result};

pub(crate) fn extract_file_bundle(files: &Bound<'_, PyAny>) -> PyResult<Vec<(String, Vec<u8>)>> {
    files.extract()
}

/// Parsed ANDI Chromatography netCDF-3 dataset.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct AndiDataset {
    pub sample_name: Option<String>,
    pub operator_name: Option<String>,
    pub source_file_reference: Option<String>,
    pub detection_method_name: Option<String>,
    pub injection_date_time_stamp: Option<String>,
    pub acquired_at: Option<String>,
    pub sampling_interval_s: Option<f64>,
    pub peaks: Vec<Peak>,
    pub traces: Vec<Trace>,
}

/// Parses an ANDI Chromatography netCDF-3 `.cdf` buffer.
#[pyfunction]
pub fn andi_parse_cdf(data: &[u8]) -> PyResult<AndiDataset> {
    let ds = py_result(andi::parse_cdf(data))?;
    Ok(AndiDataset {
        sample_name: ds.sample_name,
        operator_name: ds.operator_name,
        source_file_reference: ds.source_file_reference,
        detection_method_name: ds.detection_method_name,
        injection_date_time_stamp: ds.injection_date_time_stamp,
        acquired_at: ds.acquired_at,
        sampling_interval_s: ds.sampling_interval_s,
        peaks: ds.peaks.into_iter().map(|p| p.into_peak().into()).collect(),
        traces: ds.traces.into_iter().map(Into::into).collect(),
    })
}

/// Extracts the channel name from a Thermo Chromeleon `SignHdr` buffer.
#[pyfunction]
pub fn chromeleon_signal_name(data: &[u8]) -> Option<String> {
    chromeleon::signal_name(data)
}

/// Decodes a Thermo Chromeleon `SignHdr` `.raw` signal into `(minutes, values)`.
#[pyfunction]
pub fn chromeleon_decode_signal(data: &[u8]) -> PyResult<(Vec<f64>, Vec<f64>)> {
    py_result(chromeleon::decode_signal(data))
}

/// Lightweight header metadata from a Thermo Nicolet OMNIC `.spa` spectrum.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct OmnicHeader {
    pub title: String,
    pub y_unit: Option<String>,
    pub acquired_epoch_secs: u32,
}

/// Parsed Thermo Nicolet OMNIC `.spa` spectrum.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct OmnicSpectrum {
    pub wavenumbers: Vec<f64>,
    pub intensities: Vec<f64>,
    pub title: String,
    pub y_unit: Option<String>,
    pub acquired_epoch_secs: u32,
}

/// Parses only the header metadata of an OMNIC `.spa` spectrum buffer.
#[pyfunction]
pub fn omnic_parse_header(data: &[u8]) -> PyResult<OmnicHeader> {
    let h = py_result(omnic::parse_header(data))?;
    Ok(OmnicHeader {
        title: h.title,
        y_unit: h.y_unit,
        acquired_epoch_secs: h.acquired_epoch_secs,
    })
}

/// Parses a Thermo Nicolet OMNIC `.spa` spectrum buffer.
#[pyfunction]
pub fn omnic_parse_spa(data: &[u8]) -> PyResult<OmnicSpectrum> {
    let s = py_result(omnic::parse_spa(data))?;
    Ok(OmnicSpectrum {
        wavenumbers: s.wavenumbers,
        intensities: s.intensities,
        title: s.title,
        y_unit: s.y_unit,
        acquired_epoch_secs: s.acquired_epoch_secs,
    })
}

/// Parsed plain-text UV-Vis tabular export.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct UvTextTable {
    pub wavelengths: Vec<f64>,
    pub scans: Vec<Vec<f64>>,
}

/// Parses a plain-text UV-Vis export buffer.
#[pyfunction]
pub fn uvtext_parse(data: &[u8]) -> PyResult<UvTextTable> {
    let t = py_result(uvtext::parse_uvtext(data))?;
    Ok(UvTextTable {
        wavelengths: t.wavelengths,
        scans: t.scans,
    })
}

/// One sample row from a Thermo Finnigan `.raw` run or `.sld` sequence.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct ThermoRow {
    pub sample_name: String,
    pub sample_id: String,
    pub data_file: String,
    pub vial: String,
    pub method: String,
}

/// Parsed Thermo Finnigan `.raw` or `.sld` header and rows.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct ThermoFile {
    pub first: String,
    pub second: String,
    pub logins: Vec<String>,
    pub rows: Vec<ThermoRow>,
    pub instrument: Option<String>,
    pub scans: u32,
    pub minutes: Option<(f64, f64)>,
}

impl From<thermo::ThermoFile> for ThermoFile {
    fn from(f: thermo::ThermoFile) -> Self {
        Self {
            first: f.first,
            second: f.second,
            logins: f.logins.into_iter().collect(),
            rows: f
                .rows
                .into_iter()
                .map(|r| ThermoRow {
                    sample_name: r.sample_name,
                    sample_id: r.sample_id,
                    data_file: r.data_file,
                    vial: r.vial,
                    method: r.method,
                })
                .collect(),
            instrument: f.instrument,
            scans: f.scans,
            minutes: f.minutes,
        }
    }
}

/// Decoded Thermo Finnigan scan index entries.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct ThermoScans {
    pub minutes: Vec<f64>,
    pub tic: Vec<f64>,
    pub offsets: Vec<u64>,
    pub data_offset: u64,
}

/// Physical byte spans of user-written login slots and UTF-16LE strings in a Thermo Finnigan file.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct ThermoTextSpans {
    pub logins: Vec<(usize, usize, usize)>,
    pub strings: Vec<(usize, usize)>,
}

#[pyfunction]
pub fn thermo_parse_run(data: &[u8]) -> PyResult<ThermoFile> {
    py_result(thermo::parse_run(data)).map(Into::into)
}

#[pyfunction]
pub fn thermo_parse_sequence(data: &[u8]) -> PyResult<ThermoFile> {
    py_result(thermo::parse_sequence(data)).map(Into::into)
}

#[pyfunction]
#[pyo3(signature = (data, is_sequence=false))]
pub fn thermo_text_spans(data: &[u8], is_sequence: bool) -> PyResult<ThermoTextSpans> {
    let spans = py_result(thermo::text_spans(data, is_sequence))?;
    Ok(ThermoTextSpans {
        logins: spans.logins,
        strings: spans.strings,
    })
}

#[pyfunction]
pub fn thermo_ms_header_offset(data: &[u8]) -> PyResult<Option<u64>> {
    py_result(thermo::ms_header_offset(data))
}

#[pyfunction]
pub fn thermo_ms_header_ranges(header_offset: u64) -> PyResult<((u64, usize), (u64, usize))> {
    py_result(thermo::ms_header_ranges(header_offset))
}

#[pyfunction]
pub fn thermo_decode_scan_range(bounds: &[u8], pointers: &[u8]) -> PyResult<(u64, usize, u64)> {
    py_result(thermo::decode_scan_range(bounds, pointers))
}

#[pyfunction]
pub fn thermo_decode_scan_index(entries: &[u8], data_offset: u64) -> PyResult<ThermoScans> {
    let s = py_result(thermo::decode_scan_index(entries, data_offset))?;
    Ok(ThermoScans {
        minutes: s.minutes,
        tic: s.tic,
        offsets: s.offsets,
        data_offset: s.data_offset,
    })
}

#[pyfunction]
pub fn thermo_decode_peak_record(record: &[u8]) -> PyResult<(u32, u32)> {
    py_result(thermo::decode_peak_record(record))
}

#[pyfunction]
pub fn thermo_decode_peak_pairs(bytes: &[u8]) -> PyResult<(Vec<f64>, Vec<f64>)> {
    py_result(thermo::decode_peak_pairs(bytes))
}

/// One UV channel from a Waters MassLynx `.raw` run.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct WatersUvChannel {
    pub wavelength: f64,
    pub values: Vec<f64>,
}

/// Full 2D MS + UV polarity slice parsed from a Waters MassLynx `.raw` bundle.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct WatersPolarityRun {
    pub polarity: String,
    pub acquired: Option<String>,
    pub instrument: Option<String>,
    pub inlet_method: Option<String>,
    pub ms_method: Option<String>,
    pub rt: Vec<f64>,
    pub mz: Vec<f64>,
    pub intensity: Vec<f64>,
    pub uv_rt: Option<Vec<f64>>,
    pub uv: Vec<WatersUvChannel>,
}

/// Parses every MS polarity and UV channel from a Waters MassLynx `.raw` file bundle.
#[pyfunction]
#[pyo3(signature = (files, bin_width=0.1))]
pub fn waters_read_polarities(
    files: &Bound<'_, PyAny>,
    bin_width: f64,
) -> PyResult<BTreeMap<String, WatersPolarityRun>> {
    let owned = extract_file_bundle(files)?;
    let refs: Vec<(&str, &[u8])> = owned
        .iter()
        .map(|(name, bytes)| (name.as_str(), bytes.as_slice()))
        .collect();
    let parsed = py_result(waters::read_polarities(&refs, bin_width))?;
    Ok(parsed
        .into_iter()
        .map(|(k, v)| {
            (
                k,
                WatersPolarityRun {
                    polarity: v.polarity,
                    acquired: v.acquired,
                    instrument: v.instrument,
                    inlet_method: v.inlet_method,
                    ms_method: v.ms_method,
                    rt: v.rt,
                    mz: v.mz,
                    intensity: v.intensity,
                    uv_rt: v.uv_rt,
                    uv: v
                        .uv
                        .into_iter()
                        .map(|c| WatersUvChannel {
                            wavelength: c.wavelength,
                            values: c.values,
                        })
                        .collect(),
                },
            )
        })
        .collect())
}

/// One MS1 profile scan parsed from Agilent MassHunter `MSScan.bin`.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct MassHunterScan {
    pub minutes: f64,
    pub polarity: String,
    pub tic: f64,
    pub offset: u64,
    pub length: u64,
    pub points: u32,
    pub calibration: [f64; 10],
    pub flags: Option<u64>,
}

/// Parses Agilent MassHunter `MSScan.bin` and calibration files into ordered MS1 profile scans.
#[pyfunction]
#[pyo3(signature = (msscan_bin, msscan_xsd, msts_xml=None, masscal_bin=None, default_masscal_xml=None))]
pub fn masshunter_profile_scans(
    msscan_bin: &[u8],
    msscan_xsd: &[u8],
    msts_xml: Option<&[u8]>,
    masscal_bin: Option<&[u8]>,
    default_masscal_xml: Option<&[u8]>,
) -> PyResult<Vec<MassHunterScan>> {
    let scans = py_result(masshunter::parse_profile_scans(
        msscan_bin,
        msscan_xsd,
        msts_xml,
        masscal_bin,
        default_masscal_xml,
    ))?;
    Ok(scans
        .into_iter()
        .map(|s| MassHunterScan {
            minutes: s.minutes,
            polarity: s.polarity,
            tic: s.tic,
            offset: s.offset,
            length: s.length,
            points: s.points,
            calibration: s.calibration,
            flags: s.flags,
        })
        .collect())
}

/// Decodes an Agilent MassHunter `MSProfile.bin` RLE segment into calibrated `(m/z, intensity)`.
#[pyfunction]
#[pyo3(signature = (segment, num_points, calibration, use_flags=None))]
pub fn masshunter_decode_spectrum(
    segment: &[u8],
    num_points: u32,
    calibration: [f64; 10],
    use_flags: Option<u64>,
) -> PyResult<(Vec<f64>, Vec<f32>)> {
    py_result(masshunter::decode_profile_spectrum(
        segment,
        num_points,
        &calibration,
        use_flags,
    ))
}

/// One file embedded in a MestReNova `.mnova` file table.
#[pyclass]
#[derive(Clone)]
pub struct MnovaEmbedded {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub entry: (usize, usize),
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub end: usize,
    #[pyo3(get)]
    pub path: (usize, usize),
    #[pyo3(get)]
    pub after: (usize, usize),
    text: Vec<u8>,
}

#[pymethods]
impl MnovaEmbedded {
    #[getter]
    fn text<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.text)
    }
}

/// One parameter row from a MestReNova `.mnova` parameter table.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct MnovaRow {
    pub span: (usize, usize),
    pub label: Option<String>,
    pub key: Option<String>,
    pub formula: Option<String>,
    pub value: (usize, usize),
}

/// One Bruker experiment embedded in a MestReNova `.mnova` document.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct MnovaSpectrum {
    pub dataset: String,
    pub expno: String,
    pub position: usize,
    pub acqus: AcqusParams,
}

/// One framed 1D JPEG 2000 processed spectrum extracted from a MestReNova `.mnova` document.
#[pyclass]
#[derive(Clone)]
pub struct MnovaJ2kFrame {
    #[pyo3(get)]
    pub at: usize,
    #[pyo3(get)]
    pub position: Option<usize>,
    #[pyo3(get)]
    pub nucleus: String,
    #[pyo3(get)]
    pub low_hz: f64,
    #[pyo3(get)]
    pub mhz: f64,
    #[pyo3(get)]
    pub width_hz: f64,
    #[pyo3(get)]
    pub low_val: f32,
    #[pyo3(get)]
    pub high_val: f32,
    #[pyo3(get)]
    pub depth: u8,
    codestream: Vec<u8>,
}

#[pymethods]
impl MnovaJ2kFrame {
    #[getter]
    fn codestream<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.codestream)
    }
}

/// Parses every embedded file table in a MestReNova `.mnova` document.
#[pyfunction]
pub fn mnova_embedded(content: &[u8]) -> PyResult<Vec<MnovaEmbedded>> {
    let items = py_result(mnova::parse_embedded(content))?;
    Ok(items
        .into_iter()
        .map(|e| MnovaEmbedded {
            name: e.name,
            entry: e.entry,
            start: e.start,
            end: e.end,
            path: e.path,
            after: e.after,
            text: e.text,
        })
        .collect())
}

/// Parses every parameter table row in a MestReNova `.mnova` document.
#[pyfunction]
pub fn mnova_rows(content: &[u8]) -> PyResult<Vec<MnovaRow>> {
    let items = py_result(mnova::parse_rows(content))?;
    Ok(items
        .into_iter()
        .map(|r| MnovaRow {
            span: r.span,
            label: r.label,
            key: r.key,
            formula: r.formula,
            value: r.value,
        })
        .collect())
}

/// Finds every measurement byte span (JPEG 2000 spectra and FIDs) in a MestReNova `.mnova` document.
#[pyfunction]
pub fn mnova_measurements(content: &[u8]) -> Vec<(usize, usize)> {
    mnova::parse_measurements(content)
}

/// Extracts every embedded Bruker `acqus` experiment in a MestReNova `.mnova` document.
#[pyfunction]
pub fn mnova_spectra(content: &[u8]) -> Vec<MnovaSpectrum> {
    mnova::parse_spectra(content)
        .into_iter()
        .map(|s| MnovaSpectrum {
            dataset: s.dataset,
            expno: s.expno,
            position: s.position,
            acqus: s.acqus.into(),
        })
        .collect()
}

/// Extracts every framed 1D JPEG 2000 spectrum codestream and axis in a MestReNova `.mnova` document.
#[pyfunction]
pub fn mnova_j2k_frames(content: &[u8]) -> Vec<MnovaJ2kFrame> {
    mnova::parse_j2k_frames(content)
        .into_iter()
        .map(|f| MnovaJ2kFrame {
            at: f.at,
            position: f.position,
            nucleus: f.nucleus,
            low_hz: f.low_hz,
            mhz: f.mhz,
            width_hz: f.width_hz,
            low_val: f.low_val,
            high_val: f.high_val,
            depth: f.depth,
            codestream: f.codestream,
        })
        .collect()
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<AndiDataset>()?;
    m.add_class::<OmnicHeader>()?;
    m.add_class::<OmnicSpectrum>()?;
    m.add_class::<UvTextTable>()?;
    m.add_class::<ThermoRow>()?;
    m.add_class::<ThermoFile>()?;
    m.add_class::<ThermoScans>()?;
    m.add_class::<ThermoTextSpans>()?;
    m.add("THERMO_MAGIC", PyBytes::new(m.py(), thermo::MAGIC))?;
    m.add("THERMO_PEAK_RECORD_LEN", thermo::PEAK_RECORD_LEN)?;
    m.add_class::<WatersUvChannel>()?;
    m.add_class::<WatersPolarityRun>()?;
    m.add_class::<MassHunterScan>()?;
    m.add_class::<MnovaEmbedded>()?;
    m.add_class::<MnovaRow>()?;
    m.add_class::<MnovaSpectrum>()?;
    m.add_class::<MnovaJ2kFrame>()?;
    m.add_function(wrap_pyfunction!(andi_parse_cdf, m)?)?;
    m.add_function(wrap_pyfunction!(chromeleon_signal_name, m)?)?;
    m.add_function(wrap_pyfunction!(chromeleon_decode_signal, m)?)?;
    m.add_function(wrap_pyfunction!(omnic_parse_header, m)?)?;
    m.add_function(wrap_pyfunction!(omnic_parse_spa, m)?)?;
    m.add_function(wrap_pyfunction!(uvtext_parse, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_parse_run, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_parse_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_text_spans, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_ms_header_offset, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_ms_header_ranges, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_decode_scan_range, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_decode_scan_index, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_decode_peak_record, m)?)?;
    m.add_function(wrap_pyfunction!(thermo_decode_peak_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(waters_read_polarities, m)?)?;
    m.add_function(wrap_pyfunction!(masshunter_profile_scans, m)?)?;
    m.add_function(wrap_pyfunction!(masshunter_decode_spectrum, m)?)?;
    m.add_function(wrap_pyfunction!(mnova_embedded, m)?)?;
    m.add_function(wrap_pyfunction!(mnova_rows, m)?)?;
    m.add_function(wrap_pyfunction!(mnova_measurements, m)?)?;
    m.add_function(wrap_pyfunction!(mnova_spectra, m)?)?;
    m.add_function(wrap_pyfunction!(mnova_j2k_frames, m)?)?;
    Ok(())
}
