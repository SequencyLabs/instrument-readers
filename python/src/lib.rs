//! PyO3 bindings exposing `sequency-instruments` binary parsers to Python.

mod formats;

use std::collections::BTreeMap;
use std::path::Path;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use sequency_instruments::{
    GradientStep as RustGradientStep, InstrumentMethod as RustInstrumentMethod,
    InstrumentReading as RustInstrumentReading, Peak as RustPeak, Trace as RustTrace, bruker,
    detect as rust_detect, ole, openlab, read as rust_read, read_files as rust_read_files,
    shimadzu,
};

type Extent = (usize, usize);

pub(crate) fn py_result<T>(res: sequency_instruments::Result<T>) -> PyResult<T> {
    res.map_err(|err| PyValueError::new_err(err.to_string()))
}

/// One step in a chromatographic mobile-phase gradient program.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct GradientStep {
    pub time_min: f64,
    pub b_pct: f64,
    pub c_pct: f64,
    pub d_pct: f64,
}

#[pymethods]
impl GradientStep {
    /// Implied mobile phase line A percentage (`100 - B - C - D`, clamped at zero).
    #[must_use]
    pub fn a_pct(&self) -> f64 {
        RustGradientStep {
            time_min: self.time_min,
            b_pct: self.b_pct,
            c_pct: self.c_pct,
            d_pct: self.d_pct,
        }
        .a_pct()
    }
}

impl From<RustGradientStep> for GradientStep {
    fn from(step: RustGradientStep) -> Self {
        let RustGradientStep {
            time_min,
            b_pct,
            c_pct,
            d_pct,
        } = step;
        Self {
            time_min,
            b_pct,
            c_pct,
            d_pct,
        }
    }
}

/// Unified acquisition and processing method parameters across instrument techniques.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct InstrumentMethod {
    pub name: Option<String>,
    pub column: Option<String>,
    pub solvents: Vec<String>,
    pub gradient: Vec<GradientStep>,
    pub flow_rate_ml_min: Option<f64>,
    pub column_temp_c: Option<f64>,
    pub injection_volume_ul: Option<f64>,
    pub wavelength_nm: Option<f64>,
    pub nucleus: Option<String>,
    pub frequency_mhz: Option<f64>,
    pub pulse_program: Option<String>,
    pub scans: Option<u32>,
    pub relaxation_delay_s: Option<f64>,
}

impl From<RustInstrumentMethod> for InstrumentMethod {
    fn from(method: RustInstrumentMethod) -> Self {
        let RustInstrumentMethod {
            name,
            column,
            solvents,
            gradient,
            flow_rate_ml_min,
            column_temp_c,
            injection_volume_ul,
            wavelength_nm,
            nucleus,
            frequency_mhz,
            pulse_program,
            scans,
            relaxation_delay_s,
        } = method;
        Self {
            name,
            column,
            solvents,
            gradient: gradient.into_iter().map(Into::into).collect(),
            flow_rate_ml_min,
            column_temp_c,
            injection_volume_ul,
            wavelength_nm,
            nucleus,
            frequency_mhz,
            pulse_program,
            scans,
            relaxation_delay_s,
        }
    }
}

/// Integrated chromatographic peak from a vendor peak table.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct Peak {
    pub rt_min: Option<f64>,
    pub area: Option<f64>,
    pub area_pct: Option<f64>,
    pub height: Option<f64>,
    pub compound: Option<String>,
}

impl From<RustPeak> for Peak {
    fn from(peak: RustPeak) -> Self {
        let RustPeak {
            rt_min,
            area,
            area_pct,
            height,
            compound,
        } = peak;
        Self {
            rt_min,
            area,
            area_pct,
            height,
            compound,
        }
    }
}

/// Continuous 1D chromatogram or spectrum trace ready for plotting.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct Trace {
    pub label: String,
    pub x_unit: String,
    pub y_unit: String,
    pub points: Vec<(f64, f64)>,
}

impl Trace {
    fn to_rust(&self) -> RustTrace {
        RustTrace {
            label: self.label.clone(),
            x_unit: self.x_unit.clone(),
            y_unit: self.y_unit.clone(),
            points: self.points.clone(),
        }
    }
}

#[pymethods]
impl Trace {
    /// Downsamples points to at most `threshold` points using Largest-Triangle-Three-Buckets.
    #[must_use]
    pub fn downsample_lttb(&self, threshold: usize) -> Self {
        self.to_rust().downsample_lttb(threshold).into()
    }

    /// Downsamples points into `bins` horizontal buckets keeping each bucket's min and max.
    #[must_use]
    pub fn downsample_minmax(&self, bins: usize) -> Self {
        self.to_rust().downsample_minmax(bins).into()
    }
}

impl From<RustTrace> for Trace {
    fn from(trace: RustTrace) -> Self {
        let RustTrace {
            label,
            x_unit,
            y_unit,
            points,
        } = trace;
        Self {
            label,
            x_unit,
            y_unit,
            points,
        }
    }
}

/// Unified reading parsed from an analytical instrument run.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct InstrumentReading {
    pub technique: String,
    pub vendor_format: String,
    pub acquired_at: Option<String>,
    pub sample_name: Option<String>,
    pub instrument: Option<String>,
    pub method: InstrumentMethod,
    pub peaks: Vec<Peak>,
    pub traces: Vec<Trace>,
}

#[pymethods]
impl InstrumentReading {
    /// Returns the peak with the largest integrated area, if any peak records an area.
    #[must_use]
    pub fn main_peak(&self) -> Option<Peak> {
        self.peaks
            .iter()
            .filter(|p| p.area.is_some())
            .max_by(|a, b| {
                let lhs = a.area.unwrap_or(f64::NEG_INFINITY);
                let rhs = b.area.unwrap_or(f64::NEG_INFINITY);
                lhs.total_cmp(&rhs)
            })
            .cloned()
    }
}

impl From<RustInstrumentReading> for InstrumentReading {
    fn from(reading: RustInstrumentReading) -> Self {
        let RustInstrumentReading {
            technique,
            vendor_format,
            acquired_at,
            sample_name,
            instrument,
            method,
            peaks,
            traces,
        } = reading;
        Self {
            technique: technique.as_str().to_owned(),
            vendor_format: vendor_format.as_str().to_owned(),
            acquired_at,
            sample_name,
            instrument,
            method: method.into(),
            peaks: peaks.into_iter().map(Into::into).collect(),
            traces: traces.into_iter().map(Into::into).collect(),
        }
    }
}

/// Run summary read from a Shimadzu `.lcd` file's `File Property`, method and peak table streams.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct ShimadzuRun {
    pub injected: Option<String>,
    pub sample_name: Option<String>,
    pub data_file: Option<String>,
    pub method_file: Option<String>,
    pub instrument: String,
    pub method: InstrumentMethod,
    pub peaks: Vec<Peak>,
}

impl From<shimadzu::ShimadzuRun> for ShimadzuRun {
    fn from(run: shimadzu::ShimadzuRun) -> Self {
        let shimadzu::ShimadzuRun {
            injected,
            sample_name,
            data_file,
            method_file,
            instrument,
            method,
            peaks,
        } = run;
        Self {
            injected,
            sample_name,
            data_file,
            method_file,
            instrument,
            method: method.into(),
            peaks: peaks.into_iter().map(Into::into).collect(),
        }
    }
}

/// One stream inside an OLE2 compound file and its byte extents.
#[pyclass]
#[derive(Clone)]
pub struct OleStream {
    #[pyo3(get)]
    pub name: String,
    data: Vec<u8>,
    #[pyo3(get)]
    pub extents: Vec<Extent>,
}

#[pymethods]
impl OleStream {
    /// Raw bytes of the OLE2 stream.
    #[getter]
    fn data<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.data)
    }
}

/// Parsed OLE2 compound file layout with all streams and unallocated byte extents.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct OleLayout {
    pub streams: Vec<OleStream>,
    pub unallocated: Vec<Extent>,
}

/// Parsed Bruker TopSpin `acqus` parameters.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct AcqusParams {
    pub acquired_at: Option<String>,
    pub nucleus: String,
    pub pulse_program: String,
    pub instrument: String,
    pub solvent: String,
    pub scans: u32,
    pub relaxation_delay_s: f64,
    pub sfo1_mhz: Option<f64>,
    pub bf1_mhz: Option<f64>,
    pub sw_h_hz: Option<f64>,
    pub data_path: Option<String>,
}

impl From<bruker::AcqusParams> for AcqusParams {
    fn from(params: bruker::AcqusParams) -> Self {
        let bruker::AcqusParams {
            acquired_at,
            nucleus,
            pulse_program,
            instrument,
            solvent,
            scans,
            relaxation_delay_s,
            sfo1_mhz,
            bf1_mhz,
            sw_h_hz,
            data_path,
        } = params;
        Self {
            acquired_at,
            nucleus,
            pulse_program,
            instrument,
            solvent,
            scans,
            relaxation_delay_s,
            sfo1_mhz,
            bf1_mhz,
            sw_h_hz,
            data_path,
        }
    }
}

/// One injection entry parsed from an Agilent OpenLab `.acaml` result set.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct OpenLabInjection {
    pub id: Option<String>,
    pub data_file: String,
    pub sample_name: String,
    pub sample_type: String,
    pub vial: String,
    pub instrument: Option<String>,
    pub acquired_utc: Option<String>,
    pub method: String,
    pub processing_method: Option<String>,
    pub status: Option<String>,
    pub result_file: Option<String>,
    pub injection_volume_ul: Option<f64>,
}

impl From<openlab::OpenLabInjection> for OpenLabInjection {
    fn from(inj: openlab::OpenLabInjection) -> Self {
        let openlab::OpenLabInjection {
            id,
            data_file,
            sample_name,
            sample_type,
            vial,
            instrument,
            acquired_utc,
            method,
            processing_method,
            status,
            result_file,
            injection_volume_ul,
        } = inj;
        Self {
            id,
            data_file,
            sample_name,
            sample_type,
            vial,
            instrument,
            acquired_utc,
            method,
            processing_method,
            status,
            result_file,
            injection_volume_ul,
        }
    }
}

/// Parsed Agilent OpenLab `.acaml` result set.
#[pyclass(get_all)]
#[derive(Clone)]
pub struct AcamlResultSet {
    pub injections: Vec<OpenLabInjection>,
    pub logins: Vec<String>,
}

/// Detects the vendor format of an instrument file from magic bytes and optional filename.
#[pyfunction]
#[pyo3(signature = (data, filename=None))]
fn detect(data: &[u8], filename: Option<&str>) -> String {
    rust_detect(filename, data).as_str().to_owned()
}

/// Automatically detects and parses a single instrument file into an `InstrumentReading`.
#[pyfunction]
#[pyo3(signature = (data, filename=None))]
fn read(data: &[u8], filename: Option<&str>) -> PyResult<InstrumentReading> {
    py_result(rust_read(filename, data)).map(Into::into)
}

/// Automatically detects and parses a multi-file instrument artifact bundle in memory.
#[pyfunction]
fn read_files(files: &Bound<'_, PyAny>) -> PyResult<InstrumentReading> {
    let owned = formats::extract_file_bundle(files)?;
    let refs: Vec<(&str, &[u8])> = owned
        .iter()
        .map(|(name, bytes)| (name.as_str(), bytes.as_slice()))
        .collect();
    py_result(rust_read_files(&refs)).map(Into::into)
}

/// Parses an OLE2 compound file and returns all streams with their byte extents.
#[pyfunction]
fn ole_layout(data: &[u8]) -> PyResult<OleLayout> {
    let parsed = py_result(ole::OleFile::parse(data))?;
    let mut streams = Vec::with_capacity(parsed.streams().len());
    for info in parsed.streams() {
        let bytes = parsed.stream(&info.name).unwrap_or_default().into_owned();
        streams.push(OleStream {
            name: info.name.clone(),
            data: bytes,
            extents: info.extents.clone(),
        });
    }
    Ok(OleLayout {
        streams,
        unallocated: parsed.unallocated().to_vec(),
    })
}

/// Maps a `[start, end)` byte range inside a stream onto file byte extents.
#[pyfunction]
fn ole_place(extents: &Bound<'_, PyAny>, start: usize, end: usize) -> PyResult<Vec<Extent>> {
    let extents: Vec<Extent> = extents.extract()?;
    Ok(ole::place(&extents, start, end))
}

/// Parses a Shimadzu `.lcd` file into a unified `InstrumentReading`.
#[pyfunction]
fn shimadzu_read(data: &[u8]) -> PyResult<InstrumentReading> {
    py_result(shimadzu::read_lcd(data)).map(Into::into)
}

/// Reads a Shimadzu `.lcd` file's run summary without decoding its chromatograms, so a damaged trace does not hide
/// the sample, files, injection time and peaks.
#[pyfunction]
fn shimadzu_read_run(data: &[u8]) -> PyResult<ShimadzuRun> {
    py_result(shimadzu::read_run(data)).map(Into::into)
}

/// Parses a Shimadzu `.lcm` method file into a unified `InstrumentReading`.
#[pyfunction]
fn shimadzu_read_lcm(data: &[u8]) -> PyResult<InstrumentReading> {
    py_result(shimadzu::read_lcm(data)).map(Into::into)
}

/// Parses a Shimadzu `File Property` stream into key-value pairs.
#[pyfunction]
fn shimadzu_properties(stream: &[u8]) -> BTreeMap<String, String> {
    shimadzu::parse_properties(stream)
}

/// Converts Shimadzu Windows FILETIME low/high decimal strings into an ISO-8601 UTC string.
#[pyfunction]
#[pyo3(signature = (low=None, high=None))]
fn shimadzu_filetime(low: Option<&str>, high: Option<&str>) -> Option<String> {
    shimadzu::parse_filetime(low, high)
}

/// Parses a Shimadzu `PT-*` peak table stream.
#[pyfunction]
fn shimadzu_peak_table(stream: &[u8]) -> Vec<Peak> {
    shimadzu::parse_peak_table(stream)
        .into_iter()
        .map(Into::into)
        .collect()
}

/// Decodes a delta-encoded Shimadzu chromatogram sub-block.
#[pyfunction]
fn shimadzu_deltas(sub_block: &[u8]) -> PyResult<Vec<i64>> {
    py_result(shimadzu::decode_deltas(sub_block))
}

/// Parses a Shimadzu chromatogram stream into `(times_min, intensities)`.
#[pyfunction]
fn shimadzu_trace(stream: &[u8]) -> PyResult<(Vec<f64>, Vec<i64>)> {
    py_result(shimadzu::parse_trace(stream))
}

/// Extracts all chromatogram traces from a Shimadzu `.lcd` file.
#[pyfunction]
fn shimadzu_chromatograms(data: &[u8]) -> PyResult<Vec<Trace>> {
    let traces = py_result(shimadzu::chromatograms(data))?;
    Ok(traces.into_iter().map(Into::into).collect())
}

/// Parses Bruker TopSpin `acqus` JCAMP-DX parameters.
#[pyfunction]
fn bruker_parse_acqus(text: &str) -> AcqusParams {
    bruker::parse_acqus_fields(text).into()
}

/// Parses a Bruker `acqus` byte buffer into an `InstrumentReading`.
#[pyfunction]
fn bruker_read_acqus(data: &[u8]) -> PyResult<InstrumentReading> {
    py_result(bruker::read_acqus(data)).map(Into::into)
}

/// Parses Bruker processed `procs` and `1r` files (with optional `acqus`) into an `InstrumentReading`.
#[pyfunction]
#[pyo3(signature = (procs, data_1r, acqus=None))]
fn bruker_read_pdata(
    procs: &[u8],
    data_1r: &[u8],
    acqus: Option<&[u8]>,
) -> PyResult<InstrumentReading> {
    py_result(bruker::read_pdata(acqus, procs, data_1r)).map(Into::into)
}

/// Reads a Bruker TopSpin experiment directory from disk into an `InstrumentReading`.
#[pyfunction]
fn bruker_read_dir(exp_dir: &str) -> PyResult<InstrumentReading> {
    py_result(bruker::read_dir(Path::new(exp_dir))).map(Into::into)
}

/// Parses an Agilent OpenLab `.acaml` XML document into injections and operator logins.
#[pyfunction]
fn openlab_parse_acaml(acaml: &[u8]) -> PyResult<AcamlResultSet> {
    let openlab::AcamlResultSet { injections, logins } = py_result(openlab::parse_acaml(acaml))?;
    Ok(AcamlResultSet {
        injections: injections.into_iter().map(Into::into).collect(),
        logins: logins.into_iter().collect(),
    })
}

/// Reads an Agilent OpenLab `.acaml` XML document into an `InstrumentReading`.
#[pyfunction]
fn openlab_read_acaml(acaml: &[u8]) -> PyResult<InstrumentReading> {
    py_result(openlab::read_acaml(acaml)).map(Into::into)
}

/// Reads an Agilent OpenLab `.dx` zipped signal archive into an `InstrumentReading`.
#[pyfunction]
fn openlab_read_dx(dx: &[u8]) -> PyResult<InstrumentReading> {
    py_result(openlab::read_dx(dx)).map(Into::into)
}

/// Reads an Agilent OpenLab `.rx` zipped integration result archive into an `InstrumentReading`.
#[pyfunction]
fn openlab_read_rx(rx: &[u8]) -> PyResult<InstrumentReading> {
    py_result(openlab::read_rx(rx)).map(Into::into)
}

/// Extracts integrated peaks from an Agilent OpenLab `.rx` archive.
#[pyfunction]
fn openlab_parse_rx_peaks(rx: &[u8]) -> PyResult<Vec<Peak>> {
    let peaks = py_result(openlab::parse_rx_peaks(rx))?;
    Ok(peaks.into_iter().map(Into::into).collect())
}

#[pymodule(name = "sequency_instruments")]
fn sequency_instruments_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GradientStep>()?;
    m.add_class::<InstrumentMethod>()?;
    m.add_class::<Peak>()?;
    m.add_class::<Trace>()?;
    m.add_class::<InstrumentReading>()?;
    m.add_class::<ShimadzuRun>()?;
    m.add_class::<OleStream>()?;
    m.add_class::<OleLayout>()?;
    m.add_class::<AcqusParams>()?;
    m.add_class::<OpenLabInjection>()?;
    m.add_class::<AcamlResultSet>()?;
    m.add_function(wrap_pyfunction!(detect, m)?)?;
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(read_files, m)?)?;
    m.add_function(wrap_pyfunction!(ole_layout, m)?)?;
    m.add_function(wrap_pyfunction!(ole_place, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_read, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_read_run, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_read_lcm, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_properties, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_filetime, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_peak_table, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_deltas, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_trace, m)?)?;
    m.add_function(wrap_pyfunction!(shimadzu_chromatograms, m)?)?;
    m.add_function(wrap_pyfunction!(bruker_parse_acqus, m)?)?;
    m.add_function(wrap_pyfunction!(bruker_read_acqus, m)?)?;
    m.add_function(wrap_pyfunction!(bruker_read_pdata, m)?)?;
    m.add_function(wrap_pyfunction!(bruker_read_dir, m)?)?;
    m.add_function(wrap_pyfunction!(openlab_parse_acaml, m)?)?;
    m.add_function(wrap_pyfunction!(openlab_read_acaml, m)?)?;
    m.add_function(wrap_pyfunction!(openlab_read_dx, m)?)?;
    m.add_function(wrap_pyfunction!(openlab_read_rx, m)?)?;
    m.add_function(wrap_pyfunction!(openlab_parse_rx_peaks, m)?)?;
    formats::register(m)?;
    Ok(())
}
