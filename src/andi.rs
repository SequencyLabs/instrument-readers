//! Pure-Rust ANDI Chromatography (ASTM E1947) reader built on [`crate::netcdf`].
//!
//! Extracts chromatography metadata, detector traces, and integrated peak tables from classic
//! (`CDF\x01`) and 64-bit offset (`CDF\x02`) netCDF-3 files.

use serde::{Deserialize, Serialize};

use crate::binary::round_digits;
use crate::error::Result;
pub use crate::netcdf::{NcAttrValue, NcDimension, NcType, NcVariable, NetCdfFile, NetCdfVersion};
use crate::netcdf::{custom_fill_value, is_valid_number};
use crate::types::{InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat};

/// Integrated peak extracted from an ANDI Chromatography `.cdf` peak table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AndiPeak {
    /// Retention time at the peak apex in minutes.
    pub rt_min: Option<f64>,
    /// Integrated peak area in detector units.
    pub area: Option<f64>,
    /// Relative area percentage across integrated peaks.
    pub area_pct: Option<f64>,
    /// Peak height above baseline in detector units.
    pub height: Option<f64>,
    /// Peak width in minutes, when present.
    pub width_min: Option<f64>,
    /// Peak start time in minutes, when present.
    pub start_time_min: Option<f64>,
    /// Peak end time in minutes, when present.
    pub end_time_min: Option<f64>,
    /// Compound name assigned to this peak, when present.
    pub compound: Option<String>,
}

impl AndiPeak {
    /// Converts this ANDI peak into the crate's unified [`Peak`].
    #[must_use]
    pub fn into_peak(self) -> Peak {
        Peak {
            rt_min: self.rt_min,
            area: self.area,
            area_pct: self.area_pct,
            height: self.height,
            compound: self.compound,
        }
    }
}

/// Full parsed ANDI Chromatography netCDF-3 dataset.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AndiDataset {
    /// NetCDF-3 binary format version (`CDF\x01` or `CDF\x02`).
    pub version: NetCdfVersion,
    /// Sample name from the `sample_name` global attribute.
    pub sample_name: Option<String>,
    /// Operator name from `operator_name`.
    pub operator_name: Option<String>,
    /// Detection method name from `detection_method_name`.
    pub detection_method_name: Option<String>,
    /// Raw `injection_date_time_stamp` global attribute, when present.
    pub injection_date_time_stamp: Option<String>,
    /// Normalized acquisition timestamp from `dataset_date_time_stamp` or `injection_date_time_stamp`.
    pub acquired_at: Option<String>,
    /// Detector name from `detector_name`.
    pub detector_name: Option<String>,
    /// Experiment title from `experiment_title`.
    pub experiment_title: Option<String>,
    /// Source file reference from `source_file_reference`.
    pub source_file_reference: Option<String>,
    /// Detector signal unit from `detector_unit`.
    pub detector_unit: Option<String>,
    /// Acquisition delay time in seconds (`actual_delay_time`).
    pub delay_time_s: Option<f64>,
    /// Total run duration in seconds (`actual_run_time_length`).
    pub run_time_length_s: Option<f64>,
    /// Sampling interval in seconds (`actual_sampling_interval`).
    pub sampling_interval_s: Option<f64>,
    /// Integrated peaks from the ANDI peak table.
    pub peaks: Vec<AndiPeak>,
    /// Extracted chromatogram traces.
    pub traces: Vec<Trace>,
}

/// Parses an ANDI Chromatography netCDF-3 `.cdf` buffer into an [`AndiDataset`].
pub fn parse_cdf(data: &[u8]) -> Result<AndiDataset> {
    let nc = NetCdfFile::parse(data)?;

    let sample_name = nc.global_str("sample_name").map(ToOwned::to_owned);
    let operator_name = nc.global_str("operator_name").map(ToOwned::to_owned);
    let detection_method_name = nc
        .global_str("detection_method_name")
        .map(ToOwned::to_owned);
    let injection_date_time_stamp = nc
        .global_str("injection_date_time_stamp")
        .map(ToOwned::to_owned);
    let acquired_at = nc
        .global_str("dataset_date_time_stamp")
        .or_else(|| nc.global_str("injection_date_time_stamp"))
        .map(normalize_andi_timestamp);
    let detector_name = nc
        .global_str("detector_name")
        .map(ToOwned::to_owned)
        .or(nc.read_var_string("detector_name")?);
    let experiment_title = nc.global_str("experiment_title").map(ToOwned::to_owned);
    let source_file_reference = nc
        .global_str("source_file_reference")
        .map(ToOwned::to_owned);
    let detector_unit = nc
        .read_var_string("detector_unit")?
        .or_else(|| nc.global_str("detector_unit").map(ToOwned::to_owned))
        .or_else(|| {
            nc.variable("ordinate_values")
                .and_then(|v| v.attrs.get("units"))
                .and_then(NcAttrValue::as_str)
                .map(ToOwned::to_owned)
        });

    let delay_time_s = read_time_scalar_seconds(&nc, "actual_delay_time")?;
    let run_time_length_s = read_time_scalar_seconds(&nc, "actual_run_time_length")?;
    let sampling_interval_s = read_time_scalar_seconds(&nc, "actual_sampling_interval")?;

    let mut traces = Vec::new();
    if let Some(trace) = build_chromatogram_trace(
        &nc,
        detector_name.as_deref(),
        experiment_title.as_deref(),
        detector_unit.as_deref(),
        delay_time_s,
        run_time_length_s,
        sampling_interval_s,
    )? {
        traces.push(trace);
    }

    let peaks = parse_peak_table(&nc)?;

    Ok(AndiDataset {
        version: nc.version(),
        sample_name,
        operator_name,
        detection_method_name,
        injection_date_time_stamp,
        acquired_at,
        detector_name,
        experiment_title,
        source_file_reference,
        detector_unit,
        delay_time_s,
        run_time_length_s,
        sampling_interval_s,
        peaks,
        traces,
    })
}

/// Parses an ANDI Chromatography netCDF-3 `.cdf` file into a unified [`InstrumentReading`].
pub fn read_cdf(data: &[u8]) -> Result<InstrumentReading> {
    let dataset = parse_cdf(data)?;
    let method_name = dataset
        .experiment_title
        .clone()
        .or_else(|| dataset.source_file_reference.clone());

    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::AndiCdf,
        acquired_at: dataset.acquired_at,
        sample_name: dataset.sample_name.or(dataset.experiment_title),
        instrument: dataset.detector_name,
        method: InstrumentMethod {
            name: method_name,
            ..InstrumentMethod::default()
        },
        peaks: dataset.peaks.into_iter().map(AndiPeak::into_peak).collect(),
        traces: dataset.traces,
    })
}

fn build_chromatogram_trace(
    nc: &NetCdfFile<'_>,
    detector_name: Option<&str>,
    experiment_title: Option<&str>,
    detector_unit: Option<&str>,
    delay_time_s: Option<f64>,
    run_time_length_s: Option<f64>,
    sampling_interval_s: Option<f64>,
) -> Result<Option<Trace>> {
    let Some(ord_var) = nc.variable("ordinate_values") else {
        return Ok(None);
    };
    let Some(raw_y) = nc.read_var_f64("ordinate_values")? else {
        return Ok(None);
    };
    if raw_y.is_empty() {
        return Ok(None);
    }

    let fill = custom_fill_value(ord_var);
    let scale = ord_var
        .attrs
        .get("scale_factor")
        .and_then(NcAttrValue::as_f64)
        .unwrap_or(1.0);
    let offset = ord_var
        .attrs
        .get("add_offset")
        .and_then(NcAttrValue::as_f64)
        .unwrap_or(0.0);

    let explicit_times_min = read_optional_time_vec_minutes(nc, "raw_data_retention")?
        .or(read_optional_time_vec_minutes(nc, "retention_time")?);

    let delay_min = delay_time_s.unwrap_or(0.0) / 60.0;
    let n = raw_y.len();
    let step_min = if let Some(interval_s) = sampling_interval_s.filter(|&v| v > 0.0) {
        interval_s / 60.0
    } else if let Some(run_s) = run_time_length_s.filter(|&v| v > 0.0) {
        if n > 1 {
            (run_s / 60.0) / ((n - 1) as f64)
        } else {
            0.0
        }
    } else {
        1.0 / 60.0
    };

    let mut points = Vec::with_capacity(n);
    for (idx, &raw_val) in raw_y.iter().enumerate() {
        if !is_valid_number(raw_val, fill) {
            continue;
        }
        let y = raw_val * scale + offset;
        if !y.is_finite() {
            continue;
        }
        let x = explicit_times_min
            .as_ref()
            .and_then(|times| times.get(idx).copied().flatten())
            .unwrap_or_else(|| round_digits(delay_min + (idx as f64) * step_min, 6));
        points.push((x, y));
    }

    if points.is_empty() {
        return Ok(None);
    }

    let label = detector_name
        .or(experiment_title)
        .unwrap_or("Chromatogram")
        .to_owned();
    let y_unit = detector_unit.unwrap_or("intensity").to_owned();

    Ok(Some(Trace {
        label,
        x_unit: "min".to_owned(),
        y_unit,
        points,
    }))
}

fn parse_peak_table(nc: &NetCdfFile<'_>) -> Result<Vec<AndiPeak>> {
    let rts = read_optional_time_vec_minutes(nc, "peak_retention_time")?.unwrap_or_default();
    let areas = read_optional_num_vec(nc, "peak_area")?.unwrap_or_default();
    let area_pcts = read_optional_num_vec(nc, "peak_area_percent")?.unwrap_or_default();
    let heights = read_optional_num_vec(nc, "peak_height")?.unwrap_or_default();
    let widths = read_optional_time_vec_minutes(nc, "peak_width")?.unwrap_or_default();
    let starts = read_optional_time_vec_minutes(nc, "peak_start_time")?.unwrap_or_default();
    let ends = read_optional_time_vec_minutes(nc, "peak_end_time")?.unwrap_or_default();
    let names = nc.read_var_strings("peak_name")?.unwrap_or_default();

    let count = [
        rts.len(),
        areas.len(),
        area_pcts.len(),
        heights.len(),
        widths.len(),
        starts.len(),
        ends.len(),
        names.len(),
    ]
    .into_iter()
    .max()
    .unwrap_or(0);

    let mut peaks = Vec::with_capacity(count);
    for i in 0..count {
        let compound = names
            .get(i)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(ToOwned::to_owned);
        let peak = AndiPeak {
            rt_min: rts.get(i).copied().flatten(),
            area: areas.get(i).copied().flatten(),
            area_pct: area_pcts.get(i).copied().flatten(),
            height: heights.get(i).copied().flatten(),
            width_min: widths.get(i).copied().flatten(),
            start_time_min: starts.get(i).copied().flatten(),
            end_time_min: ends.get(i).copied().flatten(),
            compound,
        };
        if peak.rt_min.is_some()
            || peak.area.is_some()
            || peak.area_pct.is_some()
            || peak.height.is_some()
            || peak.width_min.is_some()
            || peak.start_time_min.is_some()
            || peak.end_time_min.is_some()
            || peak.compound.is_some()
        {
            peaks.push(peak);
        }
    }
    Ok(peaks)
}

fn read_time_scalar_seconds(nc: &NetCdfFile<'_>, name: &str) -> Result<Option<f64>> {
    let Some(val) = nc.read_var_scalar_f64(name)? else {
        return Ok(None);
    };
    let factor = if variable_units_are_minutes(nc, name) {
        60.0
    } else {
        1.0
    };
    Ok(Some(val * factor))
}

fn read_optional_num_vec(nc: &NetCdfFile<'_>, name: &str) -> Result<Option<Vec<Option<f64>>>> {
    let Some(var) = nc.variable(name) else {
        return Ok(None);
    };
    let fill = custom_fill_value(var);
    let Some(raw) = nc.read_var_f64(name)? else {
        return Ok(None);
    };
    Ok(Some(
        raw.into_iter()
            .map(|v| is_valid_number(v, fill).then_some(v))
            .collect(),
    ))
}

fn read_optional_time_vec_minutes(
    nc: &NetCdfFile<'_>,
    name: &str,
) -> Result<Option<Vec<Option<f64>>>> {
    let Some(values) = read_optional_num_vec(nc, name)? else {
        return Ok(None);
    };
    let divisor = if variable_units_are_minutes(nc, name) {
        1.0
    } else {
        60.0
    };
    Ok(Some(
        values
            .into_iter()
            .map(|opt| opt.map(|v| round_digits(v / divisor, 6)))
            .collect(),
    ))
}

fn variable_units_are_minutes(nc: &NetCdfFile<'_>, name: &str) -> bool {
    nc.variable(name)
        .and_then(|v| v.attrs.get("units"))
        .and_then(NcAttrValue::as_str)
        .is_some_and(|u| u.to_ascii_lowercase().starts_with("min"))
}

/// Normalizes an ASTM E1947 compact timestamp (`YYYYMMDDhhmmss[±hhmm|Z]`) to RFC 3339 when possible.
fn normalize_andi_timestamp(raw: &str) -> String {
    let trimmed = raw.trim();
    let bytes = trimmed.as_bytes();
    if bytes.len() >= 14 && bytes[..14].iter().all(u8::is_ascii_digit) {
        let year = &trimmed[0..4];
        let month = &trimmed[4..6];
        let day = &trimmed[6..8];
        let hour = &trimmed[8..10];
        let min = &trimmed[10..12];
        let sec = &trimmed[12..14];
        let rest = trimmed[14..].trim();
        let tz = if rest.is_empty() || rest.eq_ignore_ascii_case("z") {
            "Z".to_owned()
        } else if rest.len() == 5
            && (rest.starts_with('+') || rest.starts_with('-'))
            && rest[1..].bytes().all(|b| b.is_ascii_digit())
        {
            format!("{}:{}", &rest[0..3], &rest[3..5])
        } else if rest.len() == 6
            && (rest.starts_with('+') || rest.starts_with('-'))
            && rest.as_bytes()[3] == b':'
            && rest[1..3].bytes().all(|b| b.is_ascii_digit())
            && rest[4..6].bytes().all(|b| b.is_ascii_digit())
        {
            rest.to_owned()
        } else {
            return trimmed.to_owned();
        };
        return format!("{year}-{month}-{day}T{hour}:{min}:{sec}{tz}");
    }
    trimmed.to_owned()
}
