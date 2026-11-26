//! Pure-Rust reader for Waters MassLynx `.raw` LC-MS directories.
//!
//! Decodes `_HEADER.TXT` metadata, `_extern.inf` MS polarities and parameters, `_FUNCTNS.INF`
//! detector types and channel labels, `_FUNC*.IDX` scan indexes, and `_FUNC*.DAT` 2-, 4-, 6-,
//! and 8-byte intensity records into lossless acquisitions, 1D TIC/UV traces, and 2D matrices.

mod decode;
mod metadata;
mod readiness;

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentReading, Technique, Trace, VendorFormat};

use self::decode::{DecodedMatrix, decode_function_matrix, decode_function_trace, parse_idx};
use self::metadata::{collect_bundle, is_uv_function, parse_functns_inf, parse_header_txt};

pub use self::readiness::{WatersBundleFile, is_bundle_ready};

pub(crate) const IDX_RECORD_LEN: usize = 22;
pub(crate) const FUNCTNS_RECORD_LEN: usize = 416;
pub(crate) const FUNCTNS_MZ_OFFSET: usize = 160;
pub(crate) const FUNCTNS_MZ_SLOTS: usize = 32;
pub(crate) const FUNC_TYPE_MASK: u8 = 0x1F;
pub(crate) const FUNC_TYPE_MS_SCAN: u8 = 0;
pub(crate) const FUNC_TYPE_SIR: u8 = 1;
pub(crate) const FUNC_TYPE_DIODE_ARRAY: u8 = 12;
pub(crate) const MAX_MATRIX_CELLS: usize = 50_000_000;

/// Detector/acquisition function type stored in `_FUNCTNS.INF` (`record[0] & 0x1F`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WatersFunctionType {
    MsScan,
    Sir,
    DiodeArray,
    Other(u8),
}

impl WatersFunctionType {
    #[must_use]
    pub const fn from_raw(raw: u8) -> Self {
        match raw & FUNC_TYPE_MASK {
            FUNC_TYPE_MS_SCAN => Self::MsScan,
            FUNC_TYPE_SIR => Self::Sir,
            FUNC_TYPE_DIODE_ARRAY => Self::DiodeArray,
            other => Self::Other(other),
        }
    }
}

/// Ion polarity parsed from `_extern.inf`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonPolarity {
    Positive,
    Negative,
}

impl IonPolarity {
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::Positive => '+',
            Self::Negative => '-',
        }
    }

    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Positive => "+",
            Self::Negative => "-",
        }
    }
}

/// All parameters parsed for one function section in `_extern.inf` (repeated keys preserved in order).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WatersExternFunction {
    pub polarity: Option<IonPolarity>,
    pub params: Vec<(String, String)>,
}

/// Lossless header metadata parsed from `_HEADER.TXT`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct WatersHeader {
    pub sample_name: Option<String>,
    pub acquired_at: Option<String>,
    pub instrument: Option<String>,
    pub inlet_method: Option<String>,
    pub ms_method: Option<String>,
    /// Polynomial calibration coefficients keyed by 1-based function number.
    pub calibrations: BTreeMap<u32, Vec<f64>>,
}

/// Readability status of a declared Waters function's `.DAT`/`.IDX` binary pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WatersFunctionStatus {
    Complete {
        scan_count: usize,
        bytes_per_pair: usize,
    },
    Incomplete,
}

/// Lossless metadata and scan status for one declared Waters function.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersFunctionSummary {
    pub func_num: u32,
    pub func_type: Option<WatersFunctionType>,
    pub polarity: Option<IonPolarity>,
    pub channel_labels: Vec<f64>,
    pub calibration: Vec<f64>,
    pub extern_params: Vec<(String, String)>,
    pub status: WatersFunctionStatus,
}

/// Complete, lossless view of a Waters MassLynx `.raw` acquisition with zero dropped functions.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersAcquisition {
    pub header: WatersHeader,
    pub functions: BTreeMap<u32, WatersFunctionSummary>,
}

/// One UV channel extracted from a Waters MassLynx `.raw` run.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersUvChannel {
    /// Wavelength or channel label from `_FUNCTNS.INF`.
    pub wavelength: f64,
    /// Absorbance values aligned to `uv_rt`.
    pub values: Vec<f64>,
}

/// Full 2D MS + UV polarity slice parsed from a Waters MassLynx `.raw` directory.
#[derive(Debug, Clone, PartialEq)]
pub struct WatersPolarityRun {
    /// Ion polarity (`"+"` or `"-"`).
    pub polarity: String,
    /// Acquisition date and time string from `_HEADER.TXT`.
    pub acquired: Option<String>,
    /// Instrument model from `_HEADER.TXT`.
    pub instrument: Option<String>,
    /// Inlet method filename from `_HEADER.TXT`.
    pub inlet_method: Option<String>,
    /// MS method filename from `_HEADER.TXT`.
    pub ms_method: Option<String>,
    /// Retention times in minutes for the MS scans.
    pub rt: Vec<f64>,
    /// Binned m/z axis values.
    pub mz: Vec<f64>,
    /// Row-major intensity matrix (`rt.len() * mz.len()`).
    pub intensity: Vec<f64>,
    /// Retention times in minutes for the UV channels, if present.
    pub uv_rt: Option<Vec<f64>>,
    /// UV channels aligned to `uv_rt`.
    pub uv: Vec<WatersUvChannel>,
}

/// Returns `true` when `idx_bytes` describes at least one complete, readable scan whose record
/// width is consistent across all indexed scans and fits within `dat_len`.
#[must_use]
pub fn pair_has_complete_scan(idx_bytes: &[u8], dat_len: usize) -> bool {
    parse_idx(idx_bytes, dat_len).is_some()
}

/// Parses a Waters MassLynx `.raw` bundle into a lossless [`WatersAcquisition`] containing every
/// declared function's detector type, polarity, calibration, repeated parameters, and scan status
/// without silently dropping incomplete functions.
pub fn parse_acquisition(files: &[(&str, &[u8])]) -> Result<WatersAcquisition> {
    let bundle = collect_bundle(files)?;
    let func_infos = bundle
        .functns_bytes
        .map(parse_functns_inf)
        .unwrap_or_default();

    let mut declared = BTreeSet::new();
    for idx in 0..func_infos.len() {
        declared.insert((idx as u32) + 1);
    }
    declared.extend(bundle.func_dats.keys().copied());
    declared.extend(bundle.func_idxs.keys().copied());
    declared.extend(bundle.extern_funcs.keys().copied());
    declared.extend(bundle.meta.calibrations.keys().copied());

    let mut functions = BTreeMap::new();
    for func_num in declared {
        let status = match (
            bundle.func_dats.get(&func_num),
            bundle.func_idxs.get(&func_num),
        ) {
            (Some(&dat_bytes), Some(&idx_bytes)) => match parse_idx(idx_bytes, dat_bytes.len()) {
                Some((scans, fmt)) => WatersFunctionStatus::Complete {
                    scan_count: scans.len(),
                    bytes_per_pair: fmt.record_len(),
                },
                None => WatersFunctionStatus::Incomplete,
            },
            _ => WatersFunctionStatus::Incomplete,
        };
        let func_info = func_num
            .checked_sub(1)
            .and_then(|idx| func_infos.get(idx as usize));
        let extern_fn = bundle.extern_funcs.get(&func_num);
        functions.insert(
            func_num,
            WatersFunctionSummary {
                func_num,
                func_type: func_info.map(|info| WatersFunctionType::from_raw(info.func_type)),
                polarity: extern_fn.and_then(|e| e.polarity),
                channel_labels: func_info
                    .map(|info| info.labels.clone())
                    .unwrap_or_default(),
                calibration: bundle
                    .meta
                    .calibrations
                    .get(&func_num)
                    .cloned()
                    .unwrap_or_default(),
                extern_params: extern_fn.map(|e| e.params.clone()).unwrap_or_default(),
                status,
            },
        );
    }
    Ok(WatersAcquisition {
        header: bundle.meta,
        functions,
    })
}

/// Parses a Waters MassLynx `_HEADER.TXT` file into an [`InstrumentReading`] with metadata only.
pub fn read_header(data: &[u8]) -> Result<InstrumentReading> {
    let meta = parse_header_txt(data)?;
    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::WatersRaw,
        acquired_at: meta.acquired_at,
        sample_name: meta.sample_name,
        instrument: meta.instrument,
        ..InstrumentReading::default()
    })
}

/// Parses a Waters MassLynx `.raw` directory bundle into an [`InstrumentReading`].
pub fn read_bundle(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    let bundle = collect_bundle(files)?;
    let has_functns = bundle.functns_bytes.is_some();
    let func_infos = bundle
        .functns_bytes
        .map(parse_functns_inf)
        .unwrap_or_default();
    let mut traces = Vec::new();

    for (&func_num, &dat_bytes) in &bundle.func_dats {
        let Some(&idx_bytes) = bundle.func_idxs.get(&func_num) else {
            continue;
        };
        if dat_bytes.is_empty() || idx_bytes.len() < IDX_RECORD_LEN {
            continue;
        }
        let Some((scans, fmt)) = parse_idx(idx_bytes, dat_bytes.len()) else {
            continue;
        };
        let Some(points) = decode_function_trace(dat_bytes, &scans, fmt) else {
            continue;
        };
        if points.is_empty() {
            continue;
        }

        let polarity = bundle
            .extern_funcs
            .get(&func_num)
            .and_then(|e| e.polarity)
            .map(IonPolarity::as_char);
        let func_info = func_num
            .checked_sub(1)
            .and_then(|idx| func_infos.get(idx as usize));
        let (label, y_unit) = if is_uv_function(has_functns, func_info, polarity) {
            (format!("UV (Func {func_num})"), "AU")
        } else {
            match polarity {
                Some('+') => (format!("ESI+ TIC (Func {func_num})"), "counts"),
                Some('-') => (format!("ESI- TIC (Func {func_num})"), "counts"),
                _ => (format!("MS TIC (Func {func_num})"), "counts"),
            }
        };
        traces.push(Trace {
            label,
            x_unit: "min".to_owned(),
            y_unit: y_unit.to_owned(),
            points,
        });
    }

    if traces.is_empty() && bundle.meta.sample_name.is_none() && bundle.meta.instrument.is_none() {
        return Err(InstrumentError::InvalidWaters(
            "bundle contains no readable _HEADER.TXT or _FUNC*.DAT/IDX files".to_owned(),
        ));
    }

    Ok(InstrumentReading {
        technique: Technique::LcMs,
        vendor_format: VendorFormat::WatersRaw,
        acquired_at: bundle.meta.acquired_at,
        sample_name: bundle.meta.sample_name,
        instrument: bundle.meta.instrument,
        traces,
        ..InstrumentReading::default()
    })
}

/// Parses every MS polarity (`"+"` and `"-"`) and UV channel from a Waters MassLynx `.raw` bundle.
pub fn read_polarities(
    files: &[(&str, &[u8])],
    bin_width: f64,
) -> Result<BTreeMap<String, WatersPolarityRun>> {
    if !bin_width.is_finite() || bin_width <= 0.0 {
        return Err(InstrumentError::InvalidWaters(
            "bin_width must be positive and finite".to_owned(),
        ));
    }

    let bundle = collect_bundle(files)?;
    if let Some(inf) = bundle.functns_bytes {
        let expected = FUNCTNS_RECORD_LEN
            .checked_mul(bundle.func_dats.len())
            .ok_or_else(|| InstrumentError::InvalidWaters("function count overflow".to_owned()))?;
        if inf.len() != expected {
            return Err(InstrumentError::InvalidWaters(format!(
                "_FUNCTNS.INF size {} does not match {} functions",
                inf.len(),
                bundle.func_dats.len()
            )));
        }
    }

    let has_functns = bundle.functns_bytes.is_some();
    let func_infos = bundle
        .functns_bytes
        .map(parse_functns_inf)
        .unwrap_or_default();

    let mut uv_rt: Option<Vec<f64>> = None;
    let mut uv_channels: Vec<WatersUvChannel> = Vec::new();
    let mut ms_by_polarity: BTreeMap<String, DecodedMatrix> = BTreeMap::new();

    for (&func_num, &dat_bytes) in &bundle.func_dats {
        let Some(&idx_bytes) = bundle.func_idxs.get(&func_num) else {
            continue;
        };
        let Some((scans, fmt)) = parse_idx(idx_bytes, dat_bytes.len()) else {
            return Err(InstrumentError::InvalidWaters(format!(
                "unreadable index for _FUNC{func_num:03}.DAT"
            )));
        };

        let polarity = bundle
            .extern_funcs
            .get(&func_num)
            .and_then(|e| e.polarity)
            .map(IonPolarity::as_char);
        let calib = bundle
            .meta
            .calibrations
            .get(&func_num)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let func_info = func_num
            .checked_sub(1)
            .and_then(|idx| func_infos.get(idx as usize));
        let inf_labels = func_info.map(|info| info.labels.as_slice()).unwrap_or(&[]);
        let is_uv = is_uv_function(has_functns, func_info, polarity);

        if !is_uv && polarity.is_none() {
            continue;
        }

        let decoded = decode_function_matrix(
            dat_bytes,
            &scans,
            fmt,
            inf_labels,
            calib,
            if is_uv { 1.0 } else { bin_width },
            is_uv,
        )?;

        if is_uv {
            if uv_rt.is_none() {
                let nrows = decoded.rt.len();
                let ncols = decoded.ylabels.len();
                let mut chans = Vec::with_capacity(ncols);
                for (col, &wl) in decoded.ylabels.iter().enumerate() {
                    let mut vals = Vec::with_capacity(nrows);
                    for row in 0..nrows {
                        vals.push(decoded.matrix[row * ncols + col]);
                    }
                    chans.push(WatersUvChannel {
                        wavelength: wl,
                        values: vals,
                    });
                }
                uv_rt = Some(decoded.rt);
                uv_channels = chans;
            }
        } else if let Some(pol) = polarity {
            ms_by_polarity.entry(pol.to_string()).or_insert(decoded);
        }
    }

    let mut out = BTreeMap::new();
    for (polarity, ms) in ms_by_polarity {
        out.insert(
            polarity.clone(),
            WatersPolarityRun {
                polarity,
                acquired: bundle.meta.acquired_at.clone(),
                instrument: bundle.meta.instrument.clone(),
                inlet_method: bundle.meta.inlet_method.clone(),
                ms_method: bundle.meta.ms_method.clone(),
                rt: ms.rt,
                mz: ms.ylabels,
                intensity: ms.matrix,
                uv_rt: uv_rt.clone(),
                uv: uv_channels.clone(),
            },
        );
    }

    Ok(out)
}
