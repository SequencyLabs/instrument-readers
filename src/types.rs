//! Unified domain model for analytical instrument readings, methods, peaks, and traces.

use std::fmt;

use serde::{Deserialize, Serialize};

/// Analytical chemistry technique of an instrument run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Technique {
    /// High-performance liquid chromatography (or gas chromatography with liquid-style channels).
    Hplc,
    /// Liquid or gas chromatography coupled to mass spectrometry.
    LcMs,
    /// Nuclear magnetic resonance spectroscopy.
    Nmr,
    /// Infrared spectroscopy.
    Ir,
    /// Ultraviolet-visible spectrophotometry.
    Uv,
    /// Technique is not known from the file alone.
    #[default]
    Unknown,
}

impl Technique {
    /// Canonical lowercase identifier for serialization and display.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Hplc => "hplc",
            Self::LcMs => "lcms",
            Self::Nmr => "nmr",
            Self::Ir => "ir",
            Self::Uv => "uv",
            Self::Unknown => "unknown",
        }
    }
}

impl fmt::Display for Technique {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Vendor container or file format identified by magic bytes or file naming.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VendorFormat {
    /// Shimadzu LabSolutions `.lcd` OLE2 compound file.
    ShimadzuLcd,
    /// Bruker TopSpin NMR dataset (`acqus`, `procs`, `1r`).
    BrukerNmr,
    /// Agilent OpenLab CDS `.acaml` XML metadata file.
    OpenLabAcaml,
    /// Agilent OpenLab CDS `.dx` zipped raw signal container.
    OpenLabDx,
    /// Agilent OpenLab CDS `.rx` zipped integration result container.
    OpenLabRx,
    /// Generic OLE2 Compound File Binary Format container without Shimadzu streams.
    Ole2,
    /// Thermo Finnigan `.raw` or `.sld` binary file.
    ThermoRaw,
    /// ANDI Chromatography netCDF `.cdf` file.
    AndiCdf,
    /// Thermo Nicolet OMNIC `.spa` IR spectrum.
    OmnicSpa,
    /// Mestrelab MestReNova `.mnova` document.
    MestReNova,
    /// Thermo Dionex Chromeleon `.raw` signal stream.
    ChromeleonRaw,
    /// Waters MassLynx `.raw` folder (`_HEADER.TXT`).
    WatersRaw,
    /// Agilent MassHunter `.d` folder (`sample_info.xml`).
    AgilentMassHunter,
    /// Plain-text UV export (`.all`, `.spe`, `.sp2`).
    UvText,
    /// Unrecognized format.
    #[default]
    Unknown,
}

impl VendorFormat {
    /// Canonical snake-case identifier of the vendor format.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ShimadzuLcd => "shimadzu_lcd",
            Self::BrukerNmr => "bruker_nmr",
            Self::OpenLabAcaml => "openlab_acaml",
            Self::OpenLabDx => "openlab_dx",
            Self::OpenLabRx => "openlab_rx",
            Self::Ole2 => "ole2",
            Self::ThermoRaw => "thermo_raw",
            Self::AndiCdf => "andi_cdf",
            Self::OmnicSpa => "omnic_spa",
            Self::MestReNova => "mestrenova",
            Self::ChromeleonRaw => "chromeleon_raw",
            Self::WatersRaw => "waters_raw",
            Self::AgilentMassHunter => "agilent_masshunter",
            Self::UvText => "uv_text",
            Self::Unknown => "unknown",
        }
    }

    /// Default analytical technique associated with this vendor format.
    #[must_use]
    pub const fn default_technique(self) -> Technique {
        match self {
            Self::ShimadzuLcd
            | Self::OpenLabAcaml
            | Self::OpenLabDx
            | Self::OpenLabRx
            | Self::AndiCdf
            | Self::ChromeleonRaw => Technique::Hplc,
            Self::ThermoRaw | Self::WatersRaw | Self::AgilentMassHunter => Technique::LcMs,
            Self::BrukerNmr | Self::MestReNova => Technique::Nmr,
            Self::OmnicSpa => Technique::Ir,
            Self::UvText => Technique::Uv,
            Self::Ole2 | Self::Unknown => Technique::Unknown,
        }
    }
}

impl fmt::Display for VendorFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One step in a chromatographic mobile-phase gradient program.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GradientStep {
    /// Elapsed run time in minutes.
    pub time_min: f64,
    /// Mobile phase line B percentage (0 to 100).
    pub b_pct: f64,
    /// Mobile phase line C percentage (0 to 100).
    pub c_pct: f64,
    /// Mobile phase line D percentage (0 to 100).
    pub d_pct: f64,
}

impl GradientStep {
    /// Implied mobile phase line A percentage (`100 - B - C - D`, clamped at zero).
    #[must_use]
    pub fn a_pct(&self) -> f64 {
        (100.0 - self.b_pct - self.c_pct - self.d_pct).max(0.0)
    }
}

/// Unified acquisition and processing method parameters across instrument techniques.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InstrumentMethod {
    /// Method file name or identifier, if recorded.
    pub name: Option<String>,
    /// Mobile-phase solvent names (chromatography) or lock solvent (NMR).
    pub solvents: Vec<String>,
    /// Pump timetable steps for chromatographic gradients.
    pub gradient: Vec<GradientStep>,
    /// Pump flow rate in mL/min.
    pub flow_rate_ml_min: Option<f64>,
    /// Column oven temperature in degrees Celsius.
    pub column_temp_c: Option<f64>,
    /// Injected sample volume in microliters.
    pub injection_volume_ul: Option<f64>,
    /// Primary UV/PDA detection wavelength in nanometers.
    pub wavelength_nm: Option<f64>,
    /// Observed NMR nucleus, such as `"1H"`, `"13C"`, or `"19F"`.
    pub nucleus: Option<String>,
    /// Spectrometer carrier/reference frequency in megahertz.
    pub frequency_mhz: Option<f64>,
    /// NMR pulse program name, such as `"zg30"`.
    pub pulse_program: Option<String>,
    /// Number of accumulated scans (`NS`).
    pub scans: Option<u32>,
    /// Relaxation delay `D1` in seconds.
    pub relaxation_delay_s: Option<f64>,
}

/// Integrated chromatographic peak from a vendor peak table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Peak {
    /// Retention time at the peak apex in minutes.
    pub rt_min: Option<f64>,
    /// Integrated peak area in detector units.
    pub area: Option<f64>,
    /// Relative area percentage across all integrated peaks in the channel.
    pub area_pct: Option<f64>,
    /// Peak height above baseline in detector units.
    pub height: Option<f64>,
    /// Assigned compound name from the processing method, if any.
    pub compound: Option<String>,
}

/// Continuous 1D chromatogram or spectrum trace ready for plotting.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Trace {
    /// Human-readable channel or trace name, such as `"Max absorbance"` or `"1H NMR"`.
    pub label: String,
    /// Unit of the horizontal axis, such as `"min"` or `"ppm"`.
    pub x_unit: String,
    /// Unit of the vertical axis, such as `"mAU"`, `"pA"`, or `"intensity"`.
    pub y_unit: String,
    /// Ordered `(x, y)` coordinate pairs.
    pub points: Vec<(f64, f64)>,
}

impl Trace {
    /// Downsamples `self.points` to at most `threshold` points using Largest-Triangle-Three-Buckets.
    #[must_use]
    pub fn downsample_lttb(&self, threshold: usize) -> Self {
        Self {
            label: self.label.clone(),
            x_unit: self.x_unit.clone(),
            y_unit: self.y_unit.clone(),
            points: lttb(&self.points, threshold),
        }
    }

    /// Downsamples `self.points` into `bins` horizontal buckets keeping each bucket's min and max.
    ///
    /// Preserves sharp chromatographic and NMR peaks across up to `2 * bins` output points.
    #[must_use]
    pub fn downsample_minmax(&self, bins: usize) -> Self {
        Self {
            label: self.label.clone(),
            x_unit: self.x_unit.clone(),
            y_unit: self.y_unit.clone(),
            points: minmax_bins(&self.points, bins),
        }
    }
}

/// Unified reading parsed from an analytical instrument run.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InstrumentReading {
    /// Analytical technique (`Hplc`, `LcMs`, `Nmr`, `Ir`, `Uv`, or `Unknown`).
    pub technique: Technique,
    /// Detected vendor container format.
    pub vendor_format: VendorFormat,
    /// Acquisition timestamp as an ISO-8601 or vendor-normalized UTC/offset string.
    pub acquired_at: Option<String>,
    /// Sample name recorded on the instrument.
    pub sample_name: Option<String>,
    /// Instrument model or name, when available.
    pub instrument: Option<String>,
    /// Acquisition and processing method parameters.
    pub method: InstrumentMethod,
    /// Integrated peaks reported by the vendor's processing table, if present.
    pub peaks: Vec<Peak>,
    /// Extracted chromatogram or spectrum traces.
    pub traces: Vec<Trace>,
}

impl InstrumentReading {
    /// Returns the peak with the largest integrated area, matching `sequency_lab_data.instruments.main_peak`.
    #[must_use]
    pub fn main_peak(&self) -> Option<&Peak> {
        self.peaks
            .iter()
            .filter(|p| p.area.is_some())
            .max_by(|a, b| {
                let lhs = a.area.unwrap_or(f64::NEG_INFINITY);
                let rhs = b.area.unwrap_or(f64::NEG_INFINITY);
                lhs.total_cmp(&rhs)
            })
    }
}

/// Downsamples `(x, y)` points to `threshold` points via Largest-Triangle-Three-Buckets (LTTB).
#[must_use]
pub fn lttb(points: &[(f64, f64)], threshold: usize) -> Vec<(f64, f64)> {
    let n = points.len();
    if threshold == 0 || n == 0 {
        return Vec::new();
    }
    if threshold >= n {
        return points.to_vec();
    }
    if threshold == 1 {
        return vec![points[0]];
    }
    if threshold == 2 {
        return vec![points[0], points[n - 1]];
    }

    let mut sampled = Vec::with_capacity(threshold);
    sampled.push(points[0]);

    let every = (n - 2) as f64 / (threshold - 2) as f64;
    let mut a_idx = 0usize;

    for i in 0..(threshold - 2) {
        let avg_start = (((i + 1) as f64 * every).floor() as usize) + 1;
        let avg_end = ((((i + 2) as f64 * every).floor() as usize) + 1).min(n);
        let avg_slice = &points[avg_start.min(avg_end)..avg_end];

        let (avg_x, avg_y) = if avg_slice.is_empty() {
            points[n - 1]
        } else {
            let len = avg_slice.len() as f64;
            let sum_x: f64 = avg_slice.iter().map(|p| p.0).sum();
            let sum_y: f64 = avg_slice.iter().map(|p| p.1).sum();
            (sum_x / len, sum_y / len)
        };

        let range_start = ((i as f64 * every).floor() as usize) + 1;
        let range_end = ((((i + 1) as f64 * every).floor() as usize) + 1).min(n - 1);
        let (ax, ay) = points[a_idx];

        let mut best_idx = range_start;
        let mut max_area = -1.0_f64;
        for (offset, &(px, py)) in points[range_start..range_end.max(range_start + 1)]
            .iter()
            .enumerate()
        {
            let area = ((ax - avg_x) * (py - ay) - (ax - px) * (avg_y - ay)).abs();
            if area > max_area {
                max_area = area;
                best_idx = range_start + offset;
            }
        }

        sampled.push(points[best_idx]);
        a_idx = best_idx;
    }

    sampled.push(points[n - 1]);
    sampled
}

/// Downsamples `(x, y)` points into `bins` index buckets, emitting each bucket's min and max in index order.
///
/// Mirrors `sequency_lab_data.traces.envelope` so narrow spikes never disappear on Canvas 2D plots.
#[must_use]
pub fn minmax_bins(points: &[(f64, f64)], bins: usize) -> Vec<(f64, f64)> {
    let n = points.len();
    if bins == 0 || n == 0 {
        return Vec::new();
    }
    if n <= bins {
        return points.to_vec();
    }

    let mut out = Vec::with_capacity(bins.saturating_mul(2).min(n));
    for bucket in 0..bins {
        let start = (bucket * n) / bins;
        let end = ((bucket + 1) * n) / bins;
        if end <= start {
            continue;
        }

        let mut min_idx = start;
        let mut max_idx = start;
        for idx in (start + 1)..end {
            if points[idx].1.total_cmp(&points[min_idx].1).is_lt() {
                min_idx = idx;
            }
            if points[idx].1.total_cmp(&points[max_idx].1).is_gt() {
                max_idx = idx;
            }
        }

        if min_idx == max_idx {
            out.push(points[min_idx]);
        } else if min_idx < max_idx {
            out.push(points[min_idx]);
            out.push(points[max_idx]);
        } else {
            out.push(points[max_idx]);
            out.push(points[min_idx]);
        }
    }
    out
}
