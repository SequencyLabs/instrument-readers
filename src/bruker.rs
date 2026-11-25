//! Bruker TopSpin NMR reader: `acqus`/`procs` JCAMP-DX parameters and `1r` binary spectra.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::binary::{decode_text_auto, format_unix_plus_0800};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentMethod, InstrumentReading, Technique, Trace, VendorFormat};

/// Endianness of a Bruker binary `1r` processed spectrum (`BYTORDP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ByteOrder {
    /// Little-endian (`BYTORDP = 0`).
    #[default]
    LittleEndian,
    /// Big-endian (`BYTORDP = 1`).
    BigEndian,
}

/// Numeric sample representation in a Bruker `1r` binary file (`DTYPP`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum PdataType {
    /// Signed 32-bit integers (`DTYPP = 0`, four bytes per point).
    #[default]
    Int32,
    /// IEEE-754 64-bit floats (`DTYPP = 2`, eight bytes per point).
    Float64,
}

/// Parsed JCAMP-DX parameter table from an `acqus` or `procs` text file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct JcampParams {
    /// Acquisition/processing stamp (`YYYY-MM-DD HH:MM:SS +ZZZZ`), if present in `$$` comments.
    pub timestamp: Option<String>,
    /// Dataset path preceding `/acqus` in `$$` comments, if present.
    pub data_path: Option<String>,
    /// Scalar `##KEY=` / `##$KEY=` entries with surrounding `<...>` stripped.
    pub scalars: BTreeMap<String, String>,
    /// Indexed `##$KEY= (0..N)` array entries as whitespace-split tokens.
    pub arrays: BTreeMap<String, Vec<String>>,
}

/// Acquisition parameters parsed from a Bruker `acqus` JCAMP-DX file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AcqusParams {
    /// Acquisition timestamp (`YYYY-MM-DD HH:MM:SS +ZZZZ`).
    pub acquired_at: Option<String>,
    /// Observed nucleus (`NUC1`), such as `"1H"` or `"13C"`.
    pub nucleus: String,
    /// Pulse program name (`PULPROG`), such as `"zg30"`.
    pub pulse_program: String,
    /// Spectrometer name (`INSTRUM`), such as `"Avance NEO 600"`.
    pub instrument: String,
    /// Deuterated lock solvent (`SOLVENT`), such as `"CDCl3"` or `"DMSO"`.
    pub solvent: String,
    /// Number of scans (`NS`).
    pub scans: u32,
    /// Relaxation delay `D1` (`D[1]`) in seconds.
    pub relaxation_delay_s: f64,
    /// Transmitter frequency `SFO1` in MHz.
    pub sfo1_mhz: Option<f64>,
    /// Base spectrometer frequency `BF1` in MHz.
    pub bf1_mhz: Option<f64>,
    /// Spectral width `SW_h` in Hz.
    pub sw_h_hz: Option<f64>,
    /// Dataset path from the `$$ .../acqus` header comment.
    pub data_path: Option<String>,
}

/// Processing parameters parsed from a Bruker `procs` JCAMP-DX file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProcsParams {
    /// Binary byte order (`BYTORDP`: 0 = little-endian, 1 = big-endian).
    pub byte_order: ByteOrder,
    /// Binary element type (`DTYPP`: 0 = `i32`, 2 = `f64`).
    pub data_type: PdataType,
    /// Intensity scaling exponent (`NC_proc`): true intensity equals `raw * 2^NC_proc`.
    pub nc_proc: i32,
    /// Chemical shift at point 0 (high-ppm edge, `OFFSET`).
    pub offset_ppm: Option<f64>,
    /// Spectrometer reference frequency (`SF`) in MHz.
    pub sf_mhz: Option<f64>,
    /// Processed spectral width (`SW_p`) in Hz.
    pub sw_p_hz: Option<f64>,
    /// Number of real points (`SI`).
    pub si: Option<usize>,
    /// Processed axis nucleus (`AXNUC`), when recorded.
    pub axnuc: Option<String>,
}

/// Parses any Bruker JCAMP-DX parameter file (`acqus`, `procs`, `acqu2s`, `proc2s`).
#[must_use]
pub fn parse_jcamp(text: &str) -> JcampParams {
    let mut out = JcampParams::default();
    let mut active_array: Option<String> = None;

    for raw_line in text.lines() {
        let line = raw_line.trim_end_matches('\r');
        if let Some(comment) = line.strip_prefix("$$ ") {
            active_array = None;
            if out.timestamp.is_none()
                && let Some(stamp) = parse_comment_stamp(comment)
            {
                out.timestamp = Some(stamp);
            }
            if out.data_path.is_none()
                && let Some(path) = parse_acqus_path(comment)
            {
                out.data_path = Some(path);
            }
            continue;
        }
        if line.starts_with("$$") {
            active_array = None;
            continue;
        }
        if let Some(rest) = line.strip_prefix("##") {
            active_array = None;
            let Some((raw_key, raw_val)) = rest.split_once('=') else {
                continue;
            };
            let key = raw_key.strip_prefix('$').unwrap_or(raw_key).trim();
            if key.is_empty() || key == "END" {
                continue;
            }
            let val = raw_val.trim();
            if let Some(after_paren) = strip_array_spec(val) {
                let mut items: Vec<String> = after_paren
                    .split_whitespace()
                    .map(ToOwned::to_owned)
                    .collect();
                let entry = out.arrays.entry(key.to_owned()).or_default();
                entry.append(&mut items);
                active_array = Some(key.to_owned());
            } else {
                let cleaned = strip_angle_brackets(val);
                out.scalars.insert(key.to_owned(), cleaned.to_owned());
            }
            continue;
        }

        if let Some(arr_key) = active_array.as_ref()
            && let Some(vec) = out.arrays.get_mut(arr_key)
        {
            for tok in line.split_whitespace() {
                vec.push(tok.to_owned());
            }
        }
    }

    out
}

/// Parses Bruker `acqus` text into [`AcqusParams`], returning empty defaults when fields are absent.
#[must_use]
pub fn parse_acqus_fields(text: &str) -> AcqusParams {
    acqus_from_jcamp(parse_jcamp(text), false)
}

/// Parses Bruker `acqus` text into [`AcqusParams`], matching `sequency_lab_data.instruments.bruker.parse_acqus`.
pub fn parse_acqus(text: &str) -> Result<AcqusParams> {
    let jcamp = parse_jcamp(text);
    if jcamp.scalars.is_empty() && jcamp.arrays.is_empty() && jcamp.timestamp.is_none() {
        return Err(InstrumentError::InvalidBruker(
            "acqus contains no JCAMP-DX parameters".to_owned(),
        ));
    }
    Ok(acqus_from_jcamp(jcamp, true))
}

fn acqus_from_jcamp(jcamp: JcampParams, fallback_timestamp: bool) -> AcqusParams {
    let get_str = |k: &str| jcamp.scalars.get(k).cloned().unwrap_or_default();
    let get_f64 = |k: &str| {
        jcamp
            .scalars
            .get(k)
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite())
    };

    let acquired_at = jcamp
        .scalars
        .get("DATE")
        .filter(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()))
        .and_then(|v| v.parse::<i64>().ok())
        .filter(|&ts| ts > 0)
        .map(format_unix_plus_0800)
        .or_else(|| fallback_timestamp.then_some(jcamp.timestamp).flatten());

    let scans = jcamp
        .scalars
        .get("NS")
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(0);
    let relaxation_delay_s = jcamp
        .arrays
        .get("D")
        .and_then(|d| d.get(1))
        .and_then(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite())
        .unwrap_or(0.0);

    AcqusParams {
        acquired_at,
        nucleus: get_str("NUC1"),
        pulse_program: get_str("PULPROG"),
        instrument: get_str("INSTRUM"),
        solvent: get_str("SOLVENT"),
        scans,
        relaxation_delay_s,
        sfo1_mhz: get_f64("SFO1"),
        bf1_mhz: get_f64("BF1"),
        sw_h_hz: get_f64("SW_h"),
        data_path: jcamp.data_path,
    }
}

/// Parses Bruker `acqus` bytes into [`AcqusParams`].
pub fn parse_acqus_bytes(data: &[u8]) -> Result<AcqusParams> {
    let text = decode_text_auto(data);
    parse_acqus(&text)
}

/// Parses Bruker `procs` text into [`ProcsParams`].
pub fn parse_procs(text: &str) -> Result<ProcsParams> {
    let jcamp = parse_jcamp(text);
    if jcamp.scalars.is_empty() {
        return Err(InstrumentError::InvalidBruker(
            "procs contains no JCAMP-DX parameters".to_owned(),
        ));
    }

    let byte_order = match jcamp.scalars.get("BYTORDP").map(String::as_str) {
        None | Some("0") => ByteOrder::LittleEndian,
        Some("1") => ByteOrder::BigEndian,
        Some(other) => {
            return Err(InstrumentError::InvalidBruker(format!(
                "unsupported BYTORDP value: {other}"
            )));
        }
    };

    let data_type = match jcamp.scalars.get("DTYPP").map(String::as_str) {
        None | Some("0") => PdataType::Int32,
        Some("2") => PdataType::Float64,
        Some(other) => {
            return Err(InstrumentError::InvalidBruker(format!(
                "unsupported DTYPP value: {other}"
            )));
        }
    };

    let nc_proc = match jcamp.scalars.get("NC_proc") {
        None => 0,
        Some(s) => s
            .parse::<i32>()
            .map_err(|_| InstrumentError::InvalidBruker(format!("invalid NC_proc integer: {s}")))?,
    };

    let get_f64 = |k: &str| {
        jcamp
            .scalars
            .get(k)
            .and_then(|v| v.parse::<f64>().ok())
            .filter(|v| v.is_finite())
    };

    let si = jcamp
        .scalars
        .get("SI")
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0);
    let axnuc = jcamp
        .scalars
        .get("AXNUC")
        .filter(|s| !s.is_empty() && *s != "off")
        .cloned();

    Ok(ProcsParams {
        byte_order,
        data_type,
        nc_proc,
        offset_ppm: get_f64("OFFSET"),
        sf_mhz: get_f64("SF"),
        sw_p_hz: get_f64("SW_p"),
        si,
        axnuc,
    })
}

/// Parses Bruker `procs` bytes into [`ProcsParams`].
pub fn parse_procs_bytes(data: &[u8]) -> Result<ProcsParams> {
    let text = decode_text_auto(data);
    parse_procs(&text)
}

/// Decodes a binary `1r` processed real spectrum using `procs` (`BYTORDP` endianness and `NC_proc` scaling).
///
/// Chemical shifts are returned in descending `ppm` order (high `ppm` first), matching TopSpin storage
/// and `nmrglue`'s `ppm_scale()`.
pub fn read_1r(data_1r: &[u8], procs: &ProcsParams, acqus: Option<&AcqusParams>) -> Result<Trace> {
    if data_1r.is_empty() {
        return Err(InstrumentError::InvalidBruker(
            "1r binary spectrum is empty".to_owned(),
        ));
    }

    let raw_intensities = decode_1r_samples(data_1r, procs.byte_order, procs.data_type)?;
    let count = match procs.si {
        Some(si) if si > raw_intensities.len() => {
            return Err(InstrumentError::InvalidBruker(format!(
                "1r has {} points, fewer than procs SI={si}",
                raw_intensities.len()
            )));
        }
        Some(si) => si,
        None => raw_intensities.len(),
    };
    if count == 0 {
        return Err(InstrumentError::InvalidBruker(
            "1r spectrum has zero points".to_owned(),
        ));
    }

    let obs_mhz = procs
        .sf_mhz
        .or_else(|| acqus.and_then(|a| a.sfo1_mhz))
        .or_else(|| acqus.and_then(|a| a.bf1_mhz))
        .filter(|&f| f > 0.0)
        .ok_or_else(|| {
            InstrumentError::InvalidBruker("missing positive SF or SFO1 frequency".to_owned())
        })?;

    let sw_hz = acqus
        .and_then(|a| a.sw_h_hz)
        .filter(|&w| w > 0.0)
        .or_else(|| procs.sw_p_hz.filter(|&w| w > 0.0))
        .ok_or_else(|| {
            InstrumentError::InvalidBruker(
                "missing positive SW_h or SW_p spectral width".to_owned(),
            )
        })?;

    let first_ppm = match (acqus.and_then(|a| a.sfo1_mhz), procs.sf_mhz) {
        (Some(sfo1), Some(sf)) if sfo1 > 0.0 && sf > 0.0 => {
            ((sfo1 - sf) * 1.0e6) / sf + sw_hz / (2.0 * sf)
        }
        _ => procs
            .offset_ppm
            .ok_or_else(|| InstrumentError::InvalidBruker("missing OFFSET in procs".to_owned()))?,
    };

    let scale = 2.0_f64.powi(procs.nc_proc);
    let step_ppm = sw_hz / ((count as f64) * obs_mhz);

    let mut points = Vec::with_capacity(count);
    for (idx, &raw) in raw_intensities[..count].iter().enumerate() {
        let ppm = first_ppm - (idx as f64) * step_ppm;
        let intensity = raw * scale;
        points.push((ppm, intensity));
    }

    let nucleus = acqus
        .map(|a| a.nucleus.as_str())
        .filter(|s| !s.is_empty())
        .or(procs.axnuc.as_deref());
    let label = match nucleus {
        Some(nuc) => format!("{nuc} NMR"),
        None => "NMR spectrum".to_owned(),
    };

    Ok(Trace {
        label,
        x_unit: "ppm".to_owned(),
        y_unit: "intensity".to_owned(),
        points,
    })
}

/// Builds a metadata-only [`InstrumentReading`] from a Bruker `acqus` byte slice.
pub fn read_acqus(acqus_bytes: &[u8]) -> Result<InstrumentReading> {
    let acqus = parse_acqus_bytes(acqus_bytes)?;
    Ok(reading_from_parts(Some(&acqus), None, None, None))
}

/// Builds a full [`InstrumentReading`] from optional `acqus`, `procs`, and `1r` byte slices.
pub fn read_pdata(
    acqus_bytes: Option<&[u8]>,
    procs_bytes: &[u8],
    data_1r: &[u8],
) -> Result<InstrumentReading> {
    let acqus = match acqus_bytes {
        Some(b) => Some(parse_acqus_bytes(b)?),
        None => None,
    };
    let procs = parse_procs_bytes(procs_bytes)?;
    let trace = read_1r(data_1r, &procs, acqus.as_ref())?;
    Ok(reading_from_parts(
        acqus.as_ref(),
        Some(&procs),
        Some(trace),
        None,
    ))
}

/// Reads a Bruker experiment directory (`<exp>/acqus` and optional `<exp>/pdata/1/{procs,1r,title}`).
pub fn read_dir(exp_dir: &Path) -> Result<InstrumentReading> {
    let (acqus_path, pdata_dir) = if exp_dir.join("procs").is_file() {
        (
            exp_dir
                .parent()
                .and_then(Path::parent)
                .map(|p| p.join("acqus")),
            exp_dir.to_path_buf(),
        )
    } else {
        (Some(exp_dir.join("acqus")), exp_dir.join("pdata").join("1"))
    };

    let acqus_bytes = acqus_path
        .as_deref()
        .map(read_optional_file)
        .transpose()?
        .flatten();
    let title_bytes = read_optional_file(&pdata_dir.join("title"))?;
    let procs_bytes = read_optional_file(&pdata_dir.join("procs"))?;
    let data_1r = read_optional_file(&pdata_dir.join("1r"))?;

    read_experiment_bytes(
        acqus_bytes.as_deref(),
        procs_bytes.as_deref(),
        data_1r.as_deref(),
        title_bytes.as_deref(),
    )?
    .ok_or_else(|| {
        InstrumentError::InvalidBruker(format!(
            "no acqus or pdata/1/procs found in {}",
            exp_dir.display()
        ))
    })
}

fn read_optional_file(path: &Path) -> Result<Option<Vec<u8>>> {
    match std::fs::read(path) {
        Ok(bytes) => Ok(Some(bytes)),
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(err) => Err(err.into()),
    }
}

fn read_experiment_bytes(
    acqus_bytes: Option<&[u8]>,
    procs_bytes: Option<&[u8]>,
    data_1r: Option<&[u8]>,
    title_bytes: Option<&[u8]>,
) -> Result<Option<InstrumentReading>> {
    let acqus = acqus_bytes.map(parse_acqus_bytes).transpose()?;
    let title = title_bytes.and_then(|b| first_nonempty_line(&decode_text_auto(b)));

    if let (Some(procs_raw), Some(r_raw)) = (procs_bytes, data_1r) {
        let procs = parse_procs_bytes(procs_raw)?;
        let trace = read_1r(r_raw, &procs, acqus.as_ref())?;
        return Ok(Some(reading_from_parts(
            acqus.as_ref(),
            Some(&procs),
            Some(trace),
            title,
        )));
    }

    Ok(acqus
        .as_ref()
        .map(|acq| reading_from_parts(Some(acq), None, None, title)))
}

/// Integrates signal intensity between `lo` and `hi` ppm after subtracting median flank baseline.
///
/// Mirrors `sequency_lab_data.instruments.bruker.integral`.
#[must_use]
pub fn integral(points: &[(f64, f64)], lo: f64, hi: f64, flank: f64) -> f64 {
    let mut flanks: Vec<f64> = points
        .iter()
        .filter_map(|&(ppm, y)| {
            (((lo - flank)..lo).contains(&ppm) || (ppm > hi && ppm <= hi + flank)).then_some(y)
        })
        .collect();

    let baseline = median_in_place(&mut flanks).unwrap_or(0.0);
    points
        .iter()
        .filter(|&&(ppm, _)| ppm >= lo && ppm <= hi)
        .map(|&(_, y)| y - baseline)
        .sum()
}

/// Computes the intensity-weighted chemical shift center of the tallest signal in `[lo, hi]` ppm.
///
/// Mirrors `sequency_lab_data.instruments.bruker.multiplet_center`.
#[must_use]
pub fn multiplet_center(points: &[(f64, f64)], lo: f64, hi: f64, half_width: f64) -> Option<f64> {
    let mut top_ppm = None;
    let mut max_y = f64::NEG_INFINITY;
    for &(ppm, y) in points {
        if ppm >= lo && ppm <= hi && y > max_y {
            max_y = y;
            top_ppm = Some(ppm);
        }
    }
    let center = top_ppm?;
    let floor = 0.1 * max_y;

    let mut weighted_sum = 0.0_f64;
    let mut y_sum = 0.0_f64;
    for &(ppm, y) in points {
        if ppm >= center - half_width && ppm <= center + half_width && y > floor {
            weighted_sum += ppm * y;
            y_sum += y;
        }
    }
    (y_sum != 0.0).then_some(weighted_sum / y_sum)
}

fn reading_from_parts(
    acqus: Option<&AcqusParams>,
    procs: Option<&ProcsParams>,
    trace: Option<Trace>,
    title: Option<String>,
) -> InstrumentReading {
    let nucleus = acqus
        .map(|a| a.nucleus.clone())
        .filter(|s| !s.is_empty())
        .or_else(|| procs.and_then(|p| p.axnuc.clone()));
    let frequency_mhz = procs
        .and_then(|p| p.sf_mhz)
        .or_else(|| acqus.and_then(|a| a.sfo1_mhz))
        .or_else(|| acqus.and_then(|a| a.bf1_mhz));
    let solvents = acqus
        .map(|a| a.solvent.trim())
        .filter(|s| !s.is_empty())
        .map(|s| vec![s.to_owned()])
        .unwrap_or_default();
    let pulse_program = acqus
        .map(|a| a.pulse_program.trim())
        .filter(|s| !s.is_empty())
        .map(ToOwned::to_owned);
    let instrument = acqus
        .map(|a| a.instrument.trim())
        .filter(|s| !s.is_empty())
        .map(ToOwned::to_owned);
    let sample_name = title.or_else(|| acqus.and_then(|a| dataset_name(a.data_path.as_deref())));

    let method = InstrumentMethod {
        name: pulse_program.clone(),
        solvents,
        frequency_mhz,
        nucleus,
        pulse_program,
        scans: acqus.map(|a| a.scans).filter(|&n| n > 0),
        relaxation_delay_s: acqus.map(|a| a.relaxation_delay_s).filter(|&d| d > 0.0),
        ..InstrumentMethod::default()
    };

    InstrumentReading {
        technique: Technique::Nmr,
        vendor_format: VendorFormat::BrukerNmr,
        acquired_at: acqus.and_then(|a| a.acquired_at.clone()),
        sample_name,
        instrument,
        method,
        peaks: Vec::new(),
        traces: trace.into_iter().collect(),
    }
}

fn decode_1r_samples(data: &[u8], byte_order: ByteOrder, data_type: PdataType) -> Result<Vec<f64>> {
    match data_type {
        PdataType::Int32 => {
            let (chunks, rem) = data.as_chunks::<4>();
            if !rem.is_empty() {
                return Err(InstrumentError::InvalidBruker(format!(
                    "1r byte length {} is not a multiple of 4 for int32 data",
                    data.len()
                )));
            }
            let out = chunks
                .iter()
                .map(|&arr| {
                    let val = match byte_order {
                        ByteOrder::LittleEndian => i32::from_le_bytes(arr),
                        ByteOrder::BigEndian => i32::from_be_bytes(arr),
                    };
                    f64::from(val)
                })
                .collect();
            Ok(out)
        }
        PdataType::Float64 => {
            let (chunks, rem) = data.as_chunks::<8>();
            if !rem.is_empty() {
                return Err(InstrumentError::InvalidBruker(format!(
                    "1r byte length {} is not a multiple of 8 for float64 data",
                    data.len()
                )));
            }
            let out = chunks
                .iter()
                .map(|&arr| match byte_order {
                    ByteOrder::LittleEndian => f64::from_le_bytes(arr),
                    ByteOrder::BigEndian => f64::from_be_bytes(arr),
                })
                .collect();
            Ok(out)
        }
    }
}

fn parse_comment_stamp(comment: &str) -> Option<String> {
    let trimmed = comment.trim_start();
    if trimmed.len() < 19 {
        return None;
    }
    let dt = &trimmed[..19];
    let bytes = dt.as_bytes();
    let valid_dt = bytes[4] == b'-'
        && bytes[7] == b'-'
        && bytes[10] == b' '
        && bytes[13] == b':'
        && bytes[16] == b':'
        && bytes
            .iter()
            .enumerate()
            .all(|(i, b)| matches!(i, 4 | 7 | 10 | 13 | 16) || b.is_ascii_digit());
    if !valid_dt {
        return None;
    }

    let rest = &trimmed[19..];
    let after_frac = if let Some(frac) = rest.strip_prefix('.') {
        let digits = frac.bytes().take_while(u8::is_ascii_digit).count();
        if digits == 0 {
            return None;
        }
        &frac[digits..]
    } else {
        rest
    };

    let after_space = after_frac.strip_prefix(' ')?;
    if after_space.len() < 5 {
        return None;
    }
    let tz = &after_space[..5];
    let tz_bytes = tz.as_bytes();
    let valid_tz =
        matches!(tz_bytes[0], b'+' | b'-') && tz_bytes[1..].iter().all(u8::is_ascii_digit);
    valid_tz.then(|| format!("{dt} {tz}"))
}

fn parse_acqus_path(comment: &str) -> Option<String> {
    let trimmed = comment.trim();
    let prefix = trimmed.strip_suffix("/acqus")?;
    if prefix.is_empty() || prefix.contains(char::is_whitespace) {
        return None;
    }
    let bytes = prefix.as_bytes();
    let is_abs = prefix.starts_with('/')
        || (bytes.len() >= 3
            && bytes[0].is_ascii_alphabetic()
            && bytes[1] == b':'
            && bytes[2] == b'/');
    is_abs.then(|| prefix.to_owned())
}

fn dataset_name(data_path: Option<&str>) -> Option<String> {
    let path = data_path?;
    let mut parts = path.rsplit(['/', '\\']).filter(|s| !s.is_empty());
    let last = parts.next()?;
    if last.bytes().all(|b| b.is_ascii_digit()) {
        parts.next().map(ToOwned::to_owned)
    } else {
        Some(last.to_owned())
    }
}

fn first_nonempty_line(text: &str) -> Option<String> {
    text.lines()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .map(ToOwned::to_owned)
}

fn strip_array_spec(val: &str) -> Option<&str> {
    let rest = val.strip_prefix("(0..")?;
    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    rest[digits..].strip_prefix(')')
}

fn strip_angle_brackets(val: &str) -> &str {
    val.strip_prefix('<')
        .and_then(|s| s.strip_suffix('>'))
        .unwrap_or(val)
}

fn median_in_place(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_by(f64::total_cmp);
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        Some(values[mid])
    } else {
        Some((values[mid - 1] + values[mid]) / 2.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_ACQUS: &[u8] = b"##TITLE= Parameter file\n\
##$NUC1= <1H>\n\
##$PULPROG= <zg30>\n\
##$INSTRUM= <Avance NEO 600>\n\
##$SOLVENT= <CDCl3>\n\
##$NS= 16\n\
##$SFO1= 600.13\n\
##$SW_h= 12000.0\n\
##END=\n";

    const SAMPLE_PROCS: &[u8] = b"##TITLE= Parameter file\n\
##$BYTORDP= 0\n\
##$DTYPP= 0\n\
##$NC_proc= 0\n\
##$OFFSET= 12.0\n\
##$SF= 600.13\n\
##$SW_p= 12000.0\n\
##$SI= 2\n\
##$AXNUC= <1H>\n\
##END=\n";

    #[test]
    fn read_experiment_bytes_acqus_only() {
        let reading = read_experiment_bytes(Some(SAMPLE_ACQUS), None, None, Some(b"\n  Exp 1  \n"))
            .unwrap()
            .unwrap();
        assert_eq!(reading.sample_name.as_deref(), Some("Exp 1"));
        assert_eq!(reading.instrument.as_deref(), Some("Avance NEO 600"));
        assert_eq!(reading.method.nucleus.as_deref(), Some("1H"));
        assert_eq!(reading.method.scans, Some(16));
        assert!(reading.traces.is_empty());
    }

    #[test]
    fn read_experiment_bytes_procs_and_1r() {
        let mut data_1r = Vec::new();
        data_1r.extend_from_slice(&100_i32.to_le_bytes());
        data_1r.extend_from_slice(&250_i32.to_le_bytes());

        let reading = read_experiment_bytes(None, Some(SAMPLE_PROCS), Some(&data_1r), None)
            .unwrap()
            .unwrap();
        assert_eq!(reading.method.nucleus.as_deref(), Some("1H"));
        assert_eq!(reading.traces.len(), 1);
        assert_eq!(reading.traces[0].points.len(), 2);
        assert_eq!(reading.traces[0].points[0].1, 100.0);
        assert_eq!(reading.traces[0].points[1].1, 250.0);
    }

    #[test]
    fn read_experiment_bytes_neither() {
        let reading = read_experiment_bytes(None, None, None, Some(b"orphan title")).unwrap();
        assert!(reading.is_none());
    }
}
