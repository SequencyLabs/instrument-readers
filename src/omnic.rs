//! Pure-Rust reader for Thermo Nicolet OMNIC `.spa` infrared spectra.
//!
//! Reads the 16-byte block directory at offset 304 (key 2 = spectrum header, key 3 = float32
//! points), the title at offset 30, y-axis unit code, and acquisition timestamp at offset 296
//! (seconds since 1899-12-31T00:00:00Z).

use crate::binary::{format_unix_timestamp, read_f32_le, read_i32_le, read_u32_le};
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentReading, Technique, Trace, VendorFormat};

const MAGIC: &[u8] = b"Spectral Data File";
const DIR_OFFSET: usize = 304;
const DIR_ENTRIES: usize = 64;
const DIR_ENTRY_LEN: usize = 16;
const OMNIC_EPOCH_UNIX_SECS: i64 = -2_209_075_200;

/// Lightweight header metadata parsed from an OMNIC `.spa` file without decoding spectrum points.
#[derive(Debug, Clone, PartialEq)]
pub struct OmnicHeader {
    /// Spectrum title stored at offset 30.
    pub title: String,
    /// Canonical y-axis unit (`"reflectance_pct"`, `"transmittance_pct"`, or `"absorbance"`).
    pub y_unit: Option<String>,
    /// Acquisition timestamp as seconds since 1899-12-31T00:00:00Z.
    pub acquired_epoch_secs: u32,
}

/// Parsed OMNIC `.spa` spectrum with wavenumber axis, intensity values, metadata, and timestamp.
#[derive(Debug, Clone, PartialEq)]
pub struct OmnicSpectrum {
    /// Ascending wavenumbers in cm^-1.
    pub wavenumbers: Vec<f64>,
    /// Intensity values matching `wavenumbers`.
    pub intensities: Vec<f64>,
    /// Spectrum title stored at offset 30.
    pub title: String,
    /// Canonical y-axis unit (`"reflectance_pct"`, `"transmittance_pct"`, or `"absorbance"`).
    pub y_unit: Option<String>,
    /// Acquisition timestamp as seconds since 1899-12-31T00:00:00Z.
    pub acquired_epoch_secs: u32,
}

fn find_blocks(data: &[u8]) -> Result<(usize, Option<usize>)> {
    if !data.starts_with(MAGIC) {
        return Err(InstrumentError::InvalidOmnic(
            "not an OMNIC spectrum".to_owned(),
        ));
    }
    let mut header_off: Option<usize> = None;
    let mut points_off: Option<usize> = None;
    for i in 0..DIR_ENTRIES {
        let base = DIR_OFFSET + i * DIR_ENTRY_LEN;
        let Some(entry) = data.get(base..base + DIR_ENTRY_LEN) else {
            break;
        };
        let key = entry[0];
        let off = u32::from_le_bytes([entry[2], entry[3], entry[4], entry[5]]) as usize;
        if key == 2 && header_off.is_none() {
            header_off = Some(off);
        } else if key == 3 && points_off.is_none() {
            points_off = Some(off);
        }
    }
    let hdr = header_off.ok_or_else(|| {
        InstrumentError::InvalidOmnic("OMNIC spectrum without a header".to_owned())
    })?;
    Ok((hdr, points_off))
}

fn decode_omnic_title(data: &[u8]) -> Result<String> {
    let title_slice = data
        .get(30..285)
        .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC title".to_owned()))?;
    let title_end = title_slice
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(title_slice.len());
    let raw = &title_slice[..title_end];
    if let Ok(utf8) = std::str::from_utf8(raw) {
        return Ok(utf8.to_owned());
    }
    let (decoded, _, _) = encoding_rs::WINDOWS_1252.decode(raw);
    Ok(decoded.into_owned())
}

fn parse_header_at(data: &[u8], hdr: usize) -> Result<OmnicHeader> {
    let title = decode_omnic_title(data)?;
    let acquired_epoch_secs = read_u32_le(data, 296)
        .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC timestamp".to_owned()))?;
    let y_unit = data
        .get(hdr.checked_add(12).ok_or_else(|| {
            InstrumentError::InvalidOmnic("OMNIC header offset overflow".to_owned())
        })?)
        .and_then(|&code| match code {
            11 => Some("reflectance_pct".to_owned()),
            16 => Some("transmittance_pct".to_owned()),
            17 => Some("absorbance".to_owned()),
            _ => None,
        });
    Ok(OmnicHeader {
        title,
        y_unit,
        acquired_epoch_secs,
    })
}

/// Parses only the header metadata (title, y-axis unit, timestamp) of an OMNIC `.spa` buffer.
pub fn parse_header(data: &[u8]) -> Result<OmnicHeader> {
    let (hdr, _) = find_blocks(data)?;
    parse_header_at(data, hdr)
}

/// Parses an OMNIC `.spa` buffer into an [`OmnicSpectrum`].
pub fn parse_spa(data: &[u8]) -> Result<OmnicSpectrum> {
    let (hdr, points_off) = find_blocks(data)?;
    let pts = points_off
        .ok_or_else(|| InstrumentError::InvalidOmnic("OMNIC spectrum without points".to_owned()))?;
    let header = parse_header_at(data, hdr)?;

    let npts_i32 = read_i32_le(data, hdr + 4)
        .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC point count".to_owned()))?;
    let npts = usize::try_from(npts_i32)
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| InstrumentError::InvalidOmnic("invalid OMNIC point count".to_owned()))?;

    let wmax = f64::from(
        read_f32_le(data, hdr + 16)
            .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC wmax".to_owned()))?,
    );
    let wmin = f64::from(
        read_f32_le(data, hdr + 20)
            .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC wmin".to_owned()))?,
    );

    let byte_len = npts
        .checked_mul(4)
        .ok_or_else(|| InstrumentError::InvalidOmnic("OMNIC point count overflow".to_owned()))?;
    let end = pts
        .checked_add(byte_len)
        .ok_or_else(|| InstrumentError::InvalidOmnic("OMNIC points offset overflow".to_owned()))?;
    let raw_pts = data
        .get(pts..end)
        .ok_or_else(|| InstrumentError::InvalidOmnic("truncated OMNIC points".to_owned()))?;

    let mut wavenumbers = Vec::with_capacity(npts);
    let mut intensities = Vec::with_capacity(npts);
    let step = if npts > 1 {
        (wmin - wmax) / ((npts - 1) as f64)
    } else {
        0.0
    };

    for rev_idx in (0..npts).rev() {
        let wn = wmax + (rev_idx as f64) * step;
        let off = rev_idx * 4;
        let val = f64::from(f32::from_le_bytes([
            raw_pts[off],
            raw_pts[off + 1],
            raw_pts[off + 2],
            raw_pts[off + 3],
        ]));
        wavenumbers.push(wn);
        intensities.push(val);
    }

    Ok(OmnicSpectrum {
        wavenumbers,
        intensities,
        title: header.title,
        y_unit: header.y_unit,
        acquired_epoch_secs: header.acquired_epoch_secs,
    })
}

/// Parses an OMNIC `.spa` buffer into a unified [`InstrumentReading`].
pub fn read_spa(data: &[u8]) -> Result<InstrumentReading> {
    let spec = parse_spa(data)?;
    let acquired_at = format_unix_timestamp(
        OMNIC_EPOCH_UNIX_SECS + i64::from(spec.acquired_epoch_secs),
        0,
        'T',
        "+00:00",
    );
    let y_unit = spec
        .y_unit
        .clone()
        .unwrap_or_else(|| "intensity".to_owned());
    let points = spec.wavenumbers.into_iter().zip(spec.intensities).collect();
    let sample_name = (!spec.title.trim().is_empty()).then(|| spec.title.trim().to_owned());

    Ok(InstrumentReading {
        technique: Technique::Ir,
        vendor_format: VendorFormat::OmnicSpa,
        acquired_at: Some(acquired_at),
        sample_name,
        traces: vec![Trace {
            label: "IR".to_owned(),
            x_unit: "cm-1".to_owned(),
            y_unit,
            points,
        }],
        ..InstrumentReading::default()
    })
}
