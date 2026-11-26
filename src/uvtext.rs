//! Pure-Rust reader for plain-text UV-Vis tabular exports (`.all`, `.spe`, `.sp2`, `.sp3`, `.sp4`).

use crate::binary::decode_text_auto;
use crate::error::{InstrumentError, Result};
use crate::types::{InstrumentReading, Technique, Trace, VendorFormat};

/// Parsed tabular UV-Vis export with ascending wavelengths and one absorbance column per scan.
#[derive(Debug, Clone, PartialEq)]
pub struct UvTextTable {
    /// Ascending wavelengths in nanometers.
    pub wavelengths: Vec<f64>,
    /// Absorbance columns in column-major order (`num_scans` vectors of length `wavelengths.len()`).
    pub scans: Vec<Vec<f64>>,
}

/// Parses a plain-text UV-Vis export into a [`UvTextTable`].
pub fn parse_uvtext(data: &[u8]) -> Result<UvTextTable> {
    let text = decode_text_auto(data);
    let mut rows: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut expected_cols: Option<usize> = None;

    for line in text.lines().skip(1) {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let cols: Vec<f64> = trimmed
            .split_ascii_whitespace()
            .map(|tok| tok.parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|err| InstrumentError::InvalidUvText(format!("invalid numeric row: {err}")))?;
        if cols.iter().any(|v| !v.is_finite()) {
            return Err(InstrumentError::InvalidUvText(
                "non-finite value in UV text export".to_owned(),
            ));
        }
        if cols.len() < 3 {
            return Err(InstrumentError::InvalidUvText(
                "UV text row must have at least 3 columns".to_owned(),
            ));
        }
        let scan_cols = cols.len() - 2;
        if let Some(prev) = expected_cols {
            if prev != scan_cols {
                return Err(InstrumentError::InvalidUvText(
                    "inconsistent column count in UV text export".to_owned(),
                ));
            }
        } else {
            expected_cols = Some(scan_cols);
        }
        rows.push((cols[1], cols[2..].to_vec()));
    }

    if rows.is_empty() {
        return Err(InstrumentError::InvalidUvText(
            "UV text export contains no data rows".to_owned(),
        ));
    }

    rows.sort_by(|a, b| a.0.total_cmp(&b.0));

    let num_scans = expected_cols.unwrap_or(0);
    let mut wavelengths = Vec::with_capacity(rows.len());
    let mut scans = vec![Vec::with_capacity(rows.len()); num_scans];
    for (wl, vals) in rows {
        wavelengths.push(wl);
        for (col, val) in vals.into_iter().enumerate() {
            scans[col].push(val);
        }
    }

    Ok(UvTextTable { wavelengths, scans })
}

/// Parses a plain-text UV-Vis export into a unified [`InstrumentReading`].
pub fn read_uvtext(data: &[u8]) -> Result<InstrumentReading> {
    let table = parse_uvtext(data)?;
    let traces = table
        .scans
        .into_iter()
        .enumerate()
        .map(|(i, col)| Trace {
            label: format!("UV scan {}", i + 1),
            x_unit: "nm".to_owned(),
            y_unit: "AU".to_owned(),
            points: table.wavelengths.iter().copied().zip(col).collect(),
        })
        .collect();

    Ok(InstrumentReading {
        technique: Technique::Uv,
        vendor_format: VendorFormat::UvText,
        traces,
        ..InstrumentReading::default()
    })
}
