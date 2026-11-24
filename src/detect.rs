//! Automatic instrument format detection from filenames and magic bytes, plus top-level dispatch.

use crate::binary::{decode_text_auto, file_name};
use crate::error::{InstrumentError, Result};
use crate::ole::{OLE_MAGIC, OleFile};
use crate::types::{InstrumentReading, VendorFormat};
use crate::zip_archive::{ZIP_EOCD_MAGIC, ZIP_LOCAL_MAGIC, ZipArchive};
use crate::{bruker, openlab, shimadzu};

const MNOVA_MAGIC: &[u8] = b"Mestrelab Research S.L.\xf1\xe2\xd3\xc4";
const FINNIGAN_UTF16: &[u8] = b"F\0i\0n\0n\0i\0g\0a\0n\0";

/// Detects the [`VendorFormat`] of an instrument file using magic bytes and optional `filename`.
#[must_use]
pub fn detect(filename: Option<&str>, data: &[u8]) -> VendorFormat {
    let base = filename.map(file_name).unwrap_or("");
    let base_lower = base.to_ascii_lowercase();

    if let Some(by_magic) = detect_by_magic(&base_lower, data) {
        return by_magic;
    }
    detect_by_filename(filename.unwrap_or(""), &base_lower)
}

fn detect_by_magic(base_lower: &str, data: &[u8]) -> Option<VendorFormat> {
    if data.starts_with(&OLE_MAGIC) {
        return Some(detect_ole2(base_lower, data));
    }

    if data.starts_with(&ZIP_LOCAL_MAGIC) || data.starts_with(&ZIP_EOCD_MAGIC) {
        if base_lower.ends_with(".rx") {
            return Some(VendorFormat::OpenLabRx);
        }
        if base_lower.ends_with(".dx") {
            return Some(VendorFormat::OpenLabDx);
        }
        if let Ok(zip) = ZipArchive::parse(data) {
            let mut has_injection_acaml = false;
            let mut has_dx_member = false;
            for entry in zip.entries() {
                if entry.name == "Base/InjectionACAML" {
                    has_injection_acaml = true;
                } else if entry.name.eq_ignore_ascii_case("injection.acmd")
                    || (entry.name.len() > 3
                        && entry.name[entry.name.len() - 3..].eq_ignore_ascii_case(".ch"))
                {
                    has_dx_member = true;
                }
            }
            if has_injection_acaml {
                return Some(VendorFormat::OpenLabRx);
            }
            if has_dx_member {
                return Some(VendorFormat::OpenLabDx);
            }
        }
    }

    if data.starts_with(MNOVA_MAGIC) {
        return Some(VendorFormat::MestReNova);
    }

    if data.len() >= 18 && &data[..2] == b"\x01\xa1" && &data[2..18] == FINNIGAN_UTF16 {
        return Some(VendorFormat::ThermoRaw);
    }

    if data.starts_with(b"SignHdr\0") {
        return Some(VendorFormat::ChromeleonRaw);
    }

    if data.starts_with(b"CDF\x01") || data.starts_with(b"CDF\x02") {
        return Some(VendorFormat::AndiCdf);
    }

    if data.starts_with(b"Spectral Data File")
        || (data.len() >= 30 && &data[..2] == b"\x00\x04" && base_lower.ends_with(".spa"))
    {
        return Some(VendorFormat::OmnicSpa);
    }

    let head = &data[..data.len().min(4096)];
    if !head.is_empty() {
        let text = decode_text_auto(head);
        let trimmed = text.trim_start();
        if looks_like_jcamp(trimmed) {
            return Some(VendorFormat::BrukerNmr);
        }
        if trimmed.contains("<ACAML") || trimmed.contains("urn:schemas-agilent-com:acaml21") {
            return Some(VendorFormat::OpenLabAcaml);
        }
        if trimmed.contains("<SampleInfo>") && trimmed.contains("<Field>") {
            return Some(VendorFormat::AgilentMassHunter);
        }
        if trimmed.starts_with("$$ Version") || trimmed.contains("$$ Acquired Name:") {
            return Some(VendorFormat::WatersRaw);
        }
    }

    None
}

fn detect_by_filename(full_path: &str, base_lower: &str) -> VendorFormat {
    if base_lower.ends_with(".lcd") {
        return VendorFormat::ShimadzuLcd;
    }
    if base_lower.ends_with(".lcm") {
        return VendorFormat::ShimadzuLcm;
    }
    if base_lower.ends_with(".acaml") {
        return VendorFormat::OpenLabAcaml;
    }
    if base_lower.ends_with(".dx") {
        return VendorFormat::OpenLabDx;
    }
    if base_lower.ends_with(".rx") {
        return VendorFormat::OpenLabRx;
    }
    if matches!(
        base_lower,
        "acqus" | "acqu" | "acqu2s" | "procs" | "proc" | "proc2s" | "1r" | "1i" | "fid" | "ser"
    ) {
        return VendorFormat::BrukerNmr;
    }
    if base_lower.ends_with(".cdf") {
        return VendorFormat::AndiCdf;
    }
    if base_lower.ends_with(".spa") {
        return VendorFormat::OmnicSpa;
    }
    if base_lower.ends_with(".mnova") {
        return VendorFormat::MestReNova;
    }
    if base_lower.ends_with(".all")
        || base_lower.ends_with(".spe")
        || base_lower.ends_with(".sp2")
        || base_lower.ends_with(".sp3")
        || base_lower.ends_with(".sp4")
    {
        return VendorFormat::UvText;
    }

    let path_lower = full_path.to_ascii_lowercase();
    if base_lower == "_header.txt" || path_lower.contains(".raw/") || path_lower.contains(".raw\\")
    {
        return VendorFormat::WatersRaw;
    }
    if base_lower == "sample_info.xml" || path_lower.contains(".d/") || path_lower.contains(".d\\")
    {
        return VendorFormat::AgilentMassHunter;
    }
    if base_lower.ends_with(".raw") || base_lower.ends_with(".sld") {
        return VendorFormat::ThermoRaw;
    }

    VendorFormat::Unknown
}

fn detect_ole2(base_lower: &str, data: &[u8]) -> VendorFormat {
    if base_lower.ends_with(".lcm") {
        return VendorFormat::ShimadzuLcm;
    }
    if base_lower.ends_with(".lcd") {
        return VendorFormat::ShimadzuLcd;
    }
    OleFile::parse(data)
        .ok()
        .and_then(|ole| shimadzu::classify(&ole))
        .unwrap_or(VendorFormat::Ole2)
}

fn looks_like_jcamp(trimmed: &str) -> bool {
    if trimmed.starts_with("##TITLE=")
        || trimmed.starts_with("##JCAMPDX=")
        || trimmed.starts_with("##$")
    {
        return true;
    }
    trimmed.starts_with("$$ ") && trimmed.contains("\n##")
}

/// Automatically detects and parses a single instrument file into an [`InstrumentReading`].
pub fn read(filename: Option<&str>, data: &[u8]) -> Result<InstrumentReading> {
    let format = detect(filename, data);
    match format {
        VendorFormat::ShimadzuLcd => shimadzu::read_lcd(data),
        VendorFormat::ShimadzuLcm => {
            let mut reading = shimadzu::read_lcm(data)?;
            if reading.method.name.is_none() {
                reading.method.name = filename
                    .map(file_name)
                    .filter(|s| !s.is_empty())
                    .map(ToOwned::to_owned);
            }
            Ok(reading)
        }
        VendorFormat::OpenLabAcaml => openlab::read_acaml(data),
        VendorFormat::OpenLabDx => openlab::read_dx(data),
        VendorFormat::OpenLabRx => openlab::read_rx(data),
        VendorFormat::BrukerNmr => {
            let base = filename.map(file_name).unwrap_or("").to_ascii_lowercase();
            if matches!(base.as_str(), "procs" | "proc" | "1r" | "1i") {
                return Err(InstrumentError::InvalidBruker(
                    "Bruker processed spectra require both procs and 1r; use read_files or bruker::read_pdata".to_owned(),
                ));
            }
            bruker::read_acqus(data)
        }
        other => Err(InstrumentError::UnsupportedFormat(other.to_string())),
    }
}

/// Automatically detects and parses a multi-file instrument artifact bundle in memory.
pub fn read_files(files: &[(&str, &[u8])]) -> Result<InstrumentReading> {
    if files.is_empty() {
        return Err(InstrumentError::UnsupportedFormat(
            "empty file set".to_owned(),
        ));
    }

    let find_base = |want: &str| {
        files
            .iter()
            .find(|(p, _)| file_name(p).eq_ignore_ascii_case(want))
            .map(|(_, b)| *b)
    };

    let acqus = find_base("acqus");
    let procs = find_base("procs");
    let data_1r = find_base("1r");
    if let (Some(procs_bytes), Some(r_bytes)) = (procs, data_1r) {
        return bruker::read_pdata(acqus, procs_bytes, r_bytes);
    }
    if let Some(acqus_bytes) = acqus {
        return bruker::read_acqus(acqus_bytes);
    }

    let has_openlab = files.iter().any(|(p, b)| {
        matches!(
            detect(Some(p), b),
            VendorFormat::OpenLabAcaml | VendorFormat::OpenLabDx | VendorFormat::OpenLabRx
        )
    });
    if has_openlab {
        return openlab::read_bundle(files);
    }

    if let Some((path, lcd_bytes)) = files
        .iter()
        .find(|(p, b)| detect(Some(p), b) == VendorFormat::ShimadzuLcd)
    {
        return shimadzu::read_lcd(lcd_bytes).map_err(|err| match err {
            InstrumentError::InvalidOle(msg) => {
                InstrumentError::InvalidShimadzu(format!("{path}: {msg}"))
            }
            other => other,
        });
    }

    read(Some(files[0].0), files[0].1)
}
