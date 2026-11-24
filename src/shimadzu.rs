//! Shimadzu LabSolutions `.lcd` reader: properties, method/gradient, peak tables, and chromatograms.

use std::collections::BTreeMap;

use encoding_rs::GB18030;
use quick_xml::Reader;
use quick_xml::events::Event;

use crate::binary::{
    append_xml_ref, civil_from_days, decode_hex, decode_utf16le, file_name, read_f64_le,
    read_i32_le, read_u16_le, read_u32_le, round_digits, xml_attr,
};
use crate::error::{InstrumentError, Result};
use crate::ole::OleFile;
use crate::types::{
    GradientStep, InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat,
};

const PROCESSED_STORAGE: &str = "LSS Data Processing";
const ORIGINAL_STORAGE: &str = "LSS Data Processing Original";
const MAX_PLOT_STREAM: &str = "PDA 3D Raw Data/Max Plot";
const MAX_PLOT_STATUS: &str = "PDA 3D Raw Data/Max Plot Status";
const CHANNEL_PREFIX: &str = "LSS Raw Data/Chromatogram Ch";
const CHANNEL_STATUS: &str = "LSS Raw Data/Chromatogram Status";
const GUM_CONFIG_STREAM: &str = "GUMM_Information/ShimadzuLC.1/GUM.CONFIG";
const LC_METHOD_STREAM: &str = "GUMM_Information/ShimadzuLC.1/GUC.1.METHOD";
const PDA_METHOD_STREAM: &str = "GUMM_Information/ShimadzuPDA.1/PDA.1.METHOD";

const TRACE_HEADER_LEN: usize = 24;
const PEAK_RECORD_MIN_LEN: usize = 32;
const STATUS_RECORD_LEN: usize = 64;
const FILETIME_UNIX_DIFF_SECS: i64 = 11_644_473_600;

/// Run summary extracted directly from a Shimadzu `.lcd` compound file.
#[derive(Debug, Clone, PartialEq)]
pub struct ShimadzuRun {
    /// UTC injection time as an ISO-8601 string, if recorded.
    pub injected: Option<String>,
    /// Sample name from `File Property`, if non-empty.
    pub sample_name: Option<String>,
    /// Data file basename from `File Property`, if non-empty.
    pub data_file: Option<String>,
    /// Method file basename from `File Property`, if non-empty.
    pub method_file: Option<String>,
    /// Instrument model label, such as `"Shimadzu CBM-20A"`.
    pub instrument: String,
    /// Peaks from Ch1 (or the first available `PT-*` table).
    pub peaks: Vec<Peak>,
}

/// Parses a Shimadzu `.lcd` compound file into a [`ShimadzuRun`] summary.
pub fn read_run(data: &[u8]) -> Result<ShimadzuRun> {
    let ole = OleFile::parse(data)?;
    let prop_bytes = ole.require_stream("File Property")?;
    let props = parse_properties(&prop_bytes);
    let non_empty = |s: &str| (!s.is_empty()).then(|| s.to_owned());
    Ok(ShimadzuRun {
        injected: parse_filetime(
            props.get("dwLowDateTime").map(String::as_str),
            props.get("dwHighDateTime").map(String::as_str),
        ),
        sample_name: props.get("smpl_name").and_then(|s| non_empty(s)),
        data_file: props.get("datafile").and_then(|s| non_empty(file_name(s))),
        method_file: props
            .get("methodfile")
            .and_then(|s| non_empty(file_name(s))),
        instrument: parse_instrument(&ole),
        peaks: read_best_peaks(&ole),
    })
}

/// Parses a Shimadzu `.lcd` compound file into a unified [`InstrumentReading`].
pub fn read_lcd(data: &[u8]) -> Result<InstrumentReading> {
    let ole = OleFile::parse(data)?;
    read_lcd_ole(&ole)
}

/// Parses a Shimadzu `.lcd` from an already-opened [`OleFile`].
pub fn read_lcd_ole(ole: &OleFile<'_>) -> Result<InstrumentReading> {
    let props = match ole.stream("File Property") {
        Some(bytes) => parse_properties(&bytes),
        None => BTreeMap::new(),
    };

    let acquired_at = parse_filetime(
        props.get("dwLowDateTime").map(String::as_str),
        props.get("dwHighDateTime").map(String::as_str),
    );
    let sample_name = props
        .get("smpl_name")
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && !is_masked(s))
        .map(ToOwned::to_owned);
    let method_name = props
        .get("methodfile")
        .map(|s| file_name(s))
        .filter(|s| !s.is_empty() && !is_masked(s))
        .map(ToOwned::to_owned);
    let injection_volume_ul = props.get("inj_vol").and_then(|v| parse_inj_vol(v));

    let instrument = Some(parse_instrument(ole));
    let mut method = parse_method(ole);
    method.name = method_name;
    if method.injection_volume_ul.is_none() {
        method.injection_volume_ul = injection_volume_ul;
    }

    let peaks = read_best_peaks(ole);
    let traces = chromatograms_ole(ole)?;

    Ok(InstrumentReading {
        technique: Technique::Hplc,
        vendor_format: VendorFormat::ShimadzuLcd,
        acquired_at,
        sample_name,
        instrument,
        method,
        peaks,
        traces,
    })
}

/// Extracts key-value XML leaf fields from a Shimadzu `File Property` stream.
///
/// Decodes `@StoX@` hex strings using GB18030 and preserves the first occurrence of each tag.
#[must_use]
pub fn parse_properties(stream: &[u8]) -> BTreeMap<String, String> {
    let mut found = BTreeMap::new();
    let mut i = 0usize;
    while i < stream.len() {
        if stream[i] != b'<' {
            i += 1;
            continue;
        }
        let tag_start = i + 1;
        let mut tag_end = tag_start;
        while tag_end < stream.len() && is_word_byte(stream[tag_end]) {
            tag_end += 1;
        }
        if tag_end == tag_start || tag_end >= stream.len() || stream[tag_end] != b'>' {
            i += 1;
            continue;
        }
        let tag = &stream[tag_start..tag_end];
        let val_start = tag_end + 1;
        let mut val_end = val_start;
        while val_end < stream.len() && stream[val_end] != b'<' {
            val_end += 1;
        }
        let close_len = 2 + tag.len() + 1;
        if val_end + close_len <= stream.len()
            && stream[val_end] == b'<'
            && stream[val_end + 1] == b'/'
            && &stream[val_end + 2..val_end + 2 + tag.len()] == tag
            && stream[val_end + 2 + tag.len()] == b'>'
        {
            if let Ok(key) = std::str::from_utf8(tag) {
                let raw_val = &stream[val_start..val_end];
                found
                    .entry(key.to_owned())
                    .or_insert_with(|| decode_property_value(raw_val));
            }
            i = val_end + close_len;
        } else {
            i = tag_end + 1;
        }
    }
    found
}

fn decode_property_value(raw: &[u8]) -> String {
    if let Some(hex_part) = raw.strip_prefix(b"@StoX@")
        && let Some(bytes) = decode_hex(hex_part)
    {
        let (cow, _) = GB18030.decode_without_bom_handling(&bytes);
        return cow.into_owned();
    }
    raw.iter().map(|&b| b as char).collect()
}

/// Converts a Windows `FILETIME` split across two signed 32-bit decimal strings into ISO-8601 UTC.
#[must_use]
pub fn parse_filetime(low: Option<&str>, high: Option<&str>) -> Option<String> {
    let lo_i32: i32 = low?.trim().parse().ok()?;
    let hi_i32: i32 = high?.trim().parse().ok()?;
    if hi_i32 <= 0 {
        return None;
    }
    let ticks = ((hi_i32 as u64) << 32) | (lo_i32 as u32 as u64);
    let micros = ticks / 10;
    let total_secs = i64::try_from(micros / 1_000_000).ok()?;
    let subsec_us = micros % 1_000_000;

    let unix_secs = total_secs - FILETIME_UNIX_DIFF_SECS;
    let days = unix_secs.div_euclid(86_400);
    let secs_of_day = unix_secs.rem_euclid(86_400) as u32;
    let (year, month, day) = civil_from_days(days);
    let hour = secs_of_day / 3600;
    let minute = (secs_of_day % 3600) / 60;
    let second = secs_of_day % 60;

    if subsec_us == 0 {
        Some(format!(
            "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}+00:00"
        ))
    } else {
        Some(format!(
            "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}.{subsec_us:06}+00:00"
        ))
    }
}

/// Decodes a Shimadzu `inj_vol` property (`@FtoX@` big-endian IEEE-754 float hex or decimal).
#[must_use]
pub fn parse_inj_vol(raw: &str) -> Option<f64> {
    let trimmed = raw.trim();
    if let Some(hex) = trimmed.strip_prefix("@FtoX@") {
        let bytes: [u8; 4] = decode_hex(hex.as_bytes())?.try_into().ok()?;
        let val = f64::from(f32::from_be_bytes(bytes));
        return (val.is_finite() && val >= 0.0).then(|| round_digits(val, 4));
    }
    let val: f64 = trimmed.parse().ok()?;
    (val.is_finite() && val >= 0.0).then_some(val)
}

/// Parses one detector channel's binary `PT-*` peak table (`VER1` 20-byte or legacy 8-byte header).
#[must_use]
pub fn parse_peak_table(stream: &[u8]) -> Vec<Peak> {
    if stream.len() < 8 {
        return Vec::new();
    }
    let (count_raw, header) = if stream.starts_with(b"VER1") {
        if stream.len() < 20 {
            return Vec::new();
        }
        (read_i32_le(stream, 4).unwrap_or(0), 20usize)
    } else {
        (read_i32_le(stream, 0).unwrap_or(0), 8usize)
    };

    if count_raw <= 0 {
        return Vec::new();
    }
    let count = count_raw as usize;
    let record_size = (stream.len() - header) / count;
    if record_size < PEAK_RECORD_MIN_LEN {
        return Vec::new();
    }

    let mut rows = Vec::with_capacity(count);
    let mut total_area = 0.0_f64;
    for i in 0..count {
        let base = header + i * record_size;
        let Some(rt_ms) = read_i32_le(stream, base + 4) else {
            return Vec::new();
        };
        let Some(area) = read_f64_le(stream, base + 8) else {
            return Vec::new();
        };
        let Some(height) = read_f64_le(stream, base + 24) else {
            return Vec::new();
        };
        total_area += area;
        rows.push((rt_ms, area, height));
    }

    let denom = if total_area != 0.0 { total_area } else { 1.0 };
    rows.into_iter()
        .map(|(rt_ms, area, height)| Peak {
            rt_min: Some(round_digits(f64::from(rt_ms) / 60_000.0, 3)),
            area: Some(area),
            area_pct: Some(round_digits(area / denom * 100.0, 2)),
            height: Some(height),
            compound: None,
        })
        .collect()
}

/// Decodes one delta-encoded chromatogram sub-block into cumulative `i64` intensities.
pub fn decode_deltas(sub_block: &[u8]) -> Result<Vec<i64>> {
    let mut out = Vec::with_capacity(sub_block.len());
    let mut acc = 0i64;
    let mut i = 0usize;

    while i < sub_block.len() {
        let lead = sub_block[i];
        if lead == 0x82 {
            i += 1;
            continue;
        }
        let digit = lead >> 4;
        if digit >= 8 {
            return Err(InstrumentError::InvalidShimadzu(format!(
                "delta sign digit {digit} exceeds 4-byte encoding"
            )));
        }
        let size = 1.max(1 + (digit as usize) / 2);
        let mut raw = i64::from(lead & 0x0F);
        for offset in 1..size {
            let b = sub_block.get(i + offset).copied().unwrap_or(0);
            raw = (raw << 8) | i64::from(b);
        }
        let delta = if digit % 2 == 1 {
            raw - (1i64 << (8 * size - 4))
        } else {
            raw
        };
        acc += delta;
        out.push(acc);
        i += size;
    }

    Ok(out)
}

/// Parses a raw chromatogram stream into `(minutes, raw_values)`.
pub fn parse_trace(stream: &[u8]) -> Result<(Vec<f64>, Vec<i64>)> {
    if stream.len() < TRACE_HEADER_LEN {
        return Err(InstrumentError::InvalidShimadzu(
            "chromatogram stream shorter than 24-byte header".to_owned(),
        ));
    }
    let interval_ms = read_u32_le(stream, 4).unwrap_or(0);
    let count = read_u32_le(stream, 8).unwrap_or(0) as usize;
    let declared_len = read_i32_le(stream, 12).unwrap_or(0);
    if declared_len < 0 {
        return Err(InstrumentError::InvalidShimadzu(
            "negative chromatogram byte length".to_owned(),
        ));
    }
    let end = (declared_len as usize).min(stream.len());

    let mut values = Vec::with_capacity(count.min(65_536));
    let mut pos = TRACE_HEADER_LEN;

    while values.len() < count && pos + 2 <= end {
        let block_len = read_u16_le(stream, pos).unwrap_or(0) as usize;
        let block_end = (pos + 2 + block_len).min(stream.len());
        if pos + 2 > block_end {
            break;
        }
        let sub = decode_deltas(&stream[pos + 2..block_end])?;
        values.extend_from_slice(&sub);
        pos = pos.saturating_add(block_len + 4);
    }

    values.truncate(count);
    let step_min = f64::from(interval_ms) / 60_000.0;
    let minutes = (0..values.len()).map(|idx| idx as f64 * step_min).collect();
    Ok((minutes, values))
}

/// Extracts chromatogram traces from a Shimadzu `.lcd` byte slice.
///
/// Prefers `PDA 3D Raw Data/Max Plot` when non-empty; otherwise reads each `Chromatogram Ch<n>`.
pub fn chromatograms(data: &[u8]) -> Result<Vec<Trace>> {
    let ole = OleFile::parse(data)?;
    chromatograms_ole(&ole)
}

/// Extracts chromatogram traces from an already-parsed [`OleFile`].
pub fn chromatograms_ole(ole: &OleFile<'_>) -> Result<Vec<Trace>> {
    if let Some(pda) =
        read_channel_trace(ole, MAX_PLOT_STREAM, MAX_PLOT_STATUS, 0, "Max absorbance")?
    {
        return Ok(vec![pda]);
    }

    let mut channels: Vec<(usize, String)> = ole
        .stream_names()
        .filter_map(|name| {
            let suffix = name.strip_prefix(CHANNEL_PREFIX)?;
            if suffix.is_empty() || !suffix.bytes().all(|b| b.is_ascii_digit()) {
                return None;
            }
            let num: usize = suffix.parse().ok()?;
            (num > 0).then(|| (num, name.to_owned()))
        })
        .collect();
    channels.sort_by_key(|(n, _)| *n);

    let mut traces = Vec::new();
    for (n, stream_name) in channels {
        let label = format!("Ch{n}");
        if let Some(trace) = read_channel_trace(ole, &stream_name, CHANNEL_STATUS, n - 1, &label)? {
            traces.push(trace);
        }
    }
    Ok(traces)
}

fn read_channel_trace(
    ole: &OleFile<'_>,
    stream_name: &str,
    status_name: &str,
    status_index: usize,
    base_label: &str,
) -> Result<Option<Trace>> {
    if ole.stream_size(stream_name).unwrap_or(0) <= TRACE_HEADER_LEN {
        return Ok(None);
    }
    let Some(raw_stream) = ole.stream(stream_name) else {
        return Ok(None);
    };
    let (minutes, raw_values) = parse_trace(&raw_stream)?;
    if raw_values.is_empty() {
        return Ok(None);
    }

    let (scale, unit) = ole
        .stream(status_name)
        .and_then(|status| parse_status_record(&status, status_index))
        .unwrap_or((1.0, String::new()));

    let label = if unit.is_empty() {
        base_label.to_owned()
    } else {
        format!("{base_label} ({unit})")
    };
    let y_unit = if unit.is_empty() {
        "intensity".to_owned()
    } else {
        unit
    };
    let points = minutes
        .into_iter()
        .zip(raw_values.into_iter().map(|v| (v as f64) * scale))
        .collect();

    Ok(Some(Trace {
        label,
        x_unit: "min".to_owned(),
        y_unit,
        points,
    }))
}

fn parse_status_record(status: &[u8], index: usize) -> Option<(f64, String)> {
    let start = index.checked_mul(STATUS_RECORD_LEN)?;
    let record = status.get(start..start + STATUS_RECORD_LEN)?;
    let cf = read_f64_le(record, 16)?;
    let vf = read_f64_le(record, 32)?;
    let cf_norm = if cf > 0.0 { cf } else { 1.0 };
    let vf_norm = if vf > 0.0 { vf } else { 1.0 };
    let unit_bytes = record[40..].split(|&b| b == 0).next().unwrap_or(&[]);
    let unit: String = unit_bytes.iter().map(|&b| b as char).collect();
    Some((cf_norm / vf_norm, unit.trim().to_owned()))
}

fn read_best_peaks(ole: &OleFile<'_>) -> Vec<Peak> {
    let prefix = format!("{PROCESSED_STORAGE}/");
    let mut tables: Vec<String> = ole
        .stream_names()
        .filter_map(|path| {
            let name = path.strip_prefix(&prefix)?;
            (!name.contains('/') && name.starts_with("PT-")).then(|| name.to_owned())
        })
        .collect();
    tables.sort();

    let Some(chosen) = tables
        .iter()
        .find(|t| t.ends_with(".1"))
        .or_else(|| tables.first())
    else {
        return Vec::new();
    };

    for storage in [PROCESSED_STORAGE, ORIGINAL_STORAGE] {
        let full = format!("{storage}/{chosen}");
        if let Some(bytes) = ole.stream(&full) {
            let peaks = parse_peak_table(&bytes);
            if !peaks.is_empty() {
                return peaks;
            }
        }
    }
    Vec::new()
}

fn parse_instrument(ole: &OleFile<'_>) -> String {
    let Some(config_bytes) = ole.stream(GUM_CONFIG_STREAM) else {
        return "Shimadzu".to_owned();
    };
    let gumm = parse_gumm_xml(&decode_utf16le(&config_bytes));
    if let Some(model) = gumm.model_name.as_deref().map(str::trim)
        && !model.is_empty()
        && !is_masked(model)
    {
        return format!("Shimadzu {model}");
    }
    "Shimadzu".to_owned()
}

fn parse_method(ole: &OleFile<'_>) -> InstrumentMethod {
    let mut method = InstrumentMethod::default();

    if let Some(lc_bytes) = ole.stream(LC_METHOD_STREAM) {
        let lc = parse_gumm_xml(&decode_utf16le(&lc_bytes));
        for key in [
            "PumpA_MobilePhaseName#1A",
            "PumpA_MobilePhaseName#2A",
            "PumpA_MobilePhaseName#3A",
            "PumpA_MobilePhaseName#4A",
        ] {
            if let Some((val, _)) = lc.params.get(key) {
                let trimmed = val.trim();
                if !trimmed.is_empty() && !is_masked(trimmed) {
                    method.solvents.push(trimmed.to_owned());
                }
            }
        }

        method.flow_rate_ml_min = get_upd_float(&lc.params, "PumpA_Flow")
            .or_else(|| get_upd_float(&lc.params, "LPGE$tflow"))
            .or_else(|| get_upd_float(&lc.params, "Flow"));

        if get_upd_float(&lc.params, "OVEN$use").is_none_or(|u| u > 0.0) {
            method.column_temp_c = get_upd_float(&lc.params, "OVEN$ovent");
        }

        method.wavelength_nm = find_wavelength(&lc.params);

        let init_b = get_upd_float(&lc.params, "PumpA_BConc")
            .or_else(|| get_upd_float(&lc.params, "LPGE$bconc"));
        let init_c = get_upd_float(&lc.params, "PumpA_CConc")
            .or_else(|| get_upd_float(&lc.params, "LPGE$cconc"));
        let init_d = get_upd_float(&lc.params, "PumpA_DConc")
            .or_else(|| get_upd_float(&lc.params, "LPGE$dconc"));

        method.gradient = build_gradient(init_b, init_c, init_d, &lc.tp_events);
    }

    if let Some(pda_bytes) = ole.stream(PDA_METHOD_STREAM) {
        let pda = parse_gumm_xml(&decode_utf16le(&pda_bytes));
        if let Some(wl) = find_wavelength(&pda.params) {
            method.wavelength_nm = Some(wl);
        }
    }

    method
}

#[derive(Default)]
struct GummData {
    model_name: Option<String>,
    params: BTreeMap<String, (String, Option<String>)>,
    tp_events: Vec<(f64, i32, f64)>,
}

#[derive(Default)]
struct UpdEntry {
    id: String,
    val: String,
    kind: Option<String>,
}

#[derive(Default)]
struct TimeProgEvent {
    time_min: Option<f64>,
    code: Option<i32>,
    data: Option<f64>,
}

fn parse_gumm_xml(xml: &str) -> GummData {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().trim_text(false);

    let mut out = GummData::default();
    let mut in_model = false;
    let mut model_buf = String::new();
    let mut in_timeprog = false;
    let mut in_tpld = false;
    let mut cur_upd: Option<UpdEntry> = None;
    let mut cur_field: Option<&'static str> = None;
    let mut cur_event = TimeProgEvent::default();

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => match e.local_name().as_ref() {
                "ModelName" if out.model_name.is_none() => {
                    in_model = true;
                    model_buf.clear();
                }
                "TP" => {
                    in_timeprog = xml_attr(&e, "PgN").as_deref() == Some("SCL$TIMEPROG");
                }
                "TPLD" if in_timeprog => {
                    in_tpld = true;
                    cur_event = TimeProgEvent::default();
                }
                "UPD" => {
                    cur_upd = Some(UpdEntry {
                        id: xml_attr(&e, "ID").unwrap_or_default(),
                        val: String::new(),
                        kind: None,
                    });
                }
                "Val" if cur_upd.is_some() => cur_field = Some("Val"),
                "Type" if cur_upd.is_some() => cur_field = Some("Type"),
                _ => {}
            },
            Ok(Event::GeneralRef(r)) => {
                if in_model {
                    append_xml_ref(&mut model_buf, &r);
                } else if let (Some(upd), Some("Val")) = (cur_upd.as_mut(), cur_field) {
                    append_xml_ref(&mut upd.val, &r);
                }
            }
            Ok(Event::Text(t)) => {
                if in_model {
                    model_buf.push_str(&t);
                } else if let (Some(upd), Some(field)) = (cur_upd.as_mut(), cur_field) {
                    match field {
                        "Val" => upd.val.push_str(&t),
                        "Type" => {
                            let s = upd.kind.get_or_insert_with(String::new);
                            s.push_str(t.trim());
                        }
                        _ => {}
                    }
                }
            }
            Ok(Event::End(e)) => match e.local_name().as_ref() {
                "ModelName" if in_model => {
                    in_model = false;
                    out.model_name = Some(std::mem::take(&mut model_buf));
                }
                "Val" | "Type" => cur_field = None,
                "UPD" => {
                    if let Some(upd) = cur_upd.take() {
                        if in_tpld {
                            record_tpld_upd(&upd, &mut cur_event);
                        } else {
                            out.params.entry(upd.id).or_insert((upd.val, upd.kind));
                        }
                    }
                }
                "TPLD" if in_tpld => {
                    in_tpld = false;
                    if let (Some(t), Some(code)) = (cur_event.time_min, cur_event.code) {
                        out.tp_events.push((t, code, cur_event.data.unwrap_or(0.0)));
                    }
                }
                "TP" => in_timeprog = false,
                _ => {}
            },
            Ok(Event::Eof) | Err(_) => break,
            _ => {}
        }
    }

    out
}

fn find_wavelength(params: &BTreeMap<String, (String, Option<String>)>) -> Option<f64> {
    [
        "2DCH_Wavelength#1",
        "PDA$wavelen1",
        "SPD$wavelen1",
        "DetectorA_Wavelength1",
    ]
    .into_iter()
    .find_map(|k| get_upd_float(params, k))
    .filter(|&w| w > 0.0)
    .map(normalize_wavelength_nm)
}

fn record_tpld_upd(upd: &UpdEntry, ev: &mut TimeProgEvent) {
    match upd.id.as_str() {
        "TMPGtime" => ev.time_min = decode_upd_float(&upd.val, upd.kind.as_deref()),
        "TMPGcode" => ev.code = upd.val.trim().parse::<i32>().ok(),
        "TMPGdata" => ev.data = decode_upd_float(&upd.val, upd.kind.as_deref()),
        _ => {}
    }
}

fn build_gradient(
    init_b: Option<f64>,
    init_c: Option<f64>,
    init_d: Option<f64>,
    tp_events: &[(f64, i32, f64)],
) -> Vec<GradientStep> {
    let has_conc_events = tp_events.iter().any(|(_, c, _)| matches!(c, 5..=7));
    if !has_conc_events && init_b.is_none() && init_c.is_none() && init_d.is_none() {
        return Vec::new();
    }

    let mut cur_b = round_digits(init_b.unwrap_or(0.0), 3);
    let mut cur_c = round_digits(init_c.unwrap_or(0.0), 3);
    let mut cur_d = round_digits(init_d.unwrap_or(0.0), 3);

    let mut steps: Vec<GradientStep> = Vec::new();
    let starts_at_zero = tp_events
        .first()
        .is_some_and(|&(t, code, _)| t == 0.0 && matches!(code, 5..=7));
    if !starts_at_zero {
        steps.push(GradientStep {
            time_min: 0.0,
            b_pct: cur_b,
            c_pct: cur_c,
            d_pct: cur_d,
        });
    }

    for &(raw_time, code, raw_val) in tp_events {
        let time_min = round_digits(raw_time, 3);
        let val = round_digits(raw_val, 3);
        match code {
            5 => cur_b = val,
            6 => cur_c = val,
            7 => cur_d = val,
            97 if has_conc_events => {
                if steps.last().is_some_and(|last| time_min > last.time_min) {
                    steps.push(GradientStep {
                        time_min,
                        b_pct: cur_b,
                        c_pct: cur_c,
                        d_pct: cur_d,
                    });
                }
                continue;
            }
            _ => continue,
        }

        if let Some(last) = steps.last_mut()
            && (last.time_min - time_min).abs() < 1e-6
        {
            last.b_pct = cur_b;
            last.c_pct = cur_c;
            last.d_pct = cur_d;
        } else {
            steps.push(GradientStep {
                time_min,
                b_pct: cur_b,
                c_pct: cur_c,
                d_pct: cur_d,
            });
        }
    }

    steps
}

fn normalize_wavelength_nm(raw: f64) -> f64 {
    if raw >= 1000.0 {
        round_digits(raw / 100.0, 2)
    } else {
        round_digits(raw, 2)
    }
}

fn get_upd_float(params: &BTreeMap<String, (String, Option<String>)>, key: &str) -> Option<f64> {
    let (val, kind) = params.get(key)?;
    decode_upd_float(val, kind.as_deref())
}

fn decode_upd_float(val: &str, kind: Option<&str>) -> Option<f64> {
    let trimmed = val.trim();
    if trimmed.is_empty() {
        return None;
    }
    match kind {
        Some("4") => decode_le_hex_f32(trimmed).or_else(|| trimmed.parse::<f64>().ok()),
        Some("5") => decode_le_hex_f64(trimmed).or_else(|| trimmed.parse::<f64>().ok()),
        Some("0" | "1" | "2" | "3") => trimmed.parse::<f64>().ok(),
        _ => {
            if trimmed.len() == 8 && looks_like_hex_float(trimmed) {
                decode_le_hex_f32(trimmed)
            } else if trimmed.len() == 16 && looks_like_hex_float(trimmed) {
                decode_le_hex_f64(trimmed)
            } else {
                trimmed.parse::<f64>().ok()
            }
        }
    }
    .filter(|v| v.is_finite())
}

fn looks_like_hex_float(s: &str) -> bool {
    s.bytes().all(|b| b.is_ascii_hexdigit())
        && (s.bytes().any(|b| b.is_ascii_alphabetic()) || s.starts_with("0000"))
}

fn decode_le_hex_f32(s: &str) -> Option<f64> {
    let bytes: [u8; 4] = decode_hex(s.as_bytes())?.try_into().ok()?;
    Some(f64::from(f32::from_le_bytes(bytes)))
}

fn decode_le_hex_f64(s: &str) -> Option<f64> {
    let bytes: [u8; 8] = decode_hex(s.as_bytes())?.try_into().ok()?;
    Some(f64::from_le_bytes(bytes))
}

fn is_word_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b == b'_'
}

fn is_masked(s: &str) -> bool {
    let trimmed = s.trim();
    !trimmed.is_empty() && trimmed.bytes().all(|b| b == b'X')
}
