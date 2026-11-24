//! Shared bounds-checked binary and text decoding helpers.

use std::borrow::Cow;

use encoding_rs::{UTF_16BE, UTF_16LE};

pub(crate) fn read_u16_le(data: &[u8], offset: usize) -> Option<u16> {
    let bytes: [u8; 2] = data.get(offset..offset.checked_add(2)?)?.try_into().ok()?;
    Some(u16::from_le_bytes(bytes))
}

pub(crate) fn read_i16_be(data: &[u8], offset: usize) -> Option<i16> {
    let bytes: [u8; 2] = data.get(offset..offset.checked_add(2)?)?.try_into().ok()?;
    Some(i16::from_be_bytes(bytes))
}

pub(crate) fn read_u32_le(data: &[u8], offset: usize) -> Option<u32> {
    let bytes: [u8; 4] = data.get(offset..offset.checked_add(4)?)?.try_into().ok()?;
    Some(u32::from_le_bytes(bytes))
}

pub(crate) fn read_i32_le(data: &[u8], offset: usize) -> Option<i32> {
    let bytes: [u8; 4] = data.get(offset..offset.checked_add(4)?)?.try_into().ok()?;
    Some(i32::from_le_bytes(bytes))
}

pub(crate) fn read_i32_be(data: &[u8], offset: usize) -> Option<i32> {
    let bytes: [u8; 4] = data.get(offset..offset.checked_add(4)?)?.try_into().ok()?;
    Some(i32::from_be_bytes(bytes))
}

pub(crate) fn read_f32_be(data: &[u8], offset: usize) -> Option<f32> {
    let bytes: [u8; 4] = data.get(offset..offset.checked_add(4)?)?.try_into().ok()?;
    Some(f32::from_be_bytes(bytes))
}

pub(crate) fn read_u64_le(data: &[u8], offset: usize) -> Option<u64> {
    let bytes: [u8; 8] = data.get(offset..offset.checked_add(8)?)?.try_into().ok()?;
    Some(u64::from_le_bytes(bytes))
}

pub(crate) fn read_f64_le(data: &[u8], offset: usize) -> Option<f64> {
    let bytes: [u8; 8] = data.get(offset..offset.checked_add(8)?)?.try_into().ok()?;
    Some(f64::from_le_bytes(bytes))
}

pub(crate) fn read_f64_be(data: &[u8], offset: usize) -> Option<f64> {
    let bytes: [u8; 8] = data.get(offset..offset.checked_add(8)?)?.try_into().ok()?;
    Some(f64::from_be_bytes(bytes))
}

/// Extracts the final component of a Windows or POSIX file path.
pub(crate) fn file_name(path: &str) -> &str {
    let trimmed = path.trim_end_matches(['\\', '/']);
    if trimmed.is_empty() {
        return "";
    }
    trimmed.rsplit(['\\', '/']).next().unwrap_or("")
}

/// Rounds a floating-point value to `decimals` decimal places.
pub(crate) fn round_digits(val: f64, decimals: i32) -> f64 {
    let factor = 10.0_f64.powi(decimals);
    (val * factor).round() / factor
}

/// Decodes UTF-16LE bytes, stripping an optional BOM and trailing NUL characters.
pub(crate) fn decode_utf16le(data: &[u8]) -> String {
    let body = data.strip_prefix(&[0xFF, 0xFE]).unwrap_or(data);
    let even_len = body.len() & !1;
    let (cow, _) = UTF_16LE.decode_without_bom_handling(&body[..even_len]);
    cow.trim_end_matches('\0').to_owned()
}

/// Decodes XML bytes encoded as UTF-8, UTF-16LE, or UTF-16BE into UTF-8 text.
pub(crate) fn decode_text_auto(data: &[u8]) -> Cow<'_, str> {
    if let Some(rest) = data.strip_prefix(&[0xFF, 0xFE]) {
        let even = rest.len() & !1;
        let (cow, _) = UTF_16LE.decode_without_bom_handling(&rest[..even]);
        return Cow::Owned(cow.into_owned());
    }
    if let Some(rest) = data.strip_prefix(&[0xFE, 0xFF]) {
        let even = rest.len() & !1;
        let (cow, _) = UTF_16BE.decode_without_bom_handling(&rest[..even]);
        return Cow::Owned(cow.into_owned());
    }
    if data.len() >= 4 && data[1] == 0 && data[3] == 0 && data[0] != 0 {
        let even = data.len() & !1;
        let (cow, _) = UTF_16LE.decode_without_bom_handling(&data[..even]);
        return Cow::Owned(cow.into_owned());
    }
    let body = data.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(data);
    String::from_utf8_lossy(body)
}

/// Decodes an even-length ASCII hex slice into bytes.
pub(crate) fn decode_hex(hex: &[u8]) -> Option<Vec<u8>> {
    if !hex.len().is_multiple_of(2) {
        return None;
    }
    let mut out = Vec::with_capacity(hex.len() / 2);
    let (pairs, _) = hex.as_chunks::<2>();
    for &[hi_b, lo_b] in pairs {
        let hi = hex_nibble(hi_b)?;
        let lo = hex_nibble(lo_b)?;
        out.push((hi << 4) | lo);
    }
    Some(out)
}

pub(crate) fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

/// Appends an XML general or character entity reference from `quick-xml` into `buf`.
pub(crate) fn append_xml_ref(buf: &mut String, r: &quick_xml::events::BytesRef<'_>) {
    if let Ok(Some(ch)) = r.resolve_char_ref() {
        buf.push(ch);
        return;
    }
    match r.as_ref() {
        "amp" => buf.push('&'),
        "lt" => buf.push('<'),
        "gt" => buf.push('>'),
        "quot" => buf.push('"'),
        "apos" => buf.push('\''),
        other => {
            buf.push('&');
            buf.push_str(other);
            buf.push(';');
        }
    }
}

/// Converts a `quick-xml` local tag name into an owned `String`.
pub(crate) fn xml_local_name(local: quick_xml::name::LocalName<'_>) -> String {
    local.into_inner().to_owned()
}

/// Extracts and unescapes an XML attribute matching `name` by local name.
pub(crate) fn xml_attr(e: &quick_xml::events::BytesStart<'_>, name: &str) -> Option<String> {
    for attr in e.attributes().flatten() {
        if attr.key.local_name().as_ref() == name {
            return quick_xml::escape::unescape(&attr.value)
                .ok()
                .map(|value| value.into_owned());
        }
    }
    None
}

/// Converts days since `1970-01-01` into `(year, month, day)` in the proleptic Gregorian calendar.
pub(crate) const fn civil_from_days(days_since_unix_epoch: i64) -> (i64, u32, u32) {
    let z = days_since_unix_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if m <= 2 { y + 1 } else { y };
    (year, m as u32, d as u32)
}

/// Formats a Unix timestamp in seconds at UTC+08:00 as `YYYY-MM-DD HH:MM:SS +0800`.
pub(crate) fn format_unix_plus_0800(unix_secs: i64) -> String {
    let local = unix_secs + 8 * 3600;
    let days = local.div_euclid(86_400);
    let secs_of_day = local.rem_euclid(86_400) as u32;
    let (year, month, day) = civil_from_days(days);
    let hour = secs_of_day / 3600;
    let minute = (secs_of_day % 3600) / 60;
    let second = secs_of_day % 60;
    format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}:{second:02} +0800")
}
