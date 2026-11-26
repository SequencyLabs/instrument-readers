//! Pure-Rust reader for Agilent MassHunter `.d` LC-MS runs (`MSScan.bin`, `MSScan.xsd`,
//! `MSTS.xml`, `MSMassCal.bin`, `DefaultMassCal.xml`, and `MSProfile.bin` RLE segments).

use std::collections::BTreeMap;

use crate::binary::{
    decode_text_auto, read_f32_le, read_f64_le, read_i16_le, read_i32_le, read_i64_le, read_u32_le,
};
use crate::error::{InstrumentError, Result};
use crate::xml::{elements, parse_xml, trimmed_text};

const MAX_SPECTRUM_BLOCKS: usize = 8;
const MAX_RLE_POINTS: u32 = 0x00FF_FFFF;
const PROFILE_FORMAT_ID: i64 = 1;
const RUN_LENGTH_MARKER: u32 = 0x90;

/// One MS1 profile scan parsed from `MSScan.bin` with its 10-double mass calibration row.
#[derive(Debug, Clone, PartialEq)]
pub struct MassHunterScan {
    /// Retention time in minutes.
    pub minutes: f64,
    /// Ion polarity (`"+"` or `"-"`).
    pub polarity: String,
    /// Total ion current recorded in `MSScan.bin`.
    pub tic: f64,
    /// Byte offset of the profile segment in `MSProfile.bin`.
    pub offset: u64,
    /// Byte length of the profile segment in `MSProfile.bin`.
    pub length: u64,
    /// Number of flight-time grid points in the profile segment.
    pub points: u32,
    /// Ten mass calibration doubles (`[coeff, base, left, right, c0..c5]`).
    pub calibration: [f64; 10],
    /// Polynomial `ValueUseFlags` bitmask for this scan's calibration ID, if any.
    pub flags: Option<u64>,
}

#[derive(Clone, Copy)]
enum XsdType {
    Byte,
    Short,
    Int,
    Long,
    Float,
    Double,
}

impl XsdType {
    fn from_name(name: &str) -> Option<Self> {
        match name.rsplit(':').next().unwrap_or(name) {
            "byte" => Some(Self::Byte),
            "short" => Some(Self::Short),
            "int" => Some(Self::Int),
            "long" => Some(Self::Long),
            "float" => Some(Self::Float),
            "double" => Some(Self::Double),
            _ => None,
        }
    }

    const fn byte_size(self) -> usize {
        match self {
            Self::Byte => 1,
            Self::Short => 2,
            Self::Int | Self::Float => 4,
            Self::Long | Self::Double => 8,
        }
    }

    fn read_i64(self, data: &[u8], off: usize) -> Result<i64> {
        match self {
            Self::Byte => Ok(i64::from(
                *data
                    .get(off)
                    .ok_or_else(|| invalid("truncated MSScan.bin"))? as i8,
            )),
            Self::Short => Ok(i64::from(
                read_i16_le(data, off).ok_or_else(|| invalid("truncated MSScan.bin"))?,
            )),
            Self::Int => Ok(i64::from(
                read_i32_le(data, off).ok_or_else(|| invalid("truncated MSScan.bin"))?,
            )),
            Self::Long => read_i64_le(data, off).ok_or_else(|| invalid("truncated MSScan.bin")),
            Self::Float | Self::Double => Err(invalid("expected integer field in MSScan.bin")),
        }
    }

    fn read_f64(self, data: &[u8], off: usize) -> Result<f64> {
        match self {
            Self::Float => Ok(f64::from(
                read_f32_le(data, off).ok_or_else(|| invalid("truncated MSScan.bin"))?,
            )),
            Self::Double => read_f64_le(data, off).ok_or_else(|| invalid("truncated MSScan.bin")),
            _ => Ok(self.read_i64(data, off)? as f64),
        }
    }
}

#[derive(Clone, Copy)]
struct FieldSlot {
    offset: usize,
    ty: XsdType,
}

struct ScanLayout {
    full_rec: usize,
    block_size: usize,
    scalar_size: usize,
    scan_time: FieldSlot,
    ion_polarity: FieldSlot,
    ms_level: FieldSlot,
    tic: FieldSlot,
    calibration_id: Option<FieldSlot>,
    spectrum_format_id: FieldSlot,
    spectrum_offset: FieldSlot,
    byte_count: FieldSlot,
    point_count: FieldSlot,
}

struct TofCalibrator {
    coeff: f64,
    base: f64,
    lo: f64,
    hi: f64,
    poly: [f64; 64],
    max_order: Option<usize>,
}

impl TofCalibrator {
    fn new(calib: &[f64; 10], use_flags: Option<u64>) -> Result<Self> {
        if calib.iter().any(|v| !v.is_finite()) {
            return Err(invalid("non-finite mass calibration value"));
        }
        let mut poly = [0.0f64; 64];
        let mut max_order = None;
        if let Some(flags) = use_flags.filter(|&f| f != 0) {
            let mut n_terms = 0usize;
            for (bit, slot) in poly.iter_mut().enumerate() {
                if (flags >> bit) & 1 != 0 && n_terms < 6 {
                    *slot = calib[4 + n_terms];
                    max_order = Some(bit);
                    n_terms += 1;
                }
            }
        }
        let (lo, hi) = if calib[2] <= calib[3] {
            (calib[2], calib[3])
        } else {
            (calib[3], calib[2])
        };
        Ok(Self {
            coeff: calib[0],
            base: calib[1],
            lo,
            hi,
            poly,
            max_order,
        })
    }

    fn calibrate(&self, tof: f64) -> f64 {
        let traditional = (self.coeff * (tof - self.base)).powi(2);
        let Some(top) = self.max_order else {
            return traditional;
        };
        let clamped = tof.clamp(self.lo, self.hi);
        let mut correction = 0.0f64;
        for &c in self.poly[..=top].iter().rev() {
            correction = correction * clamped + c;
        }
        traditional - correction
    }
}

type SchemaMap = BTreeMap<String, Vec<(String, String)>>;
type CalibrationSet = (Vec<[f64; 10]>, BTreeMap<i64, u64>);
type DefaultMassCal = (BTreeMap<i64, [f64; 10]>, BTreeMap<i64, u64>);

/// Parses `MSScan.bin` (with `MSScan.xsd`, optional `MSTS.xml`, `MSMassCal.bin`, and
/// `DefaultMassCal.xml`) into ordered MS1 profile scans with per-scan calibration.
pub fn parse_profile_scans(
    msscan_bin: &[u8],
    msscan_xsd: &[u8],
    msts_xml: Option<&[u8]>,
    masscal_bin: Option<&[u8]>,
    default_masscal_xml: Option<&[u8]>,
) -> Result<Vec<MassHunterScan>> {
    let schema = parse_scan_xsd(msscan_xsd)?;
    let layout = resolve_scan_layout(&schema)?;
    let num_hint = msts_xml.and_then(count_scans_msts);
    let rec_offsets = scan_record_offsets(msscan_bin, &layout, num_hint)?;

    let mut calib_ids = Vec::with_capacity(rec_offsets.len());
    for &rec_pos in &rec_offsets {
        let cid = layout
            .calibration_id
            .map(|slot| slot.ty.read_i64(msscan_bin, rec_pos + slot.offset))
            .transpose()?;
        calib_ids.push(cid);
    }

    let (calibrations, flags_map) = load_calibration(masscal_bin, default_masscal_xml, &calib_ids)?;

    let mut scans = Vec::new();
    for (idx, &rec_pos) in rec_offsets.iter().enumerate() {
        let fmt_id = layout
            .spectrum_format_id
            .ty
            .read_i64(msscan_bin, rec_pos + layout.spectrum_format_id.offset)?;
        if fmt_id != PROFILE_FORMAT_ID {
            continue;
        }
        let ms_level = layout
            .ms_level
            .ty
            .read_i64(msscan_bin, rec_pos + layout.ms_level.offset)?;
        if ms_level != 1 {
            continue;
        }
        let polarity = match layout
            .ion_polarity
            .ty
            .read_i64(msscan_bin, rec_pos + layout.ion_polarity.offset)?
        {
            0 => "+".to_owned(),
            1 => "-".to_owned(),
            _ => continue,
        };
        let minutes = layout
            .scan_time
            .ty
            .read_f64(msscan_bin, rec_pos + layout.scan_time.offset)?;
        let tic = layout
            .tic
            .ty
            .read_f64(msscan_bin, rec_pos + layout.tic.offset)?;
        let offset = u64::try_from(
            layout
                .spectrum_offset
                .ty
                .read_i64(msscan_bin, rec_pos + layout.spectrum_offset.offset)?,
        )
        .map_err(|_| invalid("invalid SpectrumOffset"))?;
        let length = u64::try_from(
            layout
                .byte_count
                .ty
                .read_i64(msscan_bin, rec_pos + layout.byte_count.offset)?,
        )
        .map_err(|_| invalid("invalid ByteCount"))?;
        let points = u32::try_from(
            layout
                .point_count
                .ty
                .read_i64(msscan_bin, rec_pos + layout.point_count.offset)?,
        )
        .map_err(|_| invalid("invalid PointCount"))?;

        let calib = *calibrations
            .get(idx)
            .ok_or_else(|| invalid("calibration row missing for scan"))?;
        let flags = calib_ids[idx].and_then(|cid| flags_map.get(&cid).copied());

        scans.push(MassHunterScan {
            minutes,
            polarity,
            tic,
            offset,
            length,
            points,
            calibration: calib,
            flags,
        });
    }

    for pair in scans.windows(2) {
        let end = pair[0]
            .offset
            .checked_add(pair[0].length)
            .ok_or_else(|| invalid("profile scan offset overflow"))?;
        if pair[1].offset < end {
            return Err(invalid("MSScan.bin lists its profile scans out of order"));
        }
    }

    Ok(scans)
}

/// Decompresses a single run-length-encoded `MSProfile.bin` segment and calibrates its nonzero
/// `(m/z, intensity)` points.
pub fn decode_profile_spectrum(
    segment: &[u8],
    num_points: u32,
    calibration: &[f64; 10],
    use_flags: Option<u64>,
) -> Result<(Vec<f64>, Vec<f32>)> {
    if num_points > MAX_RLE_POINTS {
        return Err(invalid(
            "MSProfile.bin point count exceeds 24-bit RLE limit",
        ));
    }
    let expected_word = num_points | (RUN_LENGTH_MARKER << 24);
    let header_word = read_u32_le(segment, 16)
        .ok_or_else(|| invalid("an MSProfile.bin scan is not run-length encoded"))?;
    if header_word != expected_word {
        return Err(invalid("an MSProfile.bin scan is not run-length encoded"));
    }
    let first_tof =
        read_f64_le(segment, 0).ok_or_else(|| invalid("truncated MSProfile.bin segment header"))?;
    let step_tof =
        read_f64_le(segment, 8).ok_or_else(|| invalid("truncated MSProfile.bin segment header"))?;
    if !first_tof.is_finite() || !step_tof.is_finite() {
        return Err(invalid("non-finite TOF axis in MSProfile.bin segment"));
    }

    let calibrator = TofCalibrator::new(calibration, use_flags)?;
    let num_mz = num_points as usize;
    let comp = &segment[16..];
    let init_zero_repeat =
        read_i32_le(comp, 4).ok_or_else(|| invalid("Malformed MSProfile.bin RLE segment."))?;
    let cur_idx_i64 = -i64::from(init_zero_repeat);
    if cur_idx_i64 < 0 {
        return Err(invalid(
            "Malformed MSProfile.bin RLE segment: negative initial index.",
        ));
    }
    let mut cur_idx = usize::try_from(cur_idx_i64)
        .map_err(|_| invalid("Malformed MSProfile.bin RLE segment."))?;

    let mut mz_out = Vec::new();
    let mut inten_out = Vec::new();
    let mut offset = 8usize;
    let mut width: u8 = 4;

    while offset < comp.len() {
        let value: i64 = match width {
            1 => {
                let b = *comp
                    .get(offset)
                    .ok_or_else(|| invalid("Malformed MSProfile.bin RLE segment."))?;
                offset += 1;
                i64::from(b as i8)
            }
            2 => {
                let v = read_i16_le(comp, offset)
                    .ok_or_else(|| invalid("Malformed MSProfile.bin RLE segment."))?;
                offset += 2;
                i64::from(v)
            }
            4 => {
                let v = read_i32_le(comp, offset)
                    .ok_or_else(|| invalid("Malformed MSProfile.bin RLE segment."))?;
                offset += 4;
                i64::from(v)
            }
            _ => return Err(invalid("Malformed MSProfile.bin RLE segment.")),
        };

        if value >= 0 {
            if cur_idx >= num_mz {
                return Err(invalid("Malformed MSProfile.bin RLE segment."));
            }
            if value > 0 {
                let tof = first_tof + step_tof * (cur_idx as f64);
                mz_out.push(calibrator.calibrate(tof));
                inten_out.push(value as f32);
            }
            cur_idx += 1;
        } else {
            let mag = (-value) as u64;
            let num_zeros = usize::try_from(mag / 4)
                .map_err(|_| invalid("Malformed MSProfile.bin RLE segment."))?;
            width = match mag % 4 {
                1 => 1,
                2 => 2,
                3 => 4,
                _ => return Err(invalid("Malformed MSProfile.bin RLE segment.")),
            };
            cur_idx = cur_idx
                .checked_add(num_zeros)
                .ok_or_else(|| invalid("Malformed MSProfile.bin RLE segment."))?;
        }
    }

    Ok((mz_out, inten_out))
}

fn parse_scan_xsd(xsd: &[u8]) -> Result<SchemaMap> {
    let text = decode_text_auto(xsd);
    let doc = parse_xml(&text).map_err(|e| InstrumentError::Xml(e.to_string()))?;
    let mut schema = SchemaMap::new();
    for ct in elements(doc.root(), "complexType") {
        let Some(name) = ct.attribute("name") else {
            continue;
        };
        let mut fields = Vec::new();
        for el in elements(ct, "element") {
            if let (Some(fname), Some(ftype)) = (el.attribute("name"), el.attribute("type")) {
                fields.push((fname.to_owned(), ftype.to_owned()));
            }
        }
        schema.insert(name.to_owned(), fields);
    }
    Ok(schema)
}

fn count_scans_msts(xml: &[u8]) -> Option<usize> {
    let text = decode_text_auto(xml);
    let doc = parse_xml(&text).ok()?;
    let mut sum = 0usize;
    let mut any = false;
    for node in elements(doc.root(), "NumOfScans") {
        if let Some(txt) = trimmed_text(node)
            && let Ok(n) = txt.parse::<usize>()
        {
            sum = sum.saturating_add(n);
            any = true;
        }
    }
    any.then_some(sum)
}

fn flatten_fields(
    schema: &SchemaMap,
    type_name: &str,
    prefix: &str,
    base_offset: usize,
    depth: usize,
    out: &mut BTreeMap<String, FieldSlot>,
) -> Result<usize> {
    if depth > 16 {
        return Err(invalid("cyclic XSD type definition"));
    }
    if let Some(simple) = XsdType::from_name(type_name) {
        if !prefix.is_empty() {
            out.insert(
                prefix.to_owned(),
                FieldSlot {
                    offset: base_offset,
                    ty: simple,
                },
            );
        }
        return Ok(simple.byte_size());
    }
    let short = type_name.rsplit(':').next().unwrap_or(type_name);
    let fields = schema
        .get(short)
        .ok_or_else(|| invalid(format!("unknown XSD type {type_name}")))?;
    let mut off = base_offset;
    for (fname, ftype) in fields {
        let key = if prefix.is_empty() {
            fname.clone()
        } else {
            format!("{prefix}.{fname}")
        };
        let size = flatten_fields(schema, ftype, &key, off, depth + 1, out)?;
        off = off
            .checked_add(size)
            .ok_or_else(|| invalid("XSD type size overflow"))?;
    }
    Ok(off - base_offset)
}

fn resolve_scan_layout(schema: &SchemaMap) -> Result<ScanLayout> {
    let mut dummy = BTreeMap::new();
    let block_size = flatten_fields(schema, "SpectrumParamsType", "", 0, 0, &mut dummy)?;
    let mut slots = BTreeMap::new();
    let full_rec = flatten_fields(schema, "ScanRecordType", "", 0, 0, &mut slots)?;
    let scalar_size = full_rec
        .checked_sub(block_size)
        .ok_or_else(|| invalid("ScanRecordType smaller than SpectrumParamsType"))?;

    let req = |name: &str| -> Result<FieldSlot> {
        slots
            .get(name)
            .copied()
            .ok_or_else(|| invalid(format!("missing XSD field {name}")))
    };

    Ok(ScanLayout {
        full_rec,
        block_size,
        scalar_size,
        scan_time: req("ScanTime")?,
        ion_polarity: req("IonPolarity")?,
        ms_level: req("MSLevel")?,
        tic: req("TIC")?,
        calibration_id: slots.get("CalibrationID").copied(),
        spectrum_format_id: req("SpectrumParamValues.SpectrumFormatID")?,
        spectrum_offset: req("SpectrumParamValues.SpectrumOffset")?,
        byte_count: req("SpectrumParamValues.ByteCount")?,
        point_count: req("SpectrumParamValues.PointCount")?,
    })
}

fn scan_record_offsets(
    msscan_bin: &[u8],
    layout: &ScanLayout,
    num_hint: Option<usize>,
) -> Result<Vec<usize>> {
    let rec_start = read_u32_le(msscan_bin, 0x58)
        .ok_or_else(|| invalid("truncated MSScan.bin header"))? as usize;
    if rec_start > msscan_bin.len() {
        return Err(invalid("MSScan.bin record start past EOF"));
    }
    let body = msscan_bin.len() - rec_start;

    let candidates: Vec<usize> = (1..=MAX_SPECTRUM_BLOCKS)
        .filter_map(|n| {
            layout
                .scalar_size
                .checked_add(n.checked_mul(layout.block_size)?)
        })
        .collect();

    let mut stride: Option<usize> = None;
    if let Some(n) = num_hint
        && n > 0
        && body.is_multiple_of(n)
    {
        let hinted = body / n;
        if candidates.contains(&hinted) {
            stride = Some(hinted);
        }
    }
    if stride.is_none() {
        let divisors: Vec<usize> = candidates
            .iter()
            .copied()
            .filter(|&s| s > 0 && body.is_multiple_of(s))
            .collect();
        if !divisors.is_empty() {
            stride = divisors.into_iter().min_by_key(|&s| match num_hint {
                Some(n) => (body / s).abs_diff(n),
                None => s,
            });
        }
    }

    let step = stride.unwrap_or(layout.full_rec);
    if step == 0 || layout.full_rec == 0 {
        return Err(invalid("zero MSScan.bin record stride"));
    }

    let mut offsets = Vec::with_capacity(body / step);
    let mut pos = rec_start;
    while let Some(end) = pos.checked_add(layout.full_rec)
        && end <= msscan_bin.len()
    {
        offsets.push(pos);
        let Some(next) = pos.checked_add(step) else {
            break;
        };
        pos = next;
    }

    Ok(offsets)
}

fn load_calibration(
    masscal_bin: Option<&[u8]>,
    default_masscal_xml: Option<&[u8]>,
    calib_ids: &[Option<i64>],
) -> Result<CalibrationSet> {
    let (default_rows, flags) = default_masscal_xml
        .map(parse_default_masscal)
        .unwrap_or_default();

    if let Some(bin) = masscal_bin {
        let mut rows = Vec::with_capacity(calib_ids.len());
        for i in 0..calib_ids.len() {
            let base = 0x4Cusize
                .checked_add(
                    i.checked_mul(84)
                        .ok_or_else(|| invalid("MSMassCal offset overflow"))?,
                )
                .ok_or_else(|| invalid("MSMassCal offset overflow"))?;
            let mut row = [0.0f64; 10];
            for (j, slot) in row.iter_mut().enumerate() {
                let off = base
                    .checked_add(j * 8)
                    .ok_or_else(|| invalid("MSMassCal offset overflow"))?;
                *slot = read_f64_le(bin, off).ok_or_else(|| invalid("truncated MSMassCal.bin"))?;
            }
            rows.push(row);
        }
        return Ok((rows, flags));
    }

    let Some((_, &fallback)) = default_rows.iter().next() else {
        return Err(invalid("no mass calibration"));
    };
    let rows = calib_ids
        .iter()
        .map(|cid| {
            cid.and_then(|id| default_rows.get(&id).copied())
                .unwrap_or(fallback)
        })
        .collect();
    Ok((rows, flags))
}

fn parse_default_masscal(xml: &[u8]) -> DefaultMassCal {
    let mut rows = BTreeMap::new();
    let mut flags = BTreeMap::new();
    let text = decode_text_auto(xml);
    let Ok(doc) = parse_xml(&text) else {
        return (rows, flags);
    };
    for cal in elements(doc.root(), "DefaultCalibration") {
        let Some(cid) = cal
            .attribute("DefaultCalibrationID")
            .and_then(|s| s.trim().parse::<i64>().ok())
        else {
            continue;
        };
        let mut traditional = Vec::new();
        let mut polynomial = Vec::new();
        for step in elements(cal, "Step") {
            let formula = elements(step, "CalibrationFormula")
                .next()
                .and_then(trimmed_text);
            match formula.as_deref() {
                Some("Traditional") => {
                    traditional = elements(step, "Value")
                        .filter_map(trimmed_text)
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect();
                }
                Some("Polynomial") => {
                    polynomial = elements(step, "Value")
                        .filter_map(trimmed_text)
                        .filter_map(|s| s.parse::<f64>().ok())
                        .collect();
                    let use_flags = elements(step, "ValueUseFlags")
                        .next()
                        .and_then(trimmed_text)
                        .and_then(|s| s.parse::<u64>().ok())
                        .unwrap_or(0);
                    flags.insert(cid, use_flags);
                }
                _ => {}
            }
        }
        if traditional.len() >= 2 {
            let mut row = [0.0f64; 10];
            row[0] = traditional[0];
            row[1] = traditional[1];
            for (j, &v) in polynomial.iter().take(8).enumerate() {
                row[2 + j] = v;
            }
            rows.insert(cid, row);
        }
    }
    (rows, flags)
}

fn invalid(msg: impl Into<String>) -> InstrumentError {
    InstrumentError::InvalidMassHunter(msg.into())
}
