//! Pure-Rust netCDF-3 (`CDF\x01` and `CDF\x02`) binary container parser.
//!
//! Reads classic 32-bit offset and 64-bit offset netCDF-3 files with zero unsafe code and zero C
//! dependencies, exposing dimensions, global and variable attributes, and typed variable readers.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::binary::{read_f32_be, read_f64_be, read_i16_be, read_i32_be, read_u32_be, read_u64_be};
use crate::error::{InstrumentError, Result};

const MAGIC_CDF1: &[u8; 4] = b"CDF\x01";
const MAGIC_CDF2: &[u8; 4] = b"CDF\x02";

const TAG_ABSENT: u32 = 0;
const TAG_DIMENSION: u32 = 10;
const TAG_VARIABLE: u32 = 11;
const TAG_ATTRIBUTE: u32 = 12;

const STREAMING_NUMRECS: u32 = 0xFFFF_FFFF;

/// NetCDF-3 default fill threshold for floating-point values (`9.969209968386869e36`).
const NC_FILL_FLOAT_THRESHOLD: f64 = 9.9e36;

pub(crate) fn custom_fill_value(var: &NcVariable) -> Option<f64> {
    var.attrs
        .get("_FillValue")
        .or_else(|| var.attrs.get("missing_value"))
        .and_then(NcAttrValue::as_f64)
}

pub(crate) fn is_valid_number(val: f64, fill: Option<f64>) -> bool {
    if !val.is_finite() || val.abs() >= NC_FILL_FLOAT_THRESHOLD {
        return false;
    }
    if let Some(f) = fill
        && (val - f).abs() <= f64::EPSILON
    {
        return false;
    }
    true
}

/// NetCDF-3 binary format version.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NetCdfVersion {
    /// Classic 32-bit offset format (`CDF\x01`).
    Classic,
    /// 64-bit offset format (`CDF\x02`).
    Offset64,
}

/// NetCDF-3 primitive data type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NcType {
    /// Signed 8-bit integer (`NC_BYTE = 1`).
    Byte,
    /// 8-bit character (`NC_CHAR = 2`).
    Char,
    /// Signed 16-bit big-endian integer (`NC_SHORT = 3`).
    Short,
    /// Signed 32-bit big-endian integer (`NC_INT = 4`).
    Int,
    /// 32-bit big-endian IEEE 754 float (`NC_FLOAT = 5`).
    Float,
    /// 64-bit big-endian IEEE 754 double (`NC_DOUBLE = 6`).
    Double,
}

impl NcType {
    fn from_code(code: u32) -> Result<Self> {
        match code {
            1 => Ok(Self::Byte),
            2 => Ok(Self::Char),
            3 => Ok(Self::Short),
            4 => Ok(Self::Int),
            5 => Ok(Self::Float),
            6 => Ok(Self::Double),
            other => Err(invalid(format!("unknown netCDF type code {other}"))),
        }
    }

    const fn byte_size(self) -> usize {
        match self {
            Self::Byte | Self::Char => 1,
            Self::Short => 2,
            Self::Int | Self::Float => 4,
            Self::Double => 8,
        }
    }
}

/// Typed value stored in a netCDF-3 global or variable attribute.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum NcAttrValue {
    /// `NC_BYTE` values.
    Bytes(Vec<i8>),
    /// `NC_CHAR` text decoded as UTF-8 with trailing NULs removed.
    Chars(String),
    /// `NC_SHORT` values.
    Shorts(Vec<i16>),
    /// `NC_INT` values.
    Ints(Vec<i32>),
    /// `NC_FLOAT` values.
    Floats(Vec<f32>),
    /// `NC_DOUBLE` values.
    Doubles(Vec<f64>),
}

impl NcAttrValue {
    /// Returns the attribute as trimmed text when it is `NC_CHAR` and non-empty.
    #[must_use]
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Chars(text) => {
                let trimmed = text.trim();
                if trimmed.is_empty() {
                    None
                } else {
                    Some(trimmed)
                }
            }
            _ => None,
        }
    }

    /// Returns the first finite numeric value stored in the attribute.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        let val = match self {
            Self::Bytes(v) => v.first().map(|&x| f64::from(x)),
            Self::Shorts(v) => v.first().map(|&x| f64::from(x)),
            Self::Ints(v) => v.first().map(|&x| f64::from(x)),
            Self::Floats(v) => v.first().map(|&x| f64::from(x)),
            Self::Doubles(v) => v.first().copied(),
            Self::Chars(_) => None,
        }?;
        val.is_finite().then_some(val)
    }
}

/// Dimension declared in a netCDF-3 header.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NcDimension {
    /// Dimension name.
    pub name: String,
    /// Fixed dimension length, or `None` for the unlimited (record) dimension.
    pub length: Option<u32>,
}

/// Variable descriptor declared in a netCDF-3 header.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NcVariable {
    /// Variable name.
    pub name: String,
    /// Indices into the file's dimension table.
    pub dim_ids: Vec<usize>,
    /// Variable-scoped attributes.
    pub attrs: BTreeMap<String, NcAttrValue>,
    /// Element storage type.
    pub nc_type: NcType,
    /// Header-declared byte size (total for non-record variables, per-record slice for record variables).
    pub vsize: u32,
    /// Byte offset of the variable's first element from the start of the file.
    pub begin: u64,
    /// Whether the variable's leading dimension is the unlimited dimension.
    pub is_record: bool,
}

/// Parsed netCDF-3 binary container borrowing the original input buffer.
#[derive(Debug, Clone)]
pub struct NetCdfFile<'a> {
    version: NetCdfVersion,
    numrecs: usize,
    dims: Vec<NcDimension>,
    global_attrs: BTreeMap<String, NcAttrValue>,
    vars: Vec<NcVariable>,
    recsize: usize,
    data: &'a [u8],
}

impl<'a> NetCdfFile<'a> {
    /// Parses a netCDF-3 classic (`CDF\x01`) or 64-bit offset (`CDF\x02`) header and validates layout.
    pub fn parse(data: &'a [u8]) -> Result<Self> {
        let mut cursor = Cursor::new(data);
        let magic = cursor.read_bytes(4)?;
        let version = if magic == MAGIC_CDF1 {
            NetCdfVersion::Classic
        } else if magic == MAGIC_CDF2 {
            NetCdfVersion::Offset64
        } else {
            return Err(invalid("missing CDF\\x01 or CDF\\x02 magic header"));
        };

        let raw_numrecs = cursor.read_u32()?;
        let dims = parse_dimensions(&mut cursor)?;
        let global_attrs = parse_attributes(&mut cursor)?;
        let vars = parse_variables(&mut cursor, version, &dims)?;

        let mut recsize: usize = 0;
        let mut min_rec_begin: Option<usize> = None;
        for var in &vars {
            if !var.is_record {
                continue;
            }
            let slice_elems = product_dim_lengths(&dims, &var.dim_ids[1..])?;
            let slice_bytes = slice_elems
                .checked_mul(var.nc_type.byte_size())
                .ok_or_else(|| {
                    invalid(format!(
                        "variable '{}' record slice size overflow",
                        var.name
                    ))
                })?;
            let stride = if var.vsize == u32::MAX {
                round_up_4(slice_bytes)?
            } else {
                let declared = var.vsize as usize;
                if declared < slice_bytes {
                    return Err(invalid(format!(
                        "variable '{}' vsize ({declared}) smaller than record slice ({slice_bytes})",
                        var.name
                    )));
                }
                declared
            };
            recsize = recsize
                .checked_add(stride)
                .ok_or_else(|| invalid("record size overflow"))?;
            let begin = usize::try_from(var.begin).map_err(|_| {
                invalid(format!(
                    "variable '{}' offset exceeds addressable range",
                    var.name
                ))
            })?;
            min_rec_begin = Some(match min_rec_begin {
                Some(prev) => prev.min(begin),
                None => begin,
            });
        }

        let numrecs = if raw_numrecs == STREAMING_NUMRECS {
            match (min_rec_begin, recsize) {
                (Some(begin), stride) if stride > 0 => {
                    if begin > data.len() {
                        return Err(invalid("streaming record offset beyond end of file"));
                    }
                    (data.len() - begin) / stride
                }
                _ => 0,
            }
        } else {
            raw_numrecs as usize
        };

        Ok(Self {
            version,
            numrecs,
            dims,
            global_attrs,
            vars,
            recsize,
            data,
        })
    }

    /// NetCDF-3 format version.
    #[must_use]
    pub const fn version(&self) -> NetCdfVersion {
        self.version
    }

    /// Number of records in the unlimited dimension.
    #[must_use]
    pub const fn numrecs(&self) -> usize {
        self.numrecs
    }

    /// Declared dimensions.
    #[must_use]
    pub fn dimensions(&self) -> &[NcDimension] {
        &self.dims
    }

    /// Global attributes by name.
    #[must_use]
    pub fn global_attrs(&self) -> &BTreeMap<String, NcAttrValue> {
        &self.global_attrs
    }

    /// Declared variables.
    #[must_use]
    pub fn variables(&self) -> &[NcVariable] {
        &self.vars
    }

    /// Looks up a global string attribute by name.
    #[must_use]
    pub fn global_str(&self, name: &str) -> Option<&str> {
        self.global_attrs.get(name).and_then(NcAttrValue::as_str)
    }

    /// Looks up a variable descriptor by name.
    #[must_use]
    pub fn variable(&self, name: &str) -> Option<&NcVariable> {
        self.vars.iter().find(|v| v.name == name)
    }

    /// Reads all elements of a numeric variable as `f64`, supporting both record and non-record layouts.
    pub fn read_var_f64(&self, name: &str) -> Result<Option<Vec<f64>>> {
        let Some(var) = self.variable(name) else {
            return Ok(None);
        };
        if var.nc_type == NcType::Char {
            return Err(invalid(format!(
                "variable '{name}' is NC_CHAR, expected numeric"
            )));
        }

        let elem_size = var.nc_type.byte_size();
        if var.is_record {
            let inner = product_dim_lengths(&self.dims, &var.dim_ids[1..])?;
            let slice_bytes = inner
                .checked_mul(elem_size)
                .ok_or_else(|| invalid(format!("variable '{name}' byte length overflow")))?;
            let total_elems = self
                .numrecs
                .checked_mul(inner)
                .ok_or_else(|| invalid(format!("variable '{name}' element count overflow")))?;
            if total_elems == 0 {
                return Ok(Some(Vec::new()));
            }
            let begin = self.checked_record_bounds(var, slice_bytes)?;
            let mut out = Vec::with_capacity(total_elems);
            for rec in 0..self.numrecs {
                let offset = begin + rec * self.recsize;
                decode_numeric_slice(
                    &self.data[offset..offset + slice_bytes],
                    var.nc_type,
                    inner,
                    &mut out,
                )?;
            }
            Ok(Some(out))
        } else {
            let total_elems = product_dim_lengths(&self.dims, &var.dim_ids)?;
            let total_bytes = total_elems
                .checked_mul(elem_size)
                .ok_or_else(|| invalid(format!("variable '{name}' byte length overflow")))?;
            let slice = self.non_record_slice(var, total_bytes)?;
            let mut out = Vec::with_capacity(total_elems);
            decode_numeric_slice(slice, var.nc_type, total_elems, &mut out)?;
            Ok(Some(out))
        }
    }

    /// Reads the first finite non-fill numeric value of a variable.
    pub fn read_var_scalar_f64(&self, name: &str) -> Result<Option<f64>> {
        let Some(var) = self.variable(name) else {
            return Ok(None);
        };
        let fill = custom_fill_value(var);
        let Some(values) = self.read_var_f64(name)? else {
            return Ok(None);
        };
        Ok(values
            .into_iter()
            .next()
            .filter(|&v| is_valid_number(v, fill)))
    }

    /// Reads a 1D `NC_CHAR` variable as a trimmed string (truncated at the first NUL byte).
    pub fn read_var_string(&self, name: &str) -> Result<Option<String>> {
        let Some(var) = self.variable(name) else {
            return Ok(None);
        };
        if var.nc_type != NcType::Char {
            return Err(invalid(format!("variable '{name}' is not NC_CHAR")));
        }
        let raw = self.read_var_bytes(var)?;
        let text = decode_c_string(&raw);
        if text.is_empty() {
            Ok(None)
        } else {
            Ok(Some(text))
        }
    }

    /// Reads a 2D `NC_CHAR` variable (array of fixed-width strings) into a vector of strings.
    pub fn read_var_strings(&self, name: &str) -> Result<Option<Vec<String>>> {
        let Some(var) = self.variable(name) else {
            return Ok(None);
        };
        if var.nc_type != NcType::Char {
            return Err(invalid(format!("variable '{name}' is not NC_CHAR")));
        }
        if var.dim_ids.len() == 1 {
            return Ok(Some(vec![decode_c_string(&self.read_var_bytes(var)?)]));
        }
        if var.dim_ids.len() != 2 {
            return Err(invalid(format!(
                "variable '{name}' has rank {}, expected 1 or 2 for string array",
                var.dim_ids.len()
            )));
        }
        let rows = self.dim_len(var.dim_ids[0]);
        let row_len = self.dim_len(var.dim_ids[1]);
        if rows == 0 || row_len == 0 {
            return Ok(Some(Vec::new()));
        }
        let raw = self.read_var_bytes(var)?;
        let expected = rows
            .checked_mul(row_len)
            .ok_or_else(|| invalid(format!("variable '{name}' string array size overflow")))?;
        if raw.len() < expected {
            return Err(invalid(format!("variable '{name}' truncated string data")));
        }
        let mut out = Vec::with_capacity(rows);
        for chunk in raw[..expected].chunks_exact(row_len) {
            out.push(decode_c_string(chunk));
        }
        Ok(Some(out))
    }

    fn dim_len(&self, dim_id: usize) -> usize {
        self.dims[dim_id]
            .length
            .map_or(self.numrecs, |len| len as usize)
    }

    fn read_var_bytes(&self, var: &NcVariable) -> Result<Vec<u8>> {
        let elem_size = var.nc_type.byte_size();
        if var.is_record {
            let inner = product_dim_lengths(&self.dims, &var.dim_ids[1..])?;
            let slice_bytes = inner
                .checked_mul(elem_size)
                .ok_or_else(|| invalid(format!("variable '{}' byte length overflow", var.name)))?;
            let total_bytes = self.numrecs.checked_mul(slice_bytes).ok_or_else(|| {
                invalid(format!(
                    "variable '{}' total byte length overflow",
                    var.name
                ))
            })?;
            if total_bytes == 0 {
                return Ok(Vec::new());
            }
            let begin = self.checked_record_bounds(var, slice_bytes)?;
            let mut out = Vec::with_capacity(total_bytes);
            for rec in 0..self.numrecs {
                let offset = begin + rec * self.recsize;
                out.extend_from_slice(&self.data[offset..offset + slice_bytes]);
            }
            Ok(out)
        } else {
            let total_elems = product_dim_lengths(&self.dims, &var.dim_ids)?;
            let total_bytes = total_elems
                .checked_mul(elem_size)
                .ok_or_else(|| invalid(format!("variable '{}' byte length overflow", var.name)))?;
            Ok(self.non_record_slice(var, total_bytes)?.to_vec())
        }
    }

    fn non_record_slice(&self, var: &NcVariable, total_bytes: usize) -> Result<&'a [u8]> {
        let begin = usize::try_from(var.begin).map_err(|_| {
            invalid(format!(
                "variable '{}' offset exceeds addressable range",
                var.name
            ))
        })?;
        let end = begin
            .checked_add(total_bytes)
            .ok_or_else(|| invalid(format!("variable '{}' end offset overflow", var.name)))?;
        self.data.get(begin..end).ok_or_else(|| {
            invalid(format!(
                "variable '{}' data range {begin}..{end} exceeds file length {}",
                var.name,
                self.data.len()
            ))
        })
    }

    fn checked_record_bounds(&self, var: &NcVariable, slice_bytes: usize) -> Result<usize> {
        let begin = usize::try_from(var.begin).map_err(|_| {
            invalid(format!(
                "variable '{}' offset exceeds addressable range",
                var.name
            ))
        })?;
        if self.numrecs == 0 || slice_bytes == 0 {
            return Ok(begin);
        }
        let last_stride = (self.numrecs - 1)
            .checked_mul(self.recsize)
            .ok_or_else(|| invalid(format!("variable '{}' record stride overflow", var.name)))?;
        let last_start = begin
            .checked_add(last_stride)
            .ok_or_else(|| invalid(format!("variable '{}' record offset overflow", var.name)))?;
        let last_end = last_start
            .checked_add(slice_bytes)
            .ok_or_else(|| invalid(format!("variable '{}' record end overflow", var.name)))?;
        if last_end > self.data.len() {
            return Err(invalid(format!(
                "variable '{}' record data ends at {last_end}, exceeding file length {}",
                var.name,
                self.data.len()
            )));
        }
        Ok(begin)
    }
}

fn parse_dimensions(cursor: &mut Cursor<'_>) -> Result<Vec<NcDimension>> {
    let tag = cursor.read_u32()?;
    let count = cursor.read_u32()? as usize;
    if tag == TAG_ABSENT {
        if count != 0 {
            return Err(invalid("ABSENT dimension list has non-zero count"));
        }
        return Ok(Vec::new());
    }
    if tag != TAG_DIMENSION {
        return Err(invalid(format!("unexpected dimension list tag {tag}")));
    }
    cursor.require_at_least(count, 8, "dimension list")?;

    let mut dims = Vec::with_capacity(count);
    let mut saw_unlimited = false;
    for _ in 0..count {
        let name = cursor.read_name()?;
        let raw_len = cursor.read_u32()?;
        let length = if raw_len == 0 {
            if saw_unlimited {
                return Err(invalid("multiple unlimited dimensions declared"));
            }
            saw_unlimited = true;
            None
        } else {
            Some(raw_len)
        };
        dims.push(NcDimension { name, length });
    }
    Ok(dims)
}

fn parse_attributes(cursor: &mut Cursor<'_>) -> Result<BTreeMap<String, NcAttrValue>> {
    let tag = cursor.read_u32()?;
    let count = cursor.read_u32()? as usize;
    if tag == TAG_ABSENT {
        if count != 0 {
            return Err(invalid("ABSENT attribute list has non-zero count"));
        }
        return Ok(BTreeMap::new());
    }
    if tag != TAG_ATTRIBUTE {
        return Err(invalid(format!("unexpected attribute list tag {tag}")));
    }
    cursor.require_at_least(count, 12, "attribute list")?;

    let mut attrs = BTreeMap::new();
    for _ in 0..count {
        let name = cursor.read_name()?;
        let nc_type = NcType::from_code(cursor.read_u32()?)?;
        let nelems = cursor.read_u32()? as usize;
        let byte_len = nelems
            .checked_mul(nc_type.byte_size())
            .ok_or_else(|| invalid(format!("attribute '{name}' byte length overflow")))?;
        let raw = cursor.read_padded_bytes(byte_len)?;
        let value = decode_attr_value(raw, nc_type, nelems)?;
        attrs.insert(name, value);
    }
    Ok(attrs)
}

fn parse_variables(
    cursor: &mut Cursor<'_>,
    version: NetCdfVersion,
    dims: &[NcDimension],
) -> Result<Vec<NcVariable>> {
    let tag = cursor.read_u32()?;
    let count = cursor.read_u32()? as usize;
    if tag == TAG_ABSENT {
        if count != 0 {
            return Err(invalid("ABSENT variable list has non-zero count"));
        }
        return Ok(Vec::new());
    }
    if tag != TAG_VARIABLE {
        return Err(invalid(format!("unexpected variable list tag {tag}")));
    }
    cursor.require_at_least(count, 24, "variable list")?;

    let mut vars = Vec::with_capacity(count);
    for _ in 0..count {
        let name = cursor.read_name()?;
        let rank = cursor.read_u32()? as usize;
        cursor.require_at_least(rank, 4, "variable dim_ids")?;
        let mut dim_ids = Vec::with_capacity(rank);
        for dim_pos in 0..rank {
            let dim_id = cursor.read_u32()? as usize;
            let Some(dim) = dims.get(dim_id) else {
                return Err(invalid(format!(
                    "variable '{name}' references out-of-range dimension id {dim_id}"
                )));
            };
            if dim_pos > 0 && dim.length.is_none() {
                return Err(invalid(format!(
                    "variable '{name}' uses unlimited dimension at non-leading index {dim_pos}"
                )));
            }
            dim_ids.push(dim_id);
        }
        let attrs = parse_attributes(cursor)?;
        let nc_type = NcType::from_code(cursor.read_u32()?)?;
        let vsize = cursor.read_u32()?;
        let begin = match version {
            NetCdfVersion::Classic => u64::from(cursor.read_u32()?),
            NetCdfVersion::Offset64 => cursor.read_u64()?,
        };
        let is_record = dim_ids.first().is_some_and(|&id| dims[id].length.is_none());
        vars.push(NcVariable {
            name,
            dim_ids,
            attrs,
            nc_type,
            vsize,
            begin,
            is_record,
        });
    }
    Ok(vars)
}

fn decode_attr_value(raw: &[u8], nc_type: NcType, nelems: usize) -> Result<NcAttrValue> {
    match nc_type {
        NcType::Byte => Ok(NcAttrValue::Bytes(
            raw.iter().map(|&b| i8::from_be_bytes([b])).collect(),
        )),
        NcType::Char => Ok(NcAttrValue::Chars(decode_c_string(raw))),
        NcType::Short => {
            let mut out = Vec::with_capacity(nelems);
            for i in 0..nelems {
                let val = read_i16_be(raw, i * 2)
                    .ok_or_else(|| invalid("truncated NC_SHORT attribute"))?;
                out.push(val);
            }
            Ok(NcAttrValue::Shorts(out))
        }
        NcType::Int => {
            let mut out = Vec::with_capacity(nelems);
            for i in 0..nelems {
                let val =
                    read_i32_be(raw, i * 4).ok_or_else(|| invalid("truncated NC_INT attribute"))?;
                out.push(val);
            }
            Ok(NcAttrValue::Ints(out))
        }
        NcType::Float => {
            let mut out = Vec::with_capacity(nelems);
            for i in 0..nelems {
                let val = read_f32_be(raw, i * 4)
                    .ok_or_else(|| invalid("truncated NC_FLOAT attribute"))?;
                out.push(val);
            }
            Ok(NcAttrValue::Floats(out))
        }
        NcType::Double => {
            let mut out = Vec::with_capacity(nelems);
            for i in 0..nelems {
                let val = read_f64_be(raw, i * 8)
                    .ok_or_else(|| invalid("truncated NC_DOUBLE attribute"))?;
                out.push(val);
            }
            Ok(NcAttrValue::Doubles(out))
        }
    }
}

fn decode_numeric_slice(
    raw: &[u8],
    nc_type: NcType,
    nelems: usize,
    out: &mut Vec<f64>,
) -> Result<()> {
    match nc_type {
        NcType::Byte => {
            for &b in raw.iter().take(nelems) {
                out.push(f64::from(i8::from_be_bytes([b])));
            }
        }
        NcType::Short => {
            for i in 0..nelems {
                let val = read_i16_be(raw, i * 2)
                    .ok_or_else(|| invalid("truncated NC_SHORT variable data"))?;
                out.push(f64::from(val));
            }
        }
        NcType::Int => {
            for i in 0..nelems {
                let val = read_i32_be(raw, i * 4)
                    .ok_or_else(|| invalid("truncated NC_INT variable data"))?;
                out.push(f64::from(val));
            }
        }
        NcType::Float => {
            for i in 0..nelems {
                let val = read_f32_be(raw, i * 4)
                    .ok_or_else(|| invalid("truncated NC_FLOAT variable data"))?;
                out.push(f64::from(val));
            }
        }
        NcType::Double => {
            for i in 0..nelems {
                let val = read_f64_be(raw, i * 8)
                    .ok_or_else(|| invalid("truncated NC_DOUBLE variable data"))?;
                out.push(val);
            }
        }
        NcType::Char => {
            return Err(invalid("cannot decode NC_CHAR as numeric"));
        }
    }
    Ok(())
}

fn product_dim_lengths(dims: &[NcDimension], dim_ids: &[usize]) -> Result<usize> {
    let mut total = 1usize;
    for &id in dim_ids {
        let len =
            dims.get(id).and_then(|d| d.length).ok_or_else(|| {
                invalid("unexpected unlimited or invalid dimension in fixed shape")
            })? as usize;
        total = total
            .checked_mul(len)
            .ok_or_else(|| invalid("dimension product overflow"))?;
    }
    Ok(total)
}

fn round_up_4(n: usize) -> Result<usize> {
    let rem = n % 4;
    if rem == 0 {
        Ok(n)
    } else {
        n.checked_add(4 - rem)
            .ok_or_else(|| invalid("4-byte alignment overflow"))
    }
}

fn decode_c_string(bytes: &[u8]) -> String {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let slice = &bytes[..end];
    if let Ok(utf8) = std::str::from_utf8(slice) {
        return utf8.trim().to_owned();
    }
    let (decoded, _, _) = encoding_rs::GB18030.decode(slice);
    decoded.trim().to_owned()
}

pub(crate) fn invalid(msg: impl Into<String>) -> InstrumentError {
    InstrumentError::InvalidAndi(msg.into())
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    const fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    const fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    fn require_at_least(&self, count: usize, min_elem_bytes: usize, label: &str) -> Result<()> {
        let min_bytes = count
            .checked_mul(min_elem_bytes)
            .ok_or_else(|| invalid(format!("{label} count overflow")))?;
        if min_bytes > self.remaining() {
            return Err(invalid(format!(
                "{label} declares {count} items, exceeding remaining {} bytes",
                self.remaining()
            )));
        }
        Ok(())
    }

    fn read_bytes(&mut self, len: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(len)
            .ok_or_else(|| invalid("header cursor overflow"))?;
        let slice = self.data.get(self.pos..end).ok_or_else(|| {
            invalid(format!(
                "unexpected end of file at offset {} (needed {len} bytes)",
                self.pos
            ))
        })?;
        self.pos = end;
        Ok(slice)
    }

    fn read_padded_bytes(&mut self, len: usize) -> Result<&'a [u8]> {
        let padded = round_up_4(len)?;
        let full = self.read_bytes(padded)?;
        Ok(&full[..len])
    }

    fn read_u32(&mut self) -> Result<u32> {
        let val = read_u32_be(self.data, self.pos).ok_or_else(|| {
            invalid(format!(
                "unexpected end of file reading u32 at offset {}",
                self.pos
            ))
        })?;
        self.pos += 4;
        Ok(val)
    }

    fn read_u64(&mut self) -> Result<u64> {
        let val = read_u64_be(self.data, self.pos).ok_or_else(|| {
            invalid(format!(
                "unexpected end of file reading u64 at offset {}",
                self.pos
            ))
        })?;
        self.pos += 8;
        Ok(val)
    }

    fn read_name(&mut self) -> Result<String> {
        let len = self.read_u32()? as usize;
        let bytes = self.read_padded_bytes(len)?;
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }
}
