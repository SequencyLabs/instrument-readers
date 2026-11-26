//! Integration tests for the pure-Rust ANDI Chromatography netCDF-3 (`.cdf`) reader.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use sequency_instruments::andi::{
    NcAttrValue, NcType, NetCdfFile, NetCdfVersion, parse_cdf, read_cdf,
};
use sequency_instruments::{InstrumentError, Technique, VendorFormat, detect, read, read_files};

const TAG_ABSENT: u32 = 0;
const TAG_DIMENSION: u32 = 10;
const TAG_VARIABLE: u32 = 11;
const TAG_ATTRIBUTE: u32 = 12;

enum AttrSpec<'a> {
    Bytes(&'a [i8]),
    Chars(&'a str),
    Shorts(&'a [i16]),
    Ints(&'a [i32]),
    Floats(&'a [f32]),
    Doubles(&'a [f64]),
}

struct VarSpec<'a> {
    name: &'a str,
    dim_ids: Vec<u32>,
    attrs: Vec<(&'a str, AttrSpec<'a>)>,
    nc_type: NcType,
    /// Non-record variables store all raw big-endian bytes in `slices[0]`.
    /// Record variables store one raw big-endian slice per record in `slices`.
    slices: Vec<Vec<u8>>,
    is_record: bool,
}

struct NetCdfBuilder<'a> {
    version: NetCdfVersion,
    numrecs: u32,
    dims: Vec<(&'a str, u32)>,
    global_attrs: Vec<(&'a str, AttrSpec<'a>)>,
    vars: Vec<VarSpec<'a>>,
}

impl<'a> NetCdfBuilder<'a> {
    fn new(version: NetCdfVersion) -> Self {
        Self {
            version,
            numrecs: 0,
            dims: Vec::new(),
            global_attrs: Vec::new(),
            vars: Vec::new(),
        }
    }

    fn numrecs(mut self, numrecs: u32) -> Self {
        self.numrecs = numrecs;
        self
    }

    fn dim(&mut self, name: &'a str, len: u32) -> u32 {
        let id = self.dims.len() as u32;
        self.dims.push((name, len));
        id
    }

    fn global_attr(&mut self, name: &'a str, val: AttrSpec<'a>) {
        self.global_attrs.push((name, val));
    }

    fn add_var(
        &mut self,
        name: &'a str,
        dim_ids: Vec<u32>,
        attrs: Vec<(&'a str, AttrSpec<'a>)>,
        nc_type: NcType,
        raw_bytes: Vec<u8>,
    ) {
        self.vars.push(VarSpec {
            name,
            dim_ids,
            attrs,
            nc_type,
            slices: vec![raw_bytes],
            is_record: false,
        });
    }

    fn add_record_var(
        &mut self,
        name: &'a str,
        dim_ids: Vec<u32>,
        attrs: Vec<(&'a str, AttrSpec<'a>)>,
        nc_type: NcType,
        record_slices: Vec<Vec<u8>>,
    ) {
        self.vars.push(VarSpec {
            name,
            dim_ids,
            attrs,
            nc_type,
            slices: record_slices,
            is_record: true,
        });
    }

    fn build(&self) -> Vec<u8> {
        let mut header = Vec::new();
        match self.version {
            NetCdfVersion::Classic => header.extend_from_slice(b"CDF\x01"),
            NetCdfVersion::Offset64 => header.extend_from_slice(b"CDF\x02"),
        }
        header.extend_from_slice(&self.numrecs.to_be_bytes());

        if self.dims.is_empty() {
            header.extend_from_slice(&TAG_ABSENT.to_be_bytes());
            header.extend_from_slice(&0u32.to_be_bytes());
        } else {
            header.extend_from_slice(&TAG_DIMENSION.to_be_bytes());
            header.extend_from_slice(&(self.dims.len() as u32).to_be_bytes());
            for &(name, len) in &self.dims {
                write_name(&mut header, name);
                header.extend_from_slice(&len.to_be_bytes());
            }
        }

        write_attrs(&mut header, &self.global_attrs);

        let mut begin_patch_positions = Vec::with_capacity(self.vars.len());
        if self.vars.is_empty() {
            header.extend_from_slice(&TAG_ABSENT.to_be_bytes());
            header.extend_from_slice(&0u32.to_be_bytes());
        } else {
            header.extend_from_slice(&TAG_VARIABLE.to_be_bytes());
            header.extend_from_slice(&(self.vars.len() as u32).to_be_bytes());
            for var in &self.vars {
                write_name(&mut header, var.name);
                header.extend_from_slice(&(var.dim_ids.len() as u32).to_be_bytes());
                for &dim_id in &var.dim_ids {
                    header.extend_from_slice(&dim_id.to_be_bytes());
                }
                write_attrs(&mut header, &var.attrs);
                header.extend_from_slice(&nc_type_code(var.nc_type).to_be_bytes());
                let raw_vsize = var.slices.first().map_or(0, Vec::len);
                let vsize = pad4(raw_vsize) as u32;
                header.extend_from_slice(&vsize.to_be_bytes());
                begin_patch_positions.push(header.len());
                match self.version {
                    NetCdfVersion::Classic => header.extend_from_slice(&0u32.to_be_bytes()),
                    NetCdfVersion::Offset64 => header.extend_from_slice(&0u64.to_be_bytes()),
                }
            }
        }

        let mut out = header;
        let mut begins = vec![0u64; self.vars.len()];

        for (idx, var) in self.vars.iter().enumerate() {
            if var.is_record {
                continue;
            }
            begins[idx] = out.len() as u64;
            let bytes = var.slices.first().map_or(&[][..], Vec::as_slice);
            out.extend_from_slice(bytes);
            pad_to_4(&mut out);
        }

        let rec_vars: Vec<usize> = self
            .vars
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.is_record.then_some(i))
            .collect();
        if !rec_vars.is_empty() {
            let record_count = rec_vars
                .iter()
                .map(|&i| self.vars[i].slices.len())
                .max()
                .unwrap_or(0);
            let mut rec_offset = out.len() as u64;
            for &idx in &rec_vars {
                begins[idx] = rec_offset;
                let slice_len = self.vars[idx].slices.first().map_or(0, Vec::len);
                rec_offset += pad4(slice_len) as u64;
            }
            for rec in 0..record_count {
                for &idx in &rec_vars {
                    let bytes = self.vars[idx]
                        .slices
                        .get(rec)
                        .map_or(&[][..], Vec::as_slice);
                    out.extend_from_slice(bytes);
                    pad_to_4(&mut out);
                }
            }
        }

        for (idx, patch_pos) in begin_patch_positions.into_iter().enumerate() {
            match self.version {
                NetCdfVersion::Classic => {
                    let off = (begins[idx] as u32).to_be_bytes();
                    out[patch_pos..patch_pos + 4].copy_from_slice(&off);
                }
                NetCdfVersion::Offset64 => {
                    let off = begins[idx].to_be_bytes();
                    out[patch_pos..patch_pos + 8].copy_from_slice(&off);
                }
            }
        }

        out
    }
}

fn nc_type_code(ty: NcType) -> u32 {
    match ty {
        NcType::Byte => 1,
        NcType::Char => 2,
        NcType::Short => 3,
        NcType::Int => 4,
        NcType::Float => 5,
        NcType::Double => 6,
    }
}

fn pad4(len: usize) -> usize {
    let rem = len % 4;
    if rem == 0 { len } else { len + (4 - rem) }
}

fn pad_to_4(buf: &mut Vec<u8>) {
    let rem = buf.len() % 4;
    if rem != 0 {
        buf.resize(buf.len() + (4 - rem), 0);
    }
}

fn write_name(buf: &mut Vec<u8>, name: &str) {
    let bytes = name.as_bytes();
    buf.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
    buf.extend_from_slice(bytes);
    pad_to_4(buf);
}

fn write_attrs(buf: &mut Vec<u8>, attrs: &[(&str, AttrSpec<'_>)]) {
    if attrs.is_empty() {
        buf.extend_from_slice(&TAG_ABSENT.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        return;
    }
    buf.extend_from_slice(&TAG_ATTRIBUTE.to_be_bytes());
    buf.extend_from_slice(&(attrs.len() as u32).to_be_bytes());
    for (name, val) in attrs {
        write_name(buf, name);
        match val {
            AttrSpec::Bytes(vals) => {
                buf.extend_from_slice(&1u32.to_be_bytes());
                buf.extend_from_slice(&(vals.len() as u32).to_be_bytes());
                for &v in *vals {
                    buf.push(v.to_be_bytes()[0]);
                }
                pad_to_4(buf);
            }
            AttrSpec::Chars(text) => {
                buf.extend_from_slice(&2u32.to_be_bytes());
                let bytes = text.as_bytes();
                buf.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
                buf.extend_from_slice(bytes);
                pad_to_4(buf);
            }
            AttrSpec::Shorts(vals) => {
                buf.extend_from_slice(&3u32.to_be_bytes());
                buf.extend_from_slice(&(vals.len() as u32).to_be_bytes());
                for &v in *vals {
                    buf.extend_from_slice(&v.to_be_bytes());
                }
                pad_to_4(buf);
            }
            AttrSpec::Ints(vals) => {
                buf.extend_from_slice(&4u32.to_be_bytes());
                buf.extend_from_slice(&(vals.len() as u32).to_be_bytes());
                for &v in *vals {
                    buf.extend_from_slice(&v.to_be_bytes());
                }
                pad_to_4(buf);
            }
            AttrSpec::Floats(vals) => {
                buf.extend_from_slice(&5u32.to_be_bytes());
                buf.extend_from_slice(&(vals.len() as u32).to_be_bytes());
                for &v in *vals {
                    buf.extend_from_slice(&v.to_be_bytes());
                }
                pad_to_4(buf);
            }
            AttrSpec::Doubles(vals) => {
                buf.extend_from_slice(&6u32.to_be_bytes());
                buf.extend_from_slice(&(vals.len() as u32).to_be_bytes());
                for &v in *vals {
                    buf.extend_from_slice(&v.to_be_bytes());
                }
                pad_to_4(buf);
            }
        }
    }
}

fn f32_bytes(vals: &[f32]) -> Vec<u8> {
    vals.iter().flat_map(|v| v.to_be_bytes()).collect()
}

fn f64_bytes(vals: &[f64]) -> Vec<u8> {
    vals.iter().flat_map(|v| v.to_be_bytes()).collect()
}

fn i16_bytes(vals: &[i16]) -> Vec<u8> {
    vals.iter().flat_map(|v| v.to_be_bytes()).collect()
}

fn fixed_strings(names: &[&str], width: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(names.len() * width);
    for &name in names {
        let mut slot = vec![0u8; width];
        let bytes = name.as_bytes();
        let copy_len = bytes.len().min(width);
        slot[..copy_len].copy_from_slice(&bytes[..copy_len]);
        out.extend_from_slice(&slot);
    }
    out
}

fn build_sample_cdf1() -> Vec<u8> {
    let mut b = NetCdfBuilder::new(NetCdfVersion::Classic);
    let point_number = b.dim("point_number", 5);
    let peak_number = b.dim("peak_number", 2);
    let str32 = b.dim("_32_byte_string", 32);
    let unit_len = b.dim("_16_byte_string", 16);

    b.global_attr("sample_name", AttrSpec::Chars("RXN-042-crude"));
    b.global_attr(
        "dataset_date_time_stamp",
        AttrSpec::Chars("20260415143000+0800"),
    );
    b.global_attr("detector_name", AttrSpec::Chars("Agilent 1260 DAD 254nm"));
    b.global_attr("experiment_title", AttrSpec::Chars("C18_10min_grad"));
    b.global_attr("source_file_reference", AttrSpec::Chars("042_crude.D"));

    b.add_var(
        "actual_delay_time",
        vec![],
        vec![],
        NcType::Float,
        f32_bytes(&[30.0]),
    );
    b.add_var(
        "actual_run_time_length",
        vec![],
        vec![],
        NcType::Float,
        f32_bytes(&[120.0]),
    );
    b.add_var(
        "actual_sampling_interval",
        vec![],
        vec![],
        NcType::Float,
        f32_bytes(&[30.0]),
    );
    b.add_var(
        "detector_unit",
        vec![unit_len],
        vec![],
        NcType::Char,
        fixed_strings(&["mAU"], 16),
    );
    b.add_var(
        "ordinate_values",
        vec![point_number],
        vec![],
        NcType::Float,
        f32_bytes(&[1.5, 12.0, 240.5, 18.25, 2.0]),
    );

    b.add_var(
        "peak_retention_time",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[90.0, 150.0]),
    );
    b.add_var(
        "peak_area",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[12_500.0, 37_500.0]),
    );
    b.add_var(
        "peak_area_percent",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[25.0, 75.0]),
    );
    b.add_var(
        "peak_height",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[80.0, 240.5]),
    );
    b.add_var(
        "peak_width",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[6.0, 9.0]),
    );
    b.add_var(
        "peak_start_time",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[84.0, 141.0]),
    );
    b.add_var(
        "peak_end_time",
        vec![peak_number],
        vec![],
        NcType::Float,
        f32_bytes(&[96.0, 159.0]),
    );
    b.add_var(
        "peak_name",
        vec![peak_number, str32],
        vec![],
        NcType::Char,
        fixed_strings(&["Starting Material", "Product A"], 32),
    );

    b.build()
}

#[test]
fn parses_cdf1_chromatogram_peaks_and_metadata_via_detect_and_read() {
    let bytes = build_sample_cdf1();
    assert_eq!(detect(Some("run.cdf"), &bytes), VendorFormat::AndiCdf);
    assert_eq!(detect(None, &bytes), VendorFormat::AndiCdf);

    let dataset = parse_cdf(&bytes).unwrap();
    assert_eq!(dataset.version, NetCdfVersion::Classic);
    assert_eq!(dataset.sample_name.as_deref(), Some("RXN-042-crude"));
    assert_eq!(
        dataset.acquired_at.as_deref(),
        Some("2026-04-15T14:30:00+08:00")
    );
    assert_eq!(
        dataset.detector_name.as_deref(),
        Some("Agilent 1260 DAD 254nm")
    );
    assert_eq!(dataset.experiment_title.as_deref(), Some("C18_10min_grad"));
    assert_eq!(
        dataset.source_file_reference.as_deref(),
        Some("042_crude.D")
    );
    assert_eq!(dataset.detector_unit.as_deref(), Some("mAU"));
    assert_eq!(dataset.delay_time_s, Some(30.0));
    assert_eq!(dataset.run_time_length_s, Some(120.0));
    assert_eq!(dataset.sampling_interval_s, Some(30.0));

    assert_eq!(dataset.traces.len(), 1);
    let trace = &dataset.traces[0];
    assert_eq!(trace.label, "Agilent 1260 DAD 254nm");
    assert_eq!(trace.x_unit, "min");
    assert_eq!(trace.y_unit, "mAU");
    assert_eq!(
        trace.points,
        vec![
            (0.5, 1.5),
            (1.0, 12.0),
            (1.5, 240.5),
            (2.0, 18.25),
            (2.5, 2.0),
        ]
    );

    assert_eq!(dataset.peaks.len(), 2);
    assert_eq!(dataset.peaks[0].rt_min, Some(1.5));
    assert_eq!(dataset.peaks[0].area, Some(12_500.0));
    assert_eq!(dataset.peaks[0].area_pct, Some(25.0));
    assert_eq!(dataset.peaks[0].height, Some(80.0));
    assert_eq!(dataset.peaks[0].width_min, Some(0.1));
    assert_eq!(dataset.peaks[0].start_time_min, Some(1.4));
    assert_eq!(dataset.peaks[0].end_time_min, Some(1.6));
    assert_eq!(
        dataset.peaks[0].compound.as_deref(),
        Some("Starting Material")
    );

    assert_eq!(dataset.peaks[1].rt_min, Some(2.5));
    assert_eq!(dataset.peaks[1].compound.as_deref(), Some("Product A"));

    let reading = read(Some("run.cdf"), &bytes).unwrap();
    assert_eq!(reading.technique, Technique::Hplc);
    assert_eq!(reading.vendor_format, VendorFormat::AndiCdf);
    assert_eq!(reading.sample_name.as_deref(), Some("RXN-042-crude"));
    assert_eq!(
        reading.instrument.as_deref(),
        Some("Agilent 1260 DAD 254nm")
    );
    assert_eq!(reading.method.name.as_deref(), Some("C18_10min_grad"));
    assert_eq!(reading.peaks.len(), 2);
    assert_eq!(
        reading.main_peak().and_then(|p| p.compound.as_deref()),
        Some("Product A")
    );

    let bundle_reading = read_files(&[
        ("bundle/notes.txt", b"operator note"),
        ("bundle/run.cdf", bytes.as_slice()),
    ])
    .unwrap();
    assert_eq!(bundle_reading, reading);
}

#[test]
fn parses_cdf2_record_variables_all_types_scaling_and_streaming_numrecs() {
    let mut b = NetCdfBuilder::new(NetCdfVersion::Offset64).numrecs(0xFFFF_FFFF);
    let rec_dim = b.dim("point_number", 0);
    let peak_dim = b.dim("peak_number", 2);

    b.global_attr("sample_name", AttrSpec::Chars("CDF2-Sample"));
    b.global_attr(
        "dataset_date_time_stamp",
        AttrSpec::Chars("2026-05-01T08:15:00Z"),
    );
    b.global_attr("detector_name", AttrSpec::Chars("FID-A"));
    b.global_attr("byte_flags", AttrSpec::Bytes(&[-3, 0, 7]));
    b.global_attr("short_ids", AttrSpec::Shorts(&[-12, 345, 9]));
    b.global_attr("int_counts", AttrSpec::Ints(&[100_000, -42]));
    b.global_attr("float_gain", AttrSpec::Floats(&[1.25]));
    b.global_attr("double_ref", AttrSpec::Doubles(&[12.345_678_9]));

    b.add_var(
        "actual_delay_time",
        vec![],
        vec![("units", AttrSpec::Chars("minutes"))],
        NcType::Double,
        f64_bytes(&[0.25]),
    );
    b.add_var(
        "actual_run_time_length",
        vec![],
        vec![("units", AttrSpec::Chars("minutes"))],
        NcType::Double,
        f64_bytes(&[0.75]),
    );

    b.add_var(
        "peak_retention_time",
        vec![peak_dim],
        vec![
            ("units", AttrSpec::Chars("min")),
            ("_FillValue", AttrSpec::Doubles(&[-999.0])),
        ],
        NcType::Double,
        f64_bytes(&[0.5, -999.0]),
    );
    b.add_var(
        "peak_area",
        vec![peak_dim],
        vec![],
        NcType::Double,
        f64_bytes(&[4200.0, 9.969_209_968_386_869e36]),
    );

    b.add_record_var(
        "ordinate_values",
        vec![rec_dim],
        vec![
            ("units", AttrSpec::Chars("pA")),
            ("scale_factor", AttrSpec::Doubles(&[0.5])),
            ("add_offset", AttrSpec::Doubles(&[10.0])),
        ],
        NcType::Short,
        vec![
            i16_bytes(&[2]),
            i16_bytes(&[20]),
            i16_bytes(&[100]),
            i16_bytes(&[4]),
        ],
    );
    b.add_record_var(
        "point_flag",
        vec![rec_dim],
        vec![],
        NcType::Byte,
        vec![vec![1], vec![2], vec![3], vec![4]],
    );

    let bytes = b.build();
    assert_eq!(detect(None, &bytes), VendorFormat::AndiCdf);

    let nc = NetCdfFile::parse(&bytes).unwrap();
    assert_eq!(nc.version(), NetCdfVersion::Offset64);
    assert_eq!(nc.numrecs(), 4);
    assert_eq!(
        nc.global_attrs().get("byte_flags"),
        Some(&NcAttrValue::Bytes(vec![-3, 0, 7]))
    );
    assert_eq!(
        nc.global_attrs().get("short_ids"),
        Some(&NcAttrValue::Shorts(vec![-12, 345, 9]))
    );
    assert_eq!(
        nc.global_attrs().get("int_counts"),
        Some(&NcAttrValue::Ints(vec![100_000, -42]))
    );
    assert_eq!(
        nc.global_attrs().get("float_gain"),
        Some(&NcAttrValue::Floats(vec![1.25]))
    );
    assert_eq!(
        nc.global_attrs().get("double_ref"),
        Some(&NcAttrValue::Doubles(vec![12.345_678_9]))
    );
    assert_eq!(
        nc.read_var_f64("point_flag").unwrap(),
        Some(vec![1.0, 2.0, 3.0, 4.0])
    );

    let reading = read_cdf(&bytes).unwrap();
    assert_eq!(reading.sample_name.as_deref(), Some("CDF2-Sample"));
    assert_eq!(reading.acquired_at.as_deref(), Some("2026-05-01T08:15:00Z"));
    assert_eq!(reading.instrument.as_deref(), Some("FID-A"));
    assert_eq!(reading.traces.len(), 1);
    let trace = &reading.traces[0];
    assert_eq!(trace.y_unit, "pA");
    assert_eq!(
        trace.points,
        vec![(0.25, 11.0), (0.5, 20.0), (0.75, 60.0), (1.0, 12.0)]
    );

    assert_eq!(reading.peaks.len(), 1);
    assert_eq!(reading.peaks[0].rt_min, Some(0.5));
    assert_eq!(reading.peaks[0].area, Some(4200.0));
}

#[test]
fn rejects_malformed_truncated_and_adversarial_cdf_inputs() {
    assert!(matches!(
        read_cdf(b""),
        Err(InstrumentError::InvalidAndi(_))
    ));
    assert!(matches!(
        read_cdf(b"CDF\x03\x00\x00\x00\x00"),
        Err(InstrumentError::InvalidAndi(_))
    ));

    let valid = build_sample_cdf1();
    for cut in 0..valid.len() {
        let res = read_cdf(&valid[..cut]);
        assert!(
            matches!(res, Err(InstrumentError::InvalidAndi(_))),
            "expected InvalidAndi when truncated to {cut}/{} bytes",
            valid.len()
        );
    }

    let mut bad_absent = Vec::new();
    bad_absent.extend_from_slice(b"CDF\x01");
    bad_absent.extend_from_slice(&0u32.to_be_bytes());
    bad_absent.extend_from_slice(&0u32.to_be_bytes());
    bad_absent.extend_from_slice(&3u32.to_be_bytes());
    assert!(matches!(
        read_cdf(&bad_absent),
        Err(InstrumentError::InvalidAndi(_))
    ));

    let mut huge_dims = Vec::new();
    huge_dims.extend_from_slice(b"CDF\x01");
    huge_dims.extend_from_slice(&0u32.to_be_bytes());
    huge_dims.extend_from_slice(&TAG_DIMENSION.to_be_bytes());
    huge_dims.extend_from_slice(&u32::MAX.to_be_bytes());
    assert!(matches!(
        read_cdf(&huge_dims),
        Err(InstrumentError::InvalidAndi(_))
    ));

    let mut bad_type = Vec::new();
    bad_type.extend_from_slice(b"CDF\x01");
    bad_type.extend_from_slice(&0u32.to_be_bytes());
    bad_type.extend_from_slice(&[0u8; 8]);
    bad_type.extend_from_slice(&TAG_ATTRIBUTE.to_be_bytes());
    bad_type.extend_from_slice(&1u32.to_be_bytes());
    write_name(&mut bad_type, "x");
    bad_type.extend_from_slice(&99u32.to_be_bytes());
    bad_type.extend_from_slice(&1u32.to_be_bytes());
    bad_type.extend_from_slice(&0u32.to_be_bytes());
    assert!(matches!(
        read_cdf(&bad_type),
        Err(InstrumentError::InvalidAndi(_))
    ));

    let mut bad_dim_id = Vec::new();
    bad_dim_id.extend_from_slice(b"CDF\x01");
    bad_dim_id.extend_from_slice(&0u32.to_be_bytes());
    bad_dim_id.extend_from_slice(&[0u8; 8]);
    bad_dim_id.extend_from_slice(&[0u8; 8]);
    bad_dim_id.extend_from_slice(&TAG_VARIABLE.to_be_bytes());
    bad_dim_id.extend_from_slice(&1u32.to_be_bytes());
    write_name(&mut bad_dim_id, "v");
    bad_dim_id.extend_from_slice(&1u32.to_be_bytes());
    bad_dim_id.extend_from_slice(&4u32.to_be_bytes());
    bad_dim_id.extend_from_slice(&[0u8; 8]);
    bad_dim_id.extend_from_slice(&5u32.to_be_bytes());
    bad_dim_id.extend_from_slice(&4u32.to_be_bytes());
    bad_dim_id.extend_from_slice(&0u32.to_be_bytes());
    assert!(matches!(
        read_cdf(&bad_dim_id),
        Err(InstrumentError::InvalidAndi(_))
    ));
}
