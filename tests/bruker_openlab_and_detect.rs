//! Integration tests for Bruker NMR, Agilent OpenLab, format detection, and downsampling.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use std::io::Write;

use flate2::Compression;
use flate2::write::DeflateEncoder;
use sequency_instruments::bruker::{
    ByteOrder, integral, multiplet_center, parse_acqus, parse_procs, read_1r, read_dir, read_pdata,
};
use sequency_instruments::openlab::{inspect_dx, inspect_rx, parse_acaml, parse_rx_peaks};
use sequency_instruments::{
    InstrumentError, Technique, Trace, VendorFormat, detect, lttb, minmax_bins, read, read_files,
};

/// Builds an in-memory ZIP archive from `(name, data, deflate)` entries.
fn build_zip(entries: &[(&str, &[u8], bool)]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut central = Vec::new();

    for &(name, raw, deflate) in entries {
        let lfh_offset = out.len() as u32;
        let crc = crc32fast::hash(raw);
        let (method, payload): (u16, Vec<u8>) = if deflate {
            let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
            enc.write_all(raw).unwrap();
            (8, enc.finish().unwrap())
        } else {
            (0, raw.to_vec())
        };

        let name_bytes = name.as_bytes();
        out.extend_from_slice(&0x0403_4B50u32.to_le_bytes());
        out.extend_from_slice(&20u16.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(&method.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(&crc.to_le_bytes());
        out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        out.extend_from_slice(&(raw.len() as u32).to_le_bytes());
        out.extend_from_slice(&(name_bytes.len() as u16).to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        out.extend_from_slice(name_bytes);
        out.extend_from_slice(&payload);

        central.extend_from_slice(&0x0201_4B50u32.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes());
        central.extend_from_slice(&20u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&method.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&crc.to_le_bytes());
        central.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        central.extend_from_slice(&(raw.len() as u32).to_le_bytes());
        central.extend_from_slice(&(name_bytes.len() as u16).to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u16.to_le_bytes());
        central.extend_from_slice(&0u32.to_le_bytes());
        central.extend_from_slice(&lfh_offset.to_le_bytes());
        central.extend_from_slice(name_bytes);
    }

    let cd_offset = out.len() as u32;
    let cd_size = central.len() as u32;
    out.extend_from_slice(&central);

    out.extend_from_slice(&0x0605_4B50u32.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out.extend_from_slice(&(entries.len() as u16).to_le_bytes());
    out.extend_from_slice(&(entries.len() as u16).to_le_bytes());
    out.extend_from_slice(&cd_size.to_le_bytes());
    out.extend_from_slice(&cd_offset.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out
}

fn build_acaml(rows: &[(&str, &str, &str, &str, &str)], login: &str) -> Vec<u8> {
    let mut meta = String::new();
    let mut measured = String::new();
    for (idx, &(data, sample, kind, at, status)) in rows.iter().enumerate() {
        let n = idx + 1;
        meta.push_str(&format!(
            "<InjectionMetaData AcqMethodName=\"GCCD_HP-5_300\" DaMethodName=\"面积百分比\" \
             InjectionAcqDateTime=\"{at}\" InjectionId=\"i{n}\" InstrumentName=\"Agilent 8890 GC\" \
             InjectionVolume=\"1.5\" RawDataFileName=\"{data}.dx\" SampleName=\"{sample}\" \
             SampleType=\"{kind}\" VialNumber=\"{n}\" />"
        ));
        measured.push_str(&format!(
            "<MeasData id=\"i{n}\"><Info><CreatedBy><Username>{login}</Username></CreatedBy></Info>\
             <AcquitionStatus><Status>{status}</Status></AcquitionStatus>\
             <ExternalResultPath>{data}.rx</ExternalResultPath></MeasData>"
        ));
    }

    format!(
        "<ACAML xmlns=\"urn:schemas-agilent-com:acaml21\"><Doc><DocInfo>\
         <CreatedByUser>{login}</CreatedByUser>\
         <CustomField Name=\"InjectionMetaDataItems\"><Xml>\
         <ArrayOfInjectionMetaData xmlns=\"\">{meta}</ArrayOfInjectionMetaData>\
         </Xml></CustomField></DocInfo>\
         <Content><Injections>{measured}</Injections></Content></Doc></ACAML>"
    )
    .into_bytes()
}

fn build_rx(peaks: &[(f64, f64, f64, f64, &str)], login: &str) -> Vec<u8> {
    let mut found = String::new();
    let mut compounds = String::new();
    for (n, &(rt, area, pct, height, name)) in peaks.iter().enumerate() {
        found.push_str(&format!(
            "<Peak id=\"p{n}\"><RetentionTime val=\"{rt}\" unit=\"min\" />\
             <Area val=\"{area}\" unit=\"pA·s\" /><AreaPercent val=\"{pct}\" />\
             <Height val=\"{height}\" unit=\"pA\" /></Peak>"
        ));
        compounds.push_str(&format!(
            "<InjectionCompound><CompoundName>{name}</CompoundName><Identification><Qualified>\
             <Peaks><Peak_ID id=\"p{n}\" /></Peaks></Qualified></Identification></InjectionCompound>"
        ));
    }
    let xml = format!(
        "<ACAML xmlns=\"urn:schemas-agilent-com:acaml21\"><Doc><DocInfo>\
         <CreatedByUser>{login}</CreatedByUser></DocInfo><Content><Injections><Result>\
         <SignalResult>{found}</SignalResult>{compounds}</Result></Injections></Content></Doc></ACAML>"
    );
    build_zip(&[
        ("Base/InjectionACAML", xml.as_bytes(), true),
        ("Base/AuditTrail", b"\x0a\x06SYSTEM", false),
    ])
}

fn write_ch_utf16(buf: &mut [u8], offset: usize, text: &str) {
    let chars: Vec<u16> = text.encode_utf16().collect();
    buf[offset] = chars.len() as u8;
    for (i, &u) in chars.iter().enumerate() {
        let at = offset + 1 + i * 2;
        buf[at..at + 2].copy_from_slice(&u.to_le_bytes());
    }
}

fn build_dx(sample: &str, values: &[f64], scale: f64, signal_desc: &str, unit: &str) -> Vec<u8> {
    let mut ch = vec![0u8; 0x1800];
    ch[0..4].copy_from_slice(b"\x03179");
    ch[0x11A..0x11E].copy_from_slice(&0.0_f32.to_be_bytes());
    ch[0x11E..0x122].copy_from_slice(&60_000.0_f32.to_be_bytes());
    ch[0x127C..0x1284].copy_from_slice(&scale.to_be_bytes());
    write_ch_utf16(&mut ch, 0x35A, sample);
    write_ch_utf16(&mut ch, 0xA0E, "D:\\CDS\\TBP\\HPLC_DAD_254.amx");
    write_ch_utf16(&mut ch, 0xC11, "Agilent 1260 LC");
    write_ch_utf16(&mut ch, 0x104C, unit);
    write_ch_utf16(&mut ch, 0x1075, signal_desc);

    for &v in values {
        ch.extend_from_slice(&v.to_le_bytes());
    }

    let acmd = format!(
        "<ACMD><InjectionInfo><SampleName>{sample}</SampleName>\
         <AcquisitionMethod>D:\\CDS\\TBP\\HPLC_DAD_254.amx</AcquisitionMethod></InjectionInfo></ACMD>"
    );

    build_zip(&[
        ("12af0742.CH", &ch, true),
        ("injection.acmd", acmd.as_bytes(), true),
    ])
}

#[test]
fn bruker_acqus_and_1r_respect_bytordp_endianness_and_nc_proc_scaling() {
    let acqus_text = concat!(
        "##TITLE= Parameter file, TopSpin 4.5.0\n",
        "##$INSTRUM= <Avance NEO 600>\n",
        "##$NUC1= <1H>\n",
        "##$PULPROG= <noesyigld1d>\n",
        "##$SOLVENT= <CDCl3>\n",
        "##$NS= 16\n",
        "##$D= (0..3)\n",
        "0 10.0 0.002 0.0\n",
        "$$ 2026-05-29 14:52:03.540 +0800  root@CZC2439XZZ\n",
        "$$ D:/nmr_data/DB/20260529/operator/20260529-DB-TBP-260474-D1/1/acqus\n",
        "##END=\n"
    );

    let acqus = parse_acqus(acqus_text).unwrap();
    assert_eq!(
        acqus.acquired_at.as_deref(),
        Some("2026-05-29 14:52:03 +0800")
    );
    assert_eq!(acqus.nucleus, "1H");
    assert_eq!(acqus.pulse_program, "noesyigld1d");
    assert_eq!(acqus.instrument, "Avance NEO 600");
    assert_eq!(acqus.solvent, "CDCl3");
    assert_eq!(acqus.scans, 16);
    assert_eq!(acqus.relaxation_delay_s, 10.0);
    assert_eq!(
        acqus.data_path.as_deref(),
        Some("D:/nmr_data/DB/20260529/operator/20260529-DB-TBP-260474-D1/1")
    );

    let raw_points = [100i32, -40, 2000, 512];
    let le_1r: Vec<u8> = raw_points.iter().flat_map(|v| v.to_le_bytes()).collect();
    let be_1r: Vec<u8> = raw_points.iter().flat_map(|v| v.to_be_bytes()).collect();

    let procs_le = parse_procs(
        "##$BYTORDP= 0\n##$DTYPP= 0\n##$NC_proc= 3\n##$OFFSET= 12.0\n##$SF= 600.0\n##$SI= 4\n##$SW_p= 7200.0\n##END=\n",
    )
    .unwrap();
    let procs_be = parse_procs(
        "##$BYTORDP= 1\n##$DTYPP= 0\n##$NC_proc= -2\n##$OFFSET= 12.0\n##$SF= 600.0\n##$SI= 4\n##$SW_p= 7200.0\n##END=\n",
    )
    .unwrap();
    assert_eq!(procs_le.byte_order, ByteOrder::LittleEndian);
    assert_eq!(procs_be.byte_order, ByteOrder::BigEndian);

    let trace_le = read_1r(&le_1r, &procs_le, Some(&acqus)).unwrap();
    let trace_be = read_1r(&be_1r, &procs_be, Some(&acqus)).unwrap();

    assert_eq!(trace_le.label, "1H NMR");
    assert_eq!(trace_le.x_unit, "ppm");
    assert_eq!(
        trace_le.points,
        vec![(12.0, 800.0), (9.0, -320.0), (6.0, 16_000.0), (3.0, 4096.0)]
    );
    assert_eq!(
        trace_be.points,
        vec![(12.0, 25.0), (9.0, -10.0), (6.0, 500.0), (3.0, 128.0)]
    );
}

#[test]
fn bruker_qnmr_synthetic_spectrum_integrates_and_reads_from_directory() {
    let size = 4096usize;
    let mhz = 600.13_f64;
    let offset = 12.0_f64;
    let width_ppm = 12.0_f64;
    let sweep_hz = width_ppm * mhz;
    let sfo1 = mhz + (offset * mhz - sweep_hz / 2.0) / 1.0e6;

    let mut data_1r = Vec::with_capacity(size * 4);
    for i in 0..size {
        let ppm = offset - (i as f64) * width_ppm / (size as f64);
        let std_peak = 300_000.0 / (1.0 + ((ppm - 6.08) / 0.004).powi(2));
        let smp_peak = 100_000.0 / (1.0 + ((ppm - 7.50) / 0.004).powi(2));
        let y = (std_peak + smp_peak).round() as i32;
        data_1r.extend_from_slice(&y.to_le_bytes());
    }

    let acqus = format!(
        "##TITLE= qNMR\n##$INSTRUM= <Avance NEO 600>\n##$NUC1= <1H>\n##$PULPROG= <zg30>\n\
         ##$SOLVENT= <DMSO>\n##$NS= 32\n##$SW_h= {sweep_hz}\n##$SFO1= {sfo1}\n\
         $$ 2026-08-07 10:00:00.000 +0800  user@host\n"
    );
    let procs = format!(
        "##$BYTORDP= 0\n##$DTYPP= 0\n##$NC_proc= 0\n##$OFFSET= {offset}\n\
         ##$SF= {mhz}\n##$SI= {size}\n##$SW_p= {sweep_hz}\n##END=\n"
    );

    let reading = read_pdata(Some(acqus.as_bytes()), procs.as_bytes(), &data_1r).unwrap();
    assert_eq!(reading.technique, Technique::Nmr);
    assert_eq!(reading.vendor_format, VendorFormat::BrukerNmr);
    assert_eq!(reading.method.nucleus.as_deref(), Some("1H"));
    assert_eq!(reading.method.frequency_mhz, Some(mhz));
    assert_eq!(reading.method.solvents, vec!["DMSO"]);
    assert_eq!(reading.method.scans, Some(32));

    let pts = &reading.traces[0].points;
    assert!((pts[0].0 - 12.0).abs() < 1e-6);
    assert!(pts[0].0 > pts[size - 1].0);

    let std_area = integral(pts, 6.02, 6.14, 0.02);
    let smp_area = integral(pts, 7.44, 7.56, 0.02);
    assert!((std_area / smp_area - 3.0).abs() < 0.05);

    let center = multiplet_center(pts, 7.45, 7.55, 0.025).unwrap();
    assert!((center - 7.50).abs() < 0.005);

    let tmp = std::env::temp_dir().join(format!("sequency_bruker_test_{}", std::process::id()));
    let pdata_dir = tmp.join("1").join("pdata").join("1");
    std::fs::create_dir_all(&pdata_dir).unwrap();
    std::fs::write(tmp.join("1").join("acqus"), acqus.as_bytes()).unwrap();
    std::fs::write(pdata_dir.join("procs"), procs.as_bytes()).unwrap();
    std::fs::write(pdata_dir.join("1r"), &data_1r).unwrap();
    std::fs::write(pdata_dir.join("title"), b"TBP-260474-Q\nsecond line\n").unwrap();

    let dir_reading = read_dir(&tmp.join("1")).unwrap();
    assert_eq!(dir_reading.sample_name.as_deref(), Some("TBP-260474-Q"));
    assert_eq!(dir_reading.traces[0].points.len(), size);
    let _ = std::fs::remove_dir_all(&tmp);
}

#[test]
fn openlab_acaml_rx_and_dx_extract_metadata_peaks_and_signals() {
    let acaml = build_acaml(
        &[
            ("B1", "Blank", "Blank", "1601-01-01T00:00:00Z", "Aborted"),
            (
                "S1",
                "  TBP-250401 ",
                "Sample",
                "2025-04-01T04:21:26.904Z",
                "Completed",
            ),
        ],
        "chemist1",
    );

    let parsed = parse_acaml(&acaml).unwrap();
    assert_eq!(parsed.injections.len(), 2);
    assert!(parsed.logins.contains("chemist1"));
    assert_eq!(parsed.injections[0].acquired_utc, None);
    assert_eq!(parsed.injections[0].status.as_deref(), Some("Aborted"));
    assert_eq!(parsed.injections[1].data_file, "S1");
    assert_eq!(parsed.injections[1].sample_name, "TBP-250401");
    assert_eq!(
        parsed.injections[1].acquired_utc.as_deref(),
        Some("2025-04-01T04:21:26.904Z")
    );
    assert_eq!(parsed.injections[1].status.as_deref(), Some("Completed"));
    assert_eq!(parsed.injections[1].result_file.as_deref(), Some("S1.rx"));
    assert_eq!(parsed.injections[1].injection_volume_ul, Some(1.5));

    let acaml_str = std::str::from_utf8(&acaml).unwrap();
    let mut acaml_utf16 = vec![0xFF, 0xFE];
    acaml_utf16.extend(acaml_str.encode_utf16().flat_map(u16::to_le_bytes));
    assert_eq!(parse_acaml(&acaml_utf16).unwrap(), parsed);

    let rx = build_rx(
        &[
            (5.1204, 796.836, 61.919, 42.1256, "TBP-250401"),
            (2.38078, 203.164, 38.081, 1.5, "乙醇"),
        ],
        "SYSTEM",
    );
    let rx_info = inspect_rx(&rx).unwrap();
    assert_eq!(rx_info.entries.len(), 2);
    assert!(rx_info.entries[0].is_compressed);
    assert!(!rx_info.entries[1].is_compressed);
    let peaks = parse_rx_peaks(&rx).unwrap();
    assert_eq!(peaks.len(), 2);
    assert_eq!(
        (
            peaks[0].rt_min,
            peaks[0].area,
            peaks[0].area_pct,
            peaks[0].compound.as_deref()
        ),
        (Some(2.381), Some(203.164), Some(38.08), Some("乙醇"))
    );
    assert_eq!(
        (peaks[1].rt_min, peaks[1].compound.as_deref()),
        (Some(5.12), Some("TBP-250401"))
    );

    let dx = build_dx(
        "TBP-250401",
        &[10.0, 25.0, 5.0],
        2.0,
        "DAD1A,Sig=254.0,4.0 Ref=off",
        "mAU",
    );
    let dx_info = inspect_dx(&dx).unwrap();
    assert_eq!(dx_info.sample_name.as_deref(), Some("TBP-250401"));
    assert_eq!(dx_info.method_name.as_deref(), Some("HPLC_DAD_254.amx"));
    assert_eq!(dx_info.wavelength_nm, Some(254.0));
    assert_eq!(dx_info.traces.len(), 1);
    assert_eq!(dx_info.traces[0].y_unit, "mAU");
    assert_eq!(
        dx_info.traces[0].points,
        vec![(0.0, 20.0), (0.5, 50.0), (1.0, 10.0)]
    );

    let snapshot_dx = build_dx("SNAP", &[1.0, 2.0], 1.0, "FID1A", "pA");
    let combined = read_files(&[
        ("R.rslt/R.acaml", &acaml),
        ("R.rslt/S1.rx", &rx),
        ("R.rslt/S1.dx", &dx),
        ("R.rslt/snapshot-01.dx", &snapshot_dx),
    ])
    .unwrap();
    assert_eq!(combined.vendor_format, VendorFormat::OpenLabAcaml);
    assert_eq!(combined.sample_name.as_deref(), Some("TBP-250401"));
    assert_eq!(combined.method.wavelength_nm, Some(254.0));
    assert_eq!(combined.peaks.len(), 2);
    assert_eq!(combined.traces.len(), 1);
    assert_eq!(combined.traces[0].points[1], (0.5, 50.0));

    let mut corrupt_rx = rx.clone();
    let mid = corrupt_rx.len() / 2;
    corrupt_rx[mid] ^= 0xFF;
    assert!(matches!(
        inspect_rx(&corrupt_rx),
        Err(InstrumentError::InvalidZip(_))
    ));
}

#[test]
fn detect_identifies_formats_by_magic_bytes_and_filenames() {
    let acaml = build_acaml(
        &[("S1", "TBP-1", "Sample", "2026-01-01T00:00:00Z", "Completed")],
        "SYSTEM",
    );
    let rx = build_rx(&[(1.0, 100.0, 100.0, 10.0, "A")], "SYSTEM");
    let dx = build_dx("TBP-1", &[1.0, 2.0], 1.0, "FID1A", "pA");
    let acqus = b"##TITLE= Parameter file, TopSpin 4.5.0\n##$NUC1= <1H>\n";

    assert_eq!(detect(None, &acaml), VendorFormat::OpenLabAcaml);
    assert_eq!(detect(None, &rx), VendorFormat::OpenLabRx);
    assert_eq!(detect(None, &dx), VendorFormat::OpenLabDx);
    assert_eq!(detect(None, acqus), VendorFormat::BrukerNmr);
    assert_eq!(
        detect(Some("exp/1/pdata/1/1r"), &[0, 0, 0, 1]),
        VendorFormat::BrukerNmr
    );

    let mut finnigan = vec![0x01, 0xA1];
    for ch in "Finnigan".encode_utf16() {
        finnigan.extend_from_slice(&ch.to_le_bytes());
    }
    assert_eq!(detect(None, &finnigan), VendorFormat::ThermoRaw);
    assert_eq!(
        detect(None, b"Mestrelab Research S.L.\xf1\xe2\xd3\xc4"),
        VendorFormat::MestReNova
    );
    assert_eq!(detect(None, b"SignHdr\0extra"), VendorFormat::ChromeleonRaw);
    assert_eq!(detect(None, b"CDF\x01\x00\x00"), VendorFormat::AndiCdf);
    assert_eq!(detect(Some("TBP-1.spa"), b""), VendorFormat::OmnicSpa);
    assert_eq!(detect(Some("TBP-1.all"), b"200 0.1"), VendorFormat::UvText);
    assert_eq!(detect(None, b"random bytes"), VendorFormat::Unknown);

    let reading = read(Some("R.acaml"), &acaml).unwrap();
    let json = serde_json::to_string(&reading).unwrap();
    let roundtrip: sequency_instruments::InstrumentReading = serde_json::from_str(&json).unwrap();
    assert_eq!(roundtrip, reading);
}

#[test]
fn downsampling_preserves_endpoints_and_sharp_peaks() {
    let mut points = Vec::with_capacity(1000);
    for i in 0..1000 {
        let x = i as f64 * 0.01;
        let y = match i {
            350 => 950.0,
            351 => -120.0,
            720 => 420.0,
            _ => (x * 0.5).sin(),
        };
        points.push((x, y));
    }

    let trace = Trace {
        label: "Test".to_owned(),
        x_unit: "min".to_owned(),
        y_unit: "mAU".to_owned(),
        points: points.clone(),
    };

    let lttb_trace = trace.downsample_lttb(64);
    assert_eq!(lttb_trace.points.len(), 64);
    assert_eq!(lttb_trace.points.first(), points.first());
    assert_eq!(lttb_trace.points.last(), points.last());
    let lttb_max = lttb_trace
        .points
        .iter()
        .map(|p| p.1)
        .fold(f64::NEG_INFINITY, f64::max);
    assert_eq!(lttb_max, 950.0);

    let env_trace = trace.downsample_minmax(50);
    assert!(env_trace.points.len() <= 100);
    let env_max = env_trace
        .points
        .iter()
        .map(|p| p.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let env_min = env_trace
        .points
        .iter()
        .map(|p| p.1)
        .fold(f64::INFINITY, f64::min);
    assert_eq!(env_max, 950.0);
    assert_eq!(env_min, -120.0);
    for pair in env_trace.points.windows(2) {
        assert!(pair[0].0 <= pair[1].0);
    }

    assert!(lttb(&points, 0).is_empty());
    assert_eq!(lttb(&points, 1), vec![points[0]]);
    assert_eq!(lttb(&points, 2), vec![points[0], points[999]]);
    assert_eq!(lttb(&points[..5], 10).len(), 5);
    assert!(minmax_bins(&points, 0).is_empty());
    assert_eq!(minmax_bins(&points[..5], 10).len(), 5);
}

#[test]
fn bruker_19f_negative_offset_float64_and_validation_errors() {
    let acqus = b"##TITLE= 19F\n##$NUC1= <19F>\n##$PULPROG= <zgfhigqn.2>\n##$SOLVENT= <CDCl3>\n";
    let procs = b"##$BYTORDP= 1\n##$DTYPP= 2\n##$NC_proc= -1\n##$OFFSET= -50.0\n##$SF= 564.6\n##$SI= 3\n##$SW_p= 56460.0\n";
    let raw_f64 = [100.0_f64, -20.0_f64, 50.0_f64];
    let data_1r: Vec<u8> = raw_f64.iter().flat_map(|v| v.to_be_bytes()).collect();

    let reading = read_pdata(Some(acqus), procs, &data_1r).unwrap();
    assert_eq!(reading.method.nucleus.as_deref(), Some("19F"));
    let pts = &reading.traces[0].points;
    assert_eq!(pts.len(), 3);
    assert!((pts[0].0 - (-50.0)).abs() < 1e-9);
    assert!((pts[0].1 - 50.0).abs() < 1e-9);
    assert!(pts[1].0 < pts[0].0);

    assert!(matches!(
        read_pdata(None, procs, &[0u8; 8]),
        Err(InstrumentError::InvalidBruker(_))
    ));
    assert!(matches!(
        read_pdata(None, procs, &[0u8; 7]),
        Err(InstrumentError::InvalidBruker(_))
    ));
    assert!(matches!(
        parse_procs("##$BYTORDP= 9\n"),
        Err(InstrumentError::InvalidBruker(_))
    ));
    assert!(matches!(
        read(Some("exp/1/pdata/1/1r"), &data_1r),
        Err(InstrumentError::InvalidBruker(_))
    ));
}

#[test]
fn openlab_dx_decodes_chemstation_130_single_delta_and_181_double_delta() {
    let mut ch130 = vec![0u8; 0x1800];
    ch130[0..4].copy_from_slice(b"\x03130");
    ch130[0x11A..0x11E].copy_from_slice(&0i32.to_be_bytes());
    ch130[0x11E..0x122].copy_from_slice(&60_000i32.to_be_bytes());
    ch130[0x127C..0x1284].copy_from_slice(&0.5_f64.to_be_bytes());
    write_ch_utf16(&mut ch130, 0x35A, "TBP-130");
    write_ch_utf16(&mut ch130, 0x104C, "mAU");
    write_ch_utf16(&mut ch130, 0x1075, "VWD1A, Wavelength=210 nm");
    ch130.extend_from_slice(&[0x10, 3]);
    ch130.extend_from_slice(&i16::MIN.to_be_bytes());
    ch130.extend_from_slice(&1000i32.to_be_bytes());
    ch130.extend_from_slice(&200i16.to_be_bytes());
    ch130.extend_from_slice(&(-400i16).to_be_bytes());
    ch130.push(0x00);

    let dx130 = build_zip(&[("chan1.CH", &ch130, true)]);
    let info130 = inspect_dx(&dx130).unwrap();
    assert_eq!(info130.sample_name.as_deref(), Some("TBP-130"));
    assert_eq!(info130.wavelength_nm, Some(210.0));
    assert_eq!(
        info130.traces[0].points,
        vec![(0.0, 500.0), (0.5, 600.0), (1.0, 400.0)]
    );

    let mut ch181 = vec![0u8; 0x1800];
    ch181[0..4].copy_from_slice(b"\x03181");
    ch181[0x11A..0x11E].copy_from_slice(&0.0_f32.to_be_bytes());
    ch181[0x11E..0x122].copy_from_slice(&60_000.0_f32.to_be_bytes());
    ch181[0x127C..0x1284].copy_from_slice(&1.0_f64.to_be_bytes());
    write_ch_utf16(&mut ch181, 0x104C, "pA");
    ch181.extend_from_slice(&i16::MAX.to_be_bytes());
    ch181.extend_from_slice(&0i16.to_be_bytes());
    ch181.extend_from_slice(&100i32.to_be_bytes());
    ch181.extend_from_slice(&10i16.to_be_bytes());
    ch181.extend_from_slice(&5i16.to_be_bytes());

    let dx181 = build_zip(&[("fid1.CH", &ch181, false)]);
    let info181 = inspect_dx(&dx181).unwrap();
    assert_eq!(info181.traces[0].y_unit, "pA");
    assert_eq!(
        info181.traces[0].points,
        vec![(0.0, 100.0), (0.5, 110.0), (1.0, 125.0)]
    );
}
