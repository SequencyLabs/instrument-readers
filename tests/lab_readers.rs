#![allow(clippy::unwrap_used, clippy::expect_used)]

use sequency_instruments::{
    InstrumentError, Technique, VendorFormat, chemstation, chromeleon, detect, newstyle, read,
    read_files, waters,
};

fn build_ch_130(
    sample: &str,
    method: &str,
    instrument: &str,
    signal: &str,
    deltas: &[i16],
) -> Vec<u8> {
    let mut buf = vec![0u8; 0x1800];
    buf[0] = 3;
    buf[1..4].copy_from_slice(b"130");
    buf[0x11A..0x11E].copy_from_slice(&0i32.to_be_bytes());
    buf[0x11E..0x122].copy_from_slice(&60_000i32.to_be_bytes());
    buf[0x127C..0x1284].copy_from_slice(&0.5_f64.to_be_bytes());

    let write_u16_pascal = |dst: &mut [u8], off: usize, s: &str| {
        let units: Vec<u16> = s.encode_utf16().collect();
        dst[off] = units.len() as u8;
        for (i, u) in units.iter().enumerate() {
            dst[off + 1 + i * 2..off + 3 + i * 2].copy_from_slice(&u.to_le_bytes());
        }
    };
    write_u16_pascal(&mut buf, 0x35A, sample);
    write_u16_pascal(&mut buf, 0xA0E, method);
    write_u16_pascal(&mut buf, 0xC11, instrument);
    write_u16_pascal(&mut buf, 0x104C, "mAU");
    write_u16_pascal(&mut buf, 0x1075, signal);

    buf.push(0x10);
    buf.push(deltas.len() as u8);
    for &d in deltas {
        buf.extend_from_slice(&d.to_be_bytes());
    }
    buf
}

fn build_chemstation_ms(scans: &[(u32, &[(u16, u16)])]) -> Vec<u8> {
    let mut buf = vec![0u8; 0x200];
    buf[0..4].copy_from_slice(&[0x01, 0x32, 0x00, 0x00]);
    let label = b"MSD Spectral File";
    buf[0x04] = label.len() as u8;
    buf[0x05..0x05 + label.len()].copy_from_slice(label);

    let date = b"16 Nov 26 09:28 am +0800";
    buf[0xB2] = date.len() as u8;
    buf[0xB3..0xB3 + date.len()].copy_from_slice(date);

    let method = b"LCMS_STD.M";
    buf[0xE4] = method.len() as u8;
    buf[0xE5..0xE5 + method.len()].copy_from_slice(method);

    let start_word = ((0x200 + 2) / 2) as u16;
    buf[0x10A..0x10C].copy_from_slice(&start_word.to_be_bytes());
    buf[0x116..0x11A].copy_from_slice(&(scans.len() as u32).to_be_bytes());

    for &(time_ms, pairs) in scans {
        buf.extend_from_slice(&[0u8; 2]);
        buf.extend_from_slice(&time_ms.to_be_bytes());
        buf.extend_from_slice(&[0u8; 6]);
        buf.extend_from_slice(&(pairs.len() as u16).to_be_bytes());
        buf.extend_from_slice(&[0u8; 4]);
        for &(mz_raw, int_enc) in pairs {
            buf.extend_from_slice(&mz_raw.to_be_bytes());
            buf.extend_from_slice(&int_enc.to_be_bytes());
        }
        buf.extend_from_slice(&[0u8; 10]);
    }
    buf
}

fn encode_chromeleon_varint(val: i64, out: &mut Vec<u8>) {
    let negative = val < 0;
    let mut mag = val.unsigned_abs();
    let mut first = (mag & 0x3F) as u8;
    if negative {
        first |= 0x40;
    }
    mag >>= 6;
    if mag != 0 {
        first |= 0x80;
    }
    out.push(first);
    while mag != 0 {
        let mut byte = (mag & 0x7F) as u8;
        mag >>= 7;
        if mag != 0 {
            byte |= 0x80;
        }
        out.push(byte);
    }
}

fn build_chromeleon_signal(
    name: &str,
    per_min: f64,
    per_unit: f64,
    series: &[(i64, i64)],
) -> Vec<u8> {
    let mut hdr_body = vec![0u8; 15 + name.len()];
    hdr_body[14] = name.len() as u8;
    hdr_body[15..].copy_from_slice(name.as_bytes());
    let mut out = Vec::new();
    out.extend_from_slice(b"SignHdr\0");
    out.extend_from_slice(&((16 + hdr_body.len()) as u64).to_le_bytes());
    out.extend_from_slice(&hdr_body);

    let mut pts_body = Vec::new();
    pts_body.extend_from_slice(&per_min.to_le_bytes());
    pts_body.extend_from_slice(&per_unit.to_le_bytes());
    for &(t, v) in series {
        encode_chromeleon_varint(t, &mut pts_body);
        encode_chromeleon_varint(v, &mut pts_body);
    }
    pts_body.push(0x40);

    out.extend_from_slice(b"PtsLL2Df");
    out.extend_from_slice(&((16 + pts_body.len()) as u64).to_le_bytes());
    out.extend_from_slice(&pts_body);
    out
}

fn build_waters_idx_record(offset: u32, pair_count: u32, time_min: f32) -> [u8; 22] {
    let mut rec = [0u8; 22];
    rec[0..4].copy_from_slice(&offset.to_le_bytes());
    rec[4..8].copy_from_slice(&(pair_count & 0x3F_FFFF).to_le_bytes());
    rec[12..16].copy_from_slice(&time_min.to_le_bytes());
    rec
}

#[test]
fn newstyle_nu3000_bundle_parses_traces_peaks_and_device_metadata() {
    let data_xml = b"\xef\xbb\xbf<Data unit=\"mAu\" count=\"5\">0.10 0.50 4.25 1.10 0.05</Data>";
    let icm_xml = br#"<Root>
  <Chrom samplingRate="0.01">
    <PeakGroup count="2">
      <Peak name="N.A." AREA="1250.5" AREA_PERCENT="25.0" HEIGHT="4.25">
        <Feature StartIndex="1" SummitIndex="2" EndIndex="4" />
      </Peak>
      <Peak name="ProductA" AREA="3751.5" AREA_PERCENT="75.0" HEIGHT="12.0">
        <Feature StartIndex="60" SummitIndex="120" EndIndex="180" />
      </Peak>
    </PeakGroup>
  </Chrom>
</Root>"#;
    let sys_xml = r#"<System Version="1.0">
  <Device Category="Pump" Name="50mL半制备输液泵" Manufactory="清博华" />
  <Device Category="Detector" Name="双波长紫外可见检测器" Manufactory="清博华" />
</System>"#
        .as_bytes();

    assert_eq!(
        detect(Some("Data/ch1.Data"), data_xml),
        VendorFormat::NewstyleNu3000
    );
    assert_eq!(
        detect(Some("Data/ch1.icm"), icm_xml),
        VendorFormat::NewstyleNu3000
    );

    let reading = read_files(&[
        ("Data/run1.Data", data_xml),
        ("Data/run1.icm", icm_xml),
        ("Data/run1.sys", sys_xml),
    ])
    .unwrap();

    assert_eq!(reading.vendor_format, VendorFormat::NewstyleNu3000);
    assert_eq!(reading.technique, Technique::Hplc);
    assert!(reading.instrument.as_deref().unwrap().contains("清博华"));
    assert_eq!(reading.traces.len(), 1);
    assert_eq!(reading.traces[0].points.len(), 5);
    assert_eq!(reading.traces[0].y_unit, "mAu");
    assert!((reading.traces[0].points[2].1 - 4.25).abs() < 1e-9);
    assert_eq!(reading.peaks.len(), 2);
    assert_eq!(reading.peaks[0].compound, None);
    assert_eq!(reading.peaks[1].compound.as_deref(), Some("ProductA"));
    assert!((reading.peaks[1].rt_min.unwrap() - 0.2).abs() < 1e-6);

    let single = read(Some("ch1.Data"), data_xml).unwrap();
    assert_eq!(single.traces[0].points.len(), 5);

    assert!(matches!(
        newstyle::read_data(b"<Other>1 2</Other>"),
        Err(InstrumentError::InvalidNewstyle(_))
    ));
}

#[test]
fn chemstation_standalone_ch_ms_and_d_bundle_parse_correctly() {
    let ch = build_ch_130(
        "RXN-2026-01",
        "D:\\METHODS\\LC254.M",
        "Agilent 1200",
        "VWD1A, Wavelength=254 nm",
        &[10, 20, -5, 15],
    );
    let ms = build_chemstation_ms(&[
        (30_000, &[(4000, 100), (4200, (1 << 14) | 25)]),
        (60_000, &[(4000, (2 << 14) | 10)]),
    ]);
    let sample_xml = b"<?xml version=\"1.0\"?><Sample><Name>RXN-2026-01-XML</Name></Sample>";

    assert_eq!(detect(Some("VWD1A.ch"), &ch), VendorFormat::ChemStation);
    assert_eq!(detect(Some("MSD1.MS"), &ms), VendorFormat::ChemStation);

    let ch_only = read(Some("VWD1A.ch"), &ch).unwrap();
    assert_eq!(ch_only.technique, Technique::Hplc);
    assert_eq!(ch_only.sample_name.as_deref(), Some("RXN-2026-01"));
    assert_eq!(ch_only.method.wavelength_nm, Some(254.0));
    assert_eq!(ch_only.traces[0].points.len(), 4);

    let ms_only = read(Some("MSD1.MS"), &ms).unwrap();
    assert_eq!(ms_only.technique, Technique::LcMs);
    assert_eq!(ms_only.traces[0].points.len(), 2);
    assert!((ms_only.traces[0].points[0].0 - 0.5).abs() < 1e-9);
    assert!((ms_only.traces[0].points[0].1 - 300.0).abs() < 1e-9);
    assert!((ms_only.traces[0].points[1].1 - 640.0).abs() < 1e-9);

    let bundle = read_files(&[
        ("001-0101.D/SAMPLE.XML", sample_xml),
        ("001-0101.D/VWD1A.ch", &ch),
        ("001-0101.D/MSD1.MS", &ms),
    ])
    .unwrap();
    assert_eq!(bundle.technique, Technique::LcMs);
    assert_eq!(bundle.sample_name.as_deref(), Some("RXN-2026-01-XML"));
    assert_eq!(bundle.traces.len(), 2);

    assert!(matches!(
        chemstation::read_ms(b"bad"),
        Err(InstrumentError::InvalidChemStation(_))
    ));
}

#[test]
fn chromeleon_signal_and_bundle_decode_ptsll2df_and_skip_unsupported_aux() {
    let uv = build_chromeleon_signal(
        "UV_VIS_1",
        60.0,
        10.0,
        &[(0, 100), (1, 20), (0, -5), (0, 5)],
    );
    let mut unsupported = Vec::new();
    unsupported.extend_from_slice(b"SignHdr\0");
    let mut hdr_body = vec![0u8; 15 + 13];
    hdr_body[14] = 13;
    hdr_body[15..].copy_from_slice(b"Pump_Pressure");
    unsupported.extend_from_slice(&((16 + hdr_body.len()) as u64).to_le_bytes());
    unsupported.extend_from_slice(&hdr_body);
    unsupported.extend_from_slice(b"PtsLDiff");
    unsupported.extend_from_slice(&20u64.to_le_bytes());
    unsupported.extend_from_slice(&[0u8; 4]);

    assert_eq!(chromeleon::signal_name(&uv).as_deref(), Some("UV_VIS_1"));
    assert_eq!(
        detect(Some("095216244.raw"), &uv),
        VendorFormat::ChromeleonRaw
    );

    let reading = read(Some("095216244.raw"), &uv).unwrap();
    assert_eq!(reading.vendor_format, VendorFormat::ChromeleonRaw);
    assert_eq!(reading.traces[0].label, "UV_VIS_1");
    assert_eq!(reading.traces[0].y_unit, "mAU");
    let pts = &reading.traces[0].points;
    assert_eq!(pts.len(), 4);
    assert!((pts[0].0 - 0.0).abs() < 1e-9);
    assert!((pts[3].0 - 3.0 / 60.0).abs() < 1e-9);
    assert!((pts[0].1 - 10.0).abs() < 1e-9);
    assert!((pts[1].1 - 12.0).abs() < 1e-9);
    assert!((pts[2].1 - 13.5).abs() < 1e-9);
    assert!((pts[3].1 - 15.5).abs() < 1e-9);

    let bundle = read_files(&[("pump.raw", &unsupported), ("uv.raw", &uv)]).unwrap();
    assert_eq!(bundle.traces.len(), 1);
    assert_eq!(bundle.traces[0].label, "UV_VIS_1");

    assert!(matches!(
        chromeleon::read_signal(&unsupported),
        Err(InstrumentError::InvalidChromeleon(_))
    ));
}

#[test]
fn waters_masslynx_bundle_decodes_header_polarities_and_2_4_6_byte_functions() {
    let header = b"$$ Version: 01.00\r\n$$ Acquired Name: 20261126-TEST-01\r\n$$ Acquired Date: 26-Nov-2026\r\n$$ Acquired Time: 13:39:18\r\n$$ Instrument: SQ Detector 2\r\n";
    let extern_inf = b"Instrument Parameters\r\nPolarity\t\t\tES+\r\nInstrument Parameters\r\nPolarity\t\t\tES-\r\n";

    let mut f1_dat = Vec::new();
    f1_dat.extend_from_slice(&100i16.to_le_bytes());
    f1_dat.extend_from_slice(&0u32.to_le_bytes());
    f1_dat.extend_from_slice(&25i16.to_le_bytes());
    f1_dat.extend_from_slice(&1u32.to_le_bytes());
    let mut f1_idx = Vec::new();
    f1_idx.extend_from_slice(&build_waters_idx_record(0, 1, 0.5));
    f1_idx.extend_from_slice(&build_waters_idx_record(6, 1, 1.0));

    let mut f2_dat = Vec::new();
    let sample2: u16 = (50 << 3) | 1;
    f2_dat.extend_from_slice(&sample2.to_le_bytes());
    let f2_idx = build_waters_idx_record(0, 1, 0.75);

    let mut f3_dat = Vec::new();
    let sample4: u32 = (10 << 22) | 1024;
    f3_dat.extend_from_slice(&sample4.to_le_bytes());
    let f3_idx = build_waters_idx_record(0, 1, 0.25);

    let reading = read_files(&[
        ("run.raw/_HEADER.TXT", header),
        ("run.raw/_extern.inf", extern_inf),
        ("run.raw/_FUNC001.DAT", &f1_dat),
        ("run.raw/_FUNC001.IDX", &f1_idx),
        ("run.raw/_FUNC002.DAT", &f2_dat),
        ("run.raw/_FUNC002.IDX", &f2_idx),
        ("run.raw/_FUNC003.DAT", &f3_dat),
        ("run.raw/_FUNC003.IDX", &f3_idx),
    ])
    .unwrap();

    assert_eq!(reading.vendor_format, VendorFormat::WatersRaw);
    assert_eq!(reading.technique, Technique::LcMs);
    assert_eq!(reading.sample_name.as_deref(), Some("20261126-TEST-01"));
    assert_eq!(reading.instrument.as_deref(), Some("SQ Detector 2"));
    assert_eq!(reading.acquired_at.as_deref(), Some("26-Nov-2026 13:39:18"));
    assert_eq!(reading.traces.len(), 3);
    assert_eq!(reading.traces[0].label, "ESI+ TIC (Func 1)");
    assert!((reading.traces[0].points[0].1 - 100.0).abs() < 1e-9);
    assert!((reading.traces[0].points[1].1 - 100.0).abs() < 1e-9);
    assert_eq!(reading.traces[1].label, "ESI- TIC (Func 2)");
    assert!((reading.traces[1].points[0].1 - 200.0).abs() < 1e-9);
    assert_eq!(reading.traces[2].label, "UV (Func 3)");
    assert!((reading.traces[2].points[0].1 - 1.0).abs() < 1e-9);

    assert!(matches!(
        waters::read_header(b"not a waters header"),
        Err(InstrumentError::InvalidWaters(_))
    ));
}
