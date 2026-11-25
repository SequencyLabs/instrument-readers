//! Integration tests for the pure-Rust OLE2 Compound File reader and Shimadzu `.lcd` parser.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use encoding_rs::GB18030;
use sequency_instruments::ole::{OLE_MAGIC, OleFile, place};
use sequency_instruments::shimadzu::{
    chromatograms, decode_deltas, parse_filetime, parse_inj_vol, parse_peak_table,
    parse_properties, parse_trace, read_lcd, read_lcm, read_run,
};
use sequency_instruments::{GradientStep, InstrumentError, Technique, VendorFormat};

const END_OF_CHAIN: u32 = 0xFFFF_FFFE;
const FREE_SECT: u32 = 0xFFFF_FFFF;
const FAT_SECT: u32 = 0xFFFF_FFFD;

/// Builds a valid OLE2 compound file with interleaved 512-byte regular sectors and 64-byte mini-sectors.
///
/// Mirrors `packages/lab-data/tests/lab_files.py::compound_file` so no multi-sector stream is contiguous.
fn build_compound_file(streams: &[(&str, &[u8])]) -> Vec<u8> {
    let sector = 512usize;
    let mini = 64usize;

    let small: Vec<(&str, &[u8])> = streams
        .iter()
        .copied()
        .filter(|(_, d)| d.len() < 4096)
        .collect();
    let big: Vec<(&str, &[u8])> = streams
        .iter()
        .copied()
        .filter(|(_, d)| d.len() >= 4096)
        .collect();

    let mini_counts: Vec<usize> = small.iter().map(|(_, d)| d.len().div_ceil(mini)).collect();
    let mini_chains = interleave(&mini_counts);
    let total_mini: usize = mini_chains.iter().map(Vec::len).sum();
    let mut minifat = vec![FREE_SECT; total_mini];
    link_chains(&mini_chains, &mut minifat);

    let mut ministream = vec![0u8; total_mini * mini];
    for (chain, (_, data)) in mini_chains.iter().zip(small.iter()) {
        for (i, &m) in chain.iter().enumerate() {
            let m_idx = m as usize;
            let src_start = i * mini;
            let src_end = (src_start + mini).min(data.len());
            let slice = &data[src_start..src_end];
            ministream[m_idx * mini..m_idx * mini + slice.len()].copy_from_slice(slice);
        }
    }

    let mut storage_set = std::collections::BTreeSet::new();
    for &(path, _) in streams {
        let parts: Vec<&str> = path.split('/').collect();
        for depth in 1..parts.len() {
            storage_set.insert(parts[..depth].join("/"));
        }
    }
    let storages: Vec<String> = storage_set.into_iter().collect();

    let mut entries: Vec<String> = vec![String::new()];
    entries.extend(storages.iter().cloned());
    entries.extend(streams.iter().map(|&(p, _)| p.to_owned()));

    let mut content: Vec<&[u8]> = vec![&ministream];
    content.extend(big.iter().map(|&(_, d)| d));
    let reg_counts: Vec<usize> = content.iter().map(|d| d.len().div_ceil(sector)).collect();
    let chains = interleave(&reg_counts);
    let used: usize = chains.iter().map(Vec::len).sum();

    let dir_sectors = (entries.len() * 128).div_ceil(sector).max(1);
    let minifat_sectors = (minifat.len() * 4).div_ceil(sector);
    let mut fat_sectors = 1usize;
    while fat_sectors * (sector / 4) < used + dir_sectors + minifat_sectors + fat_sectors {
        fat_sectors += 1;
    }

    let first_dir = used;
    let first_minifat = used + dir_sectors;
    let first_fat = first_minifat + minifat_sectors;
    let mut fat = vec![FREE_SECT; fat_sectors * (sector / 4)];
    link_chains(&chains, &mut fat);
    let dir_chain: Vec<u32> = (first_dir..first_minifat).map(|n| n as u32).collect();
    let mf_chain: Vec<u32> = (first_minifat..first_fat).map(|n| n as u32).collect();
    link_chains(&[dir_chain, mf_chain], &mut fat);
    for slot in &mut fat[first_fat..first_fat + fat_sectors] {
        *slot = FAT_SECT;
    }

    let mut body = vec![0u8; (first_fat + fat_sectors) * sector];
    for (chain, data) in chains.iter().zip(content.iter()) {
        for (i, &s) in chain.iter().enumerate() {
            let s_idx = s as usize;
            let src_start = i * sector;
            let src_end = (src_start + sector).min(data.len());
            let slice = &data[src_start..src_end];
            body[s_idx * sector..s_idx * sector + slice.len()].copy_from_slice(slice);
        }
    }

    let mut starts = std::collections::BTreeMap::new();
    starts.insert(
        String::new(),
        if minifat.is_empty() {
            END_OF_CHAIN
        } else {
            chains[0][0]
        },
    );
    for (i, &(p, _)) in small.iter().enumerate() {
        let start = mini_chains[i].first().copied().unwrap_or(END_OF_CHAIN);
        starts.insert(p.to_owned(), start);
    }
    for (i, &(p, _)) in big.iter().enumerate() {
        let start = chains[1 + i].first().copied().unwrap_or(END_OF_CHAIN);
        starts.insert(p.to_owned(), start);
    }

    let mut sizes = std::collections::BTreeMap::new();
    sizes.insert(String::new(), ministream.len() as u32);
    for &(p, d) in streams {
        sizes.insert(p.to_owned(), d.len() as u32);
    }

    let index_map: std::collections::BTreeMap<&str, u32> = entries
        .iter()
        .enumerate()
        .map(|(i, p)| (p.as_str(), i as u32))
        .collect();

    let mut kids: std::collections::BTreeMap<String, Vec<String>> =
        std::collections::BTreeMap::new();
    kids.insert(String::new(), Vec::new());
    for st in &storages {
        kids.insert(st.clone(), Vec::new());
    }
    for p in &entries[1..] {
        let parent = p.rsplit_once('/').map_or("", |(par, _)| par);
        kids.get_mut(parent).unwrap().push(p.clone());
    }

    let mut directory = Vec::with_capacity(dir_sectors * sector);
    for p in &entries {
        let (parent, leaf) = match p.rsplit_once('/') {
            Some((par, l)) => (par, l),
            None => ("", p.as_str()),
        };
        let right = if p.is_empty() {
            FREE_SECT
        } else {
            let sibs = &kids[parent];
            let pos = sibs.iter().position(|x| x == p).unwrap();
            sibs.get(pos + 1)
                .map_or(FREE_SECT, |next| index_map[next.as_str()])
        };
        let child = kids
            .get(p)
            .and_then(|c| c.first())
            .map_or(FREE_SECT, |first| index_map[first.as_str()]);
        let display_name = if p.is_empty() { "Root Entry" } else { leaf };
        let kind: u8 = if p.is_empty() {
            5
        } else if kids.contains_key(p) {
            1
        } else {
            2
        };

        let mut slot = [0u8; 128];
        let mut utf16 = encode_utf16le(display_name);
        utf16.extend_from_slice(&[0, 0]);
        slot[..utf16.len()].copy_from_slice(&utf16);
        slot[64..66].copy_from_slice(&(utf16.len() as u16).to_le_bytes());
        slot[66] = kind;
        slot[67] = 1;
        slot[68..72].copy_from_slice(&FREE_SECT.to_le_bytes());
        slot[72..76].copy_from_slice(&right.to_le_bytes());
        slot[76..80].copy_from_slice(&child.to_le_bytes());
        let start_sec = if kind == 1 {
            0
        } else {
            starts.get(p).copied().unwrap_or(END_OF_CHAIN)
        };
        let sz = sizes.get(p).copied().unwrap_or(0);
        slot[116..120].copy_from_slice(&start_sec.to_le_bytes());
        slot[120..124].copy_from_slice(&sz.to_le_bytes());
        directory.extend_from_slice(&slot);
    }
    while directory.len() < dir_sectors * sector {
        let mut empty = [0u8; 128];
        empty[68..72].copy_from_slice(&FREE_SECT.to_le_bytes());
        empty[72..76].copy_from_slice(&FREE_SECT.to_le_bytes());
        empty[76..80].copy_from_slice(&FREE_SECT.to_le_bytes());
        directory.extend_from_slice(&empty);
    }

    body[first_dir * sector..first_dir * sector + directory.len()].copy_from_slice(&directory);
    for (i, &m) in minifat.iter().enumerate() {
        let off = first_minifat * sector + i * 4;
        body[off..off + 4].copy_from_slice(&m.to_le_bytes());
    }
    for (i, &f) in fat.iter().enumerate() {
        let off = first_fat * sector + i * 4;
        body[off..off + 4].copy_from_slice(&f.to_le_bytes());
    }

    let mut header = vec![0u8; 512];
    header[..8].copy_from_slice(&OLE_MAGIC);
    header[24..26].copy_from_slice(&0x003Eu16.to_le_bytes());
    header[26..28].copy_from_slice(&3u16.to_le_bytes());
    header[28..30].copy_from_slice(&0xFFFEu16.to_le_bytes());
    header[30..32].copy_from_slice(&9u16.to_le_bytes());
    header[32..34].copy_from_slice(&6u16.to_le_bytes());
    header[44..48].copy_from_slice(&(fat_sectors as u32).to_le_bytes());
    header[48..52].copy_from_slice(&(first_dir as u32).to_le_bytes());
    header[56..60].copy_from_slice(&4096u32.to_le_bytes());
    let mf_head = if minifat.is_empty() {
        END_OF_CHAIN
    } else {
        first_minifat as u32
    };
    header[60..64].copy_from_slice(&mf_head.to_le_bytes());
    header[64..68].copy_from_slice(&(minifat_sectors as u32).to_le_bytes());
    header[68..72].copy_from_slice(&END_OF_CHAIN.to_le_bytes());
    header[72..76].copy_from_slice(&0u32.to_le_bytes());

    for idx in 0..109 {
        let sec = if idx < fat_sectors {
            (first_fat + idx) as u32
        } else {
            FREE_SECT
        };
        let off = 76 + idx * 4;
        header[off..off + 4].copy_from_slice(&sec.to_le_bytes());
    }

    header.extend_from_slice(&body);
    header
}

fn interleave(counts: &[usize]) -> Vec<Vec<u32>> {
    let mut chains = vec![Vec::new(); counts.len()];
    let max_count = counts.iter().copied().max().unwrap_or(0);
    let mut next = 0u32;
    for step in 0..max_count {
        for (idx, &count) in counts.iter().enumerate() {
            if step < count {
                chains[idx].push(next);
                next += 1;
            }
        }
    }
    chains
}

fn link_chains(chains: &[Vec<u32>], table: &mut [u32]) {
    for chain in chains {
        for (i, &cur) in chain.iter().enumerate() {
            let next = chain.get(i + 1).copied().unwrap_or(END_OF_CHAIN);
            table[cur as usize] = next;
        }
    }
}

fn encode_utf16le(text: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(text.len() * 2);
    for unit in text.encode_utf16() {
        out.extend_from_slice(&unit.to_le_bytes());
    }
    out
}

fn stox_gb18030(text: &str) -> String {
    let (bytes, _, _) = GB18030.encode(text);
    let mut out = String::from("@StoX@");
    for &b in bytes.as_ref() {
        out.push_str(&format!("{b:02X}"));
    }
    out
}

fn pack_peak_record(rt_ms: i32, area: f64, height: f64, record_size: usize) -> Vec<u8> {
    let mut rec = vec![0u8; record_size];
    rec[0..4].copy_from_slice(&0x8010u32.to_le_bytes());
    rec[4..8].copy_from_slice(&rt_ms.to_le_bytes());
    rec[8..16].copy_from_slice(&area.to_le_bytes());
    rec[16..24].copy_from_slice(&0.0_f64.to_le_bytes());
    rec[24..32].copy_from_slice(&height.to_le_bytes());
    rec
}

fn pack_trace_stream(interval_ms: u32, sub_blocks: &[&[u8]]) -> Vec<u8> {
    let mut body = Vec::new();
    let mut total_points = 0u32;
    for &block in sub_blocks {
        let len = block.len() as u16;
        total_points += decode_deltas(block).unwrap().len() as u32;
        body.extend_from_slice(&len.to_le_bytes());
        body.extend_from_slice(block);
        body.extend_from_slice(&len.to_le_bytes());
    }
    let total_len = (24 + body.len()) as i32;
    let mut header = vec![0u8; 24];
    header[0..4].copy_from_slice(b"RC\0\0");
    header[4..8].copy_from_slice(&interval_ms.to_le_bytes());
    header[8..12].copy_from_slice(&total_points.to_le_bytes());
    header[12..16].copy_from_slice(&total_len.to_le_bytes());
    header.extend_from_slice(&body);
    header
}

fn pack_status_record(cf: f64, vf: f64, unit: &[u8]) -> Vec<u8> {
    let mut rec = vec![0u8; 64];
    rec[0..4].copy_from_slice(&1i32.to_le_bytes());
    rec[8..12].copy_from_slice(&1000i32.to_le_bytes());
    rec[16..24].copy_from_slice(&cf.to_le_bytes());
    rec[24..32].copy_from_slice(&1.0_f64.to_le_bytes());
    rec[32..40].copy_from_slice(&vf.to_le_bytes());
    rec[40..40 + unit.len()].copy_from_slice(unit);
    rec
}

#[test]
fn ole_reads_interleaved_mini_and_regular_streams_and_maps_extents() {
    let mini_a: Vec<u8> = (0..200u8).collect();
    let mini_b: Vec<u8> = (0..150u8).map(|b| b ^ 0xAA).collect();
    let big_a: Vec<u8> = (0..5000usize).map(|i| (i % 251) as u8).collect();
    let big_b: Vec<u8> = (0..4300usize).map(|i| ((i * 7) % 253) as u8).collect();

    let raw = build_compound_file(&[
        ("Alpha/Small1", &mini_a),
        ("Alpha/Nested/Small2", &mini_b),
        ("Beta/Large1", &big_a),
        ("Beta/Large2", &big_b),
        ("Empty", &[]),
    ]);

    let ole = OleFile::parse(&raw).unwrap();
    assert_eq!(
        ole.stream("Alpha/Small1").unwrap().as_ref(),
        mini_a.as_slice()
    );
    assert_eq!(
        ole.stream("Alpha/Nested/Small2").unwrap().as_ref(),
        mini_b.as_slice()
    );
    assert_eq!(
        ole.stream("Beta/Large1").unwrap().as_ref(),
        big_a.as_slice()
    );
    assert_eq!(
        ole.stream("Beta/Large2").unwrap().as_ref(),
        big_b.as_slice()
    );
    assert_eq!(ole.stream("Empty").unwrap().as_ref(), b"");
    assert!(!ole.exists("Missing"));

    let extents = ole.extents("Alpha/Small1").unwrap();
    assert!(extents.len() > 1);
    let mapped = place(extents, 50, 80);
    let reconstructed: Vec<u8> = mapped
        .into_iter()
        .flat_map(|(a, b)| raw[a..b].iter().copied())
        .collect();
    assert_eq!(reconstructed, &mini_a[50..80]);

    assert!(!ole.unallocated().is_empty());
}

#[test]
fn ole_rejects_looping_sector_chains_and_invalid_headers() {
    assert!(matches!(
        OleFile::parse(b"not an ole2 file"),
        Err(InstrumentError::InvalidOle(_))
    ));

    let big = vec![0x42u8; 4500];
    let mut raw = build_compound_file(&[("StreamA", &big)]);
    let first_fat_sector = u32::from_le_bytes(raw[76..80].try_into().unwrap()) as usize;
    let fat_byte_offset = (1 + first_fat_sector) * 512;
    raw[fat_byte_offset..fat_byte_offset + 4].copy_from_slice(&0u32.to_le_bytes());

    assert!(matches!(
        OleFile::parse(&raw),
        Err(InstrumentError::CorruptChain(0))
    ));
}

#[test]
fn shimadzu_file_property_and_filetime_match_labsolutions_encoding() {
    let mut stream = vec![0u8, 0u8];
    let xml = format!(
        concat!(
            "<?xml version=\"1.0\"?>\n<FileProperty><szVersion>@StoX@352E3031</szVersion></FileProperty>",
            "<?xml version=\"1.0\"?>\n<SampleInfo><smpl_name>{sample}</smpl_name>",
            "<inj_vol>@FtoX@40400000</inj_vol>",
            "<dwLowDateTime>-2085920494</dwLowDateTime><dwHighDateTime>31038886</dwHighDateTime></SampleInfo>",
            "<?xml version=\"1.0\"?>\n<SampleInfoFile><datafile>{datafile}</datafile>",
            "<methodfile>{method}</methodfile></SampleInfoFile>",
            "<?xml version=\"1.0\"?>\n<DataFileProperty><szVersion>{ver2}</szVersion></DataFileProperty>"
        ),
        sample = stox_gb18030("TBP-230538"),
        datafile = stox_gb18030("E:/2023/20230613/TBP-230538.lcd"),
        method = stox_gb18030("E:/2023/LC-通用方法-95.lcm"),
        ver2 = stox_gb18030("5.106 SP1"),
    );
    stream.extend_from_slice(xml.as_bytes());

    let props = parse_properties(&stream);
    assert_eq!(props["smpl_name"], "TBP-230538");
    assert_eq!(props["methodfile"], "E:/2023/LC-通用方法-95.lcm");
    assert_eq!(props["szVersion"], "5.01");
    assert_eq!(
        parse_filetime(
            props.get("dwLowDateTime").map(String::as_str),
            props.get("dwHighDateTime").map(String::as_str)
        )
        .as_deref(),
        Some("2023-06-13T03:24:08.331905+00:00")
    );
    assert_eq!(parse_filetime(Some("0"), Some("0")), None);
    assert_eq!(parse_inj_vol(&props["inj_vol"]), Some(3.0));
}

#[test]
fn shimadzu_peak_tables_decode_ver1_and_legacy_layouts() {
    let rows = [(550_278i32, 3000.0, 800.0), (652_800i32, 47_000.0, 9000.0)];

    let mut ver1 = b"VER1".to_vec();
    ver1.extend_from_slice(&2i32.to_le_bytes());
    ver1.extend_from_slice(&[0u8; 12]);
    for &(rt, area, height) in &rows {
        ver1.extend_from_slice(&pack_peak_record(rt, area, height, 1056));
    }

    let mut legacy = Vec::new();
    legacy.extend_from_slice(&2i32.to_le_bytes());
    legacy.extend_from_slice(&0i32.to_le_bytes());
    for &(rt, area, height) in &rows {
        legacy.extend_from_slice(&pack_peak_record(rt, area, height, 300));
    }

    for table_bytes in [&ver1, &legacy] {
        let peaks = parse_peak_table(table_bytes);
        assert_eq!(peaks.len(), 2);
        assert_eq!(
            (peaks[0].rt_min, peaks[0].area_pct, peaks[0].height),
            (Some(9.171), Some(6.0), Some(800.0))
        );
        assert_eq!(
            (peaks[1].rt_min, peaks[1].area_pct, peaks[1].height),
            (Some(10.88), Some(94.0), Some(9000.0))
        );
    }

    assert!(parse_peak_table(b"").is_empty());
    assert!(parse_peak_table(b"VER1").is_empty());
    let mut truncated = Vec::new();
    truncated.extend_from_slice(&1000i32.to_le_bytes());
    truncated.extend_from_slice(&[0u8; 44]);
    assert!(parse_peak_table(&truncated).is_empty());
}

#[test]
fn shimadzu_deltas_and_trace_sub_blocks_match_reference() {
    let sub1 = [0x05u8, 0x1E, 0x82, 0x21, 0x29, 0x3F, 0xF6];
    assert_eq!(decode_deltas(&sub1).unwrap(), vec![5, 3, 300, 290]);
    assert_eq!(
        decode_deltas(&[0x6F, 0xFF, 0xFF, 0xFF]).unwrap(),
        vec![0x0FFF_FFFF]
    );
    assert!(matches!(
        decode_deltas(&[0x80, 0, 0, 0, 0]),
        Err(InstrumentError::InvalidShimadzu(_))
    ));

    let stream = pack_trace_stream(500, &[&sub1, &[0x07]]);
    let (minutes, values) = parse_trace(&stream).unwrap();
    assert_eq!(values, vec![5, 3, 300, 290, 7]);
    for (i, &m) in minutes.iter().enumerate() {
        assert!((m - (i as f64 / 120.0)).abs() < 1e-12);
    }
}

#[test]
fn shimadzu_full_lcd_extracts_metadata_method_gradient_peaks_and_chromatograms() {
    let prop = format!(
        concat!(
            "<?xml version=\"1.0\"?><SampleInfo><smpl_name>{sample}</smpl_name>",
            "<inj_vol>@FtoX@40400000</inj_vol>",
            "<dwLowDateTime>-2085920494</dwLowDateTime><dwHighDateTime>31038886</dwHighDateTime></SampleInfo>",
            "<?xml version=\"1.0\"?><SampleInfoFile><methodfile>{method}</methodfile></SampleInfoFile>"
        ),
        sample = stox_gb18030("TBP-230538"),
        method = stox_gb18030("E:\\2023\\LC-通用方法-95.lcm"),
    );

    let config_xml = "<GUD Type=\"SystemParameter\"><ModelName>CBM-20A</ModelName></GUD>";
    let method_xml = concat!(
        "<GUD Type=\"MethodParameter\">",
        "<UPD ID=\"PumpA_MobilePhaseName#1A\"><Val>0.05%磷酸水溶液</Val><Type>8</Type></UPD>",
        "<UPD ID=\"PumpA_MobilePhaseName#2A\"><Val>乙腈</Val><Type>8</Type></UPD>",
        "<UPD ID=\"PumpA_Flow\"><Val>0000803F</Val><Type>4</Type></UPD>",
        "<UPD ID=\"OVEN$use\"><Val>1</Val><Type>2</Type></UPD>",
        "<UPD ID=\"OVEN$ovent\"><Val>0000F041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"PumpA_BConc\"><Val>0000A040</Val><Type>4</Type></UPD>",
        "<TP PgN=\"SCL$TIMEPROG\">",
        "<TPLD><UPD ID=\"TMPGtime\"><Val>00007041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"TMPGcode\"><Val>5</Val><Type>2</Type></UPD>",
        "<UPD ID=\"TMPGdata\"><Val>0000B442</Val><Type>4</Type></UPD></TPLD>",
        "<TPLD><UPD ID=\"TMPGtime\"><Val>0000B841</Val><Type>4</Type></UPD>",
        "<UPD ID=\"TMPGcode\"><Val>97</Val><Type>2</Type></UPD>",
        "<UPD ID=\"TMPGdata\"><Val>00000000</Val><Type>4</Type></UPD></TPLD>",
        "</TP></GUD>"
    );
    let pda_xml = "<GUD><UPD ID=\"2DCH_Wavelength#1\"><Val>25400</Val><Type>2</Type></UPD></GUD>";

    let signal = pack_trace_stream(240, &[&[0x2F, 0xA0]]);
    let status = pack_status_record(1.0, 1000.0, b"mAU");
    let mut pt = b"VER1".to_vec();
    pt.extend_from_slice(&2i32.to_le_bytes());
    pt.extend_from_slice(&[0u8; 12]);
    pt.extend_from_slice(&pack_peak_record(550_278, 3000.0, 800.0, 64));
    pt.extend_from_slice(&pack_peak_record(652_800, 47_000.0, 9000.0, 64));

    let cfg_u16 = encode_utf16le(config_xml);
    let mth_u16 = encode_utf16le(method_xml);
    let pda_u16 = encode_utf16le(pda_xml);

    let lcd = build_compound_file(&[
        ("File Property", prop.as_bytes()),
        ("GUMM_Information/ShimadzuLC.1/GUM.CONFIG", &cfg_u16),
        ("GUMM_Information/ShimadzuLC.1/GUC.1.METHOD", &mth_u16),
        ("GUMM_Information/ShimadzuPDA.1/PDA.1.METHOD", &pda_u16),
        ("PDA 3D Raw Data/Max Plot", &signal),
        ("PDA 3D Raw Data/Max Plot Status", &status),
        ("LSS Raw Data/Chromatogram Ch1", &signal),
        ("LSS Data Processing/PT-PDA.1", &[]),
        ("LSS Data Processing Original/PT-PDA.1", &pt),
    ]);

    let reading = read_lcd(&lcd).unwrap();
    assert_eq!(reading.technique, Technique::Hplc);
    assert_eq!(reading.vendor_format, VendorFormat::ShimadzuLcd);
    assert_eq!(
        reading.acquired_at.as_deref(),
        Some("2023-06-13T03:24:08.331905+00:00")
    );
    assert_eq!(reading.sample_name.as_deref(), Some("TBP-230538"));
    assert_eq!(reading.instrument.as_deref(), Some("Shimadzu CBM-20A"));
    assert_eq!(reading.method.name.as_deref(), Some("LC-通用方法-95.lcm"));
    assert_eq!(reading.method.solvents, vec!["0.05%磷酸水溶液", "乙腈"]);
    assert_eq!(reading.method.flow_rate_ml_min, Some(1.0));
    assert_eq!(reading.method.column_temp_c, Some(30.0));
    assert_eq!(reading.method.injection_volume_ul, Some(3.0));
    assert_eq!(reading.method.wavelength_nm, Some(254.0));
    assert_eq!(
        reading.method.gradient,
        vec![
            GradientStep {
                time_min: 0.0,
                b_pct: 5.0,
                c_pct: 0.0,
                d_pct: 0.0,
            },
            GradientStep {
                time_min: 15.0,
                b_pct: 90.0,
                c_pct: 0.0,
                d_pct: 0.0,
            },
            GradientStep {
                time_min: 23.0,
                b_pct: 90.0,
                c_pct: 0.0,
                d_pct: 0.0,
            },
        ]
    );
    assert_eq!(reading.method.gradient[0].a_pct(), 95.0);

    let main = reading.main_peak().unwrap();
    assert_eq!((main.rt_min, main.area_pct), (Some(10.88), Some(94.0)));

    assert_eq!(reading.traces.len(), 1);
    assert_eq!(reading.traces[0].label, "Max absorbance (mAU)");
    assert_eq!(reading.traces[0].x_unit, "min");
    assert_eq!(reading.traces[0].y_unit, "mAU");
    assert_eq!(reading.traces[0].points[0], (0.0, 4.0));
}

#[test]
fn shimadzu_without_pda_reads_sorted_detector_channels() {
    let signal = pack_trace_stream(240, &[&[0x2F, 0xA0]]);
    let mut status = pack_status_record(0.5, 1000.0, b"mV");
    status.extend_from_slice(&vec![0u8; 64 * 8]);
    status.extend_from_slice(&pack_status_record(1.0, 1.0, b"uRIU"));

    let lcd = build_compound_file(&[
        ("PDA 3D Raw Data/Max Plot", &[]),
        ("LSS Raw Data/Chromatogram Ch2", &[]),
        ("LSS Raw Data/Chromatogram Ch10", &signal),
        ("LSS Raw Data/Chromatogram Ch1", &signal),
        ("LSS Raw Data/Chromatogram Status", &status),
    ]);

    let traces = chromatograms(&lcd).unwrap();
    let summary: Vec<(&str, &str, f64)> = traces
        .iter()
        .map(|t| (t.label.as_str(), t.y_unit.as_str(), t.points[0].1))
        .collect();
    assert_eq!(
        summary,
        vec![("Ch1 (mV)", "mV", 2.0), ("Ch10 (uRIU)", "uRIU", 4000.0)]
    );
}

#[test]
fn shimadzu_cleaned_files_drop_masked_names_and_decode_cbm20a_lpge_gradients() {
    let prop = b"<?xml version=\"1.0\"?><SampleInfo><smpl_name>XXXXXXXXXX</smpl_name></SampleInfo>";
    let method_xml = concat!(
        "<GUD Type=\"MethodParameter\">",
        "<UPD ID=\"PumpA_MobilePhaseName#1A\"><Val>XXXX</Val><Type>8</Type></UPD>",
        "<UPD ID=\"LPGE$tflow\"><Val>0000003F</Val><Type>4</Type></UPD>",
        "<UPD ID=\"LPGE$bconc\"><Val>00002041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"LPGE$cconc\"><Val>0000A040</Val><Type>4</Type></UPD>",
        "<TP PgN=\"SCL$TIMEPROG\">",
        "<TPLD><UPD ID=\"TMPGtime\"><Val>00002041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"TMPGcode\"><Val>5</Val><Type>2</Type></UPD>",
        "<UPD ID=\"TMPGdata\"><Val>0000A042</Val><Type>4</Type></UPD></TPLD>",
        "<TPLD><UPD ID=\"TMPGtime\"><Val>00002041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"TMPGcode\"><Val>6</Val><Type>2</Type></UPD>",
        "<UPD ID=\"TMPGdata\"><Val>00002041</Val><Type>4</Type></UPD></TPLD>",
        "</TP></GUD>"
    );
    let mth_u16 = encode_utf16le(method_xml);
    let signal = pack_trace_stream(600, &[&[0x05]]);

    let lcd = build_compound_file(&[
        ("File Property", prop),
        ("GUMM_Information/ShimadzuLC.1/GUC.1.METHOD", &mth_u16),
        ("PDA 3D Raw Data/Max Plot", &signal),
    ]);

    assert_eq!(
        sequency_instruments::detect(None, &lcd),
        VendorFormat::ShimadzuLcd
    );
    let reading = read_lcd(&lcd).unwrap();
    assert_eq!(reading.sample_name, None);
    assert_eq!(read_run(&lcd).unwrap().sample_name, None);
    assert!(reading.method.solvents.is_empty());
    assert_eq!(reading.method.flow_rate_ml_min, Some(0.5));
    assert_eq!(
        reading.method.gradient,
        vec![
            GradientStep {
                time_min: 0.0,
                b_pct: 10.0,
                c_pct: 5.0,
                d_pct: 0.0,
            },
            GradientStep {
                time_min: 10.0,
                b_pct: 80.0,
                c_pct: 10.0,
                d_pct: 0.0,
            },
        ]
    );
}

#[test]
fn ole_difat_sector_chain_supports_more_than_109_fat_sectors() {
    let sector = 512usize;
    let fat_sectors = 112usize;
    let dir_sec = 0u32;
    let data_sec = 1u32;
    let difat_sec = 2u32;
    let first_fat_sec = 3u32;
    let total_sectors = 3 + fat_sectors;

    let mut body = vec![0u8; total_sectors * sector];

    let payload = vec![0x5Au8; 4200];
    let payload_sectors = payload.len().div_ceil(sector);
    let extra_data_start = total_sectors;
    body.resize((extra_data_start + payload_sectors) * sector, 0);
    body[data_sec as usize * sector..data_sec as usize * sector + sector]
        .copy_from_slice(&payload[..sector]);
    for i in 1..payload_sectors {
        let dst = (extra_data_start + i - 1) * sector;
        let src_start = i * sector;
        let src_end = (src_start + sector).min(payload.len());
        body[dst..dst + (src_end - src_start)].copy_from_slice(&payload[src_start..src_end]);
    }

    let mut root = [0u8; 128];
    let mut root_name = encode_utf16le("Root Entry");
    root_name.extend_from_slice(&[0, 0]);
    root[..root_name.len()].copy_from_slice(&root_name);
    root[64..66].copy_from_slice(&(root_name.len() as u16).to_le_bytes());
    root[66] = 5;
    root[68..72].copy_from_slice(&FREE_SECT.to_le_bytes());
    root[72..76].copy_from_slice(&FREE_SECT.to_le_bytes());
    root[76..80].copy_from_slice(&1u32.to_le_bytes());
    root[116..120].copy_from_slice(&END_OF_CHAIN.to_le_bytes());

    let mut stream_entry = [0u8; 128];
    let mut s_name = encode_utf16le("BigStream");
    s_name.extend_from_slice(&[0, 0]);
    stream_entry[..s_name.len()].copy_from_slice(&s_name);
    stream_entry[64..66].copy_from_slice(&(s_name.len() as u16).to_le_bytes());
    stream_entry[66] = 2;
    stream_entry[68..72].copy_from_slice(&FREE_SECT.to_le_bytes());
    stream_entry[72..76].copy_from_slice(&FREE_SECT.to_le_bytes());
    stream_entry[76..80].copy_from_slice(&FREE_SECT.to_le_bytes());
    stream_entry[116..120].copy_from_slice(&data_sec.to_le_bytes());
    stream_entry[120..124].copy_from_slice(&(payload.len() as u32).to_le_bytes());

    body[0..128].copy_from_slice(&root);
    body[128..256].copy_from_slice(&stream_entry);
    for slot in 2..4 {
        let base = slot * 128;
        body[base + 68..base + 72].copy_from_slice(&FREE_SECT.to_le_bytes());
        body[base + 72..base + 76].copy_from_slice(&FREE_SECT.to_le_bytes());
        body[base + 76..base + 80].copy_from_slice(&FREE_SECT.to_le_bytes());
    }

    let difat_base = difat_sec as usize * sector;
    for i in 0..127usize {
        let val = if 109 + i < fat_sectors {
            first_fat_sec + (109 + i) as u32
        } else {
            FREE_SECT
        };
        body[difat_base + i * 4..difat_base + (i + 1) * 4].copy_from_slice(&val.to_le_bytes());
    }
    body[difat_base + 127 * 4..difat_base + 128 * 4].copy_from_slice(&END_OF_CHAIN.to_le_bytes());

    let mut fat = vec![FREE_SECT; fat_sectors * 128];
    fat[dir_sec as usize] = END_OF_CHAIN;
    fat[difat_sec as usize] = 0xFFFF_FFFC;
    for i in 0..fat_sectors {
        fat[first_fat_sec as usize + i] = FAT_SECT;
    }
    let mut chain = vec![data_sec];
    for i in 1..payload_sectors {
        chain.push((extra_data_start + i - 1) as u32);
    }
    link_chains(&[chain], &mut fat);

    for (i, &entry) in fat.iter().enumerate() {
        let off = (first_fat_sec as usize * sector) + i * 4;
        body[off..off + 4].copy_from_slice(&entry.to_le_bytes());
    }

    let mut header = vec![0u8; 512];
    header[..8].copy_from_slice(&OLE_MAGIC);
    header[24..26].copy_from_slice(&0x003Eu16.to_le_bytes());
    header[26..28].copy_from_slice(&3u16.to_le_bytes());
    header[28..30].copy_from_slice(&0xFFFEu16.to_le_bytes());
    header[30..32].copy_from_slice(&9u16.to_le_bytes());
    header[32..34].copy_from_slice(&6u16.to_le_bytes());
    header[44..48].copy_from_slice(&(fat_sectors as u32).to_le_bytes());
    header[48..52].copy_from_slice(&dir_sec.to_le_bytes());
    header[56..60].copy_from_slice(&4096u32.to_le_bytes());
    header[60..64].copy_from_slice(&END_OF_CHAIN.to_le_bytes());
    header[68..72].copy_from_slice(&difat_sec.to_le_bytes());
    header[72..76].copy_from_slice(&1u32.to_le_bytes());
    for i in 0..109usize {
        let off = 76 + i * 4;
        header[off..off + 4].copy_from_slice(&(first_fat_sec + i as u32).to_le_bytes());
    }
    header.extend_from_slice(&body);

    let ole = OleFile::parse(&header).unwrap();
    assert_eq!(
        ole.stream("BigStream").unwrap().as_ref(),
        payload.as_slice()
    );
}

#[test]
fn shimadzu_lcm_method_file_detects_and_parses_method_fields() {
    let prop =
        b"<?xml version=\"1.0\"?><FileProperty><dwFileType>16777248</dwFileType></FileProperty>";
    let config_xml = "<UnitList><Unit><ModelName>CBM-20A</ModelName></Unit></UnitList>";
    let method_xml = concat!(
        "<GUD Type=\"SystemParameter\">",
        "<UPD ID=\"PumpA_MobilePhaseName#1A\"><Val>0.05%磷酸水溶液</Val><Type>8</Type></UPD>",
        "<UPD ID=\"PumpA_MobilePhaseName#2A\"><Val>乙腈</Val><Type>8</Type></UPD>",
        "<UPD ID=\"PumpA_Flow\"><Val>00000000</Val><Type>4</Type></UPD>",
        "<UPD ID=\"LPGE$tflow\"><Val>9A99993F</Val><Type>4</Type></UPD>",
        "<UPD ID=\"OVEN$use\"><Val>1</Val><Type>2</Type></UPD>",
        "<UPD ID=\"OVEN$ovent\"><Val>00004842</Val><Type>4</Type></UPD>",
        "<UPD ID=\"LPGE$cconc\"><Val>00008041</Val><Type>4</Type></UPD>",
        "<UPD ID=\"LPGE$dconc\"><Val>0000A842</Val><Type>4</Type></UPD>",
        "<TP PgN=\"SCL$TIMEPROG\">",
        "<TPLD><TPD ID=\"TMPGtime\"><Val>0000A041</Val><Type>4</Type></TPD>",
        "<TPD ID=\"TMPGcode\"><Val>6</Val><Type>18</Type></TPD>",
        "<TPD ID=\"TMPGdata\"><Val>0000A042</Val><Type>4</Type></TPD></TPLD>",
        "<TPLD><TPD ID=\"TMPGtime\"><Val>0000A041</Val><Type>4</Type></TPD>",
        "<TPD ID=\"TMPGcode\"><Val>7</Val><Type>18</Type></TPD>",
        "<TPD ID=\"TMPGdata\"><Val>0000A041</Val><Type>4</Type></TPD></TPLD>",
        "<TPLD><TPD ID=\"TMPGtime\"><Val>0000A041</Val><Type>4</Type></TPD>",
        "<TPD ID=\"TMPGcode\"><Val>97</Val><Type>18</Type></TPD>",
        "<TPD ID=\"TMPGdata\"><Val>00000000</Val><Type>4</Type></TPD></TPLD>",
        "</TP>",
        "<TP PgN=\"PUMP$PROG\">",
        "<TPLD><TPD ID=\"Flow\"><Val>00002041</Val><Type>4</Type></TPD></TPLD>",
        "</TP></GUD>"
    );
    let pda_xml = "<GUD><UPD ID=\"2DCH_Wavelength#1\"><Val>21000</Val><Type>3</Type></UPD></GUD>";

    let cfg_u16 = encode_utf16le(config_xml);
    let mth_u16 = encode_utf16le(method_xml);
    let pda_u16 = encode_utf16le(pda_xml);

    let lcm = build_compound_file(&[
        ("File Property", prop),
        ("GUMM_Information/ShimadzuLC.1/GUM.CONFIG", &cfg_u16),
        ("GUMM_Information/ShimadzuLC.1/GUC.1.METHOD", &mth_u16),
        ("GUMM_Information/ShimadzuPDA.1/PDA.1.METHOD", &pda_u16),
    ]);

    assert_eq!(
        sequency_instruments::detect(Some("HPLC-95B.lcm"), &lcm),
        VendorFormat::ShimadzuLcm
    );
    assert_eq!(
        sequency_instruments::detect(None, &lcm),
        VendorFormat::ShimadzuLcm
    );
    assert_eq!(
        sequency_instruments::detect(Some("HPLC-95B.lcm"), &[]),
        VendorFormat::ShimadzuLcm
    );

    let reading = sequency_instruments::read(Some("E:/methods/HPLC-95B.lcm"), &lcm).unwrap();
    let direct = read_lcm(&lcm).unwrap();
    assert_eq!(reading.method.name.as_deref(), Some("HPLC-95B.lcm"));
    assert_eq!(direct.method.name, None);
    assert_eq!(reading.technique, Technique::Hplc);
    assert_eq!(reading.vendor_format, VendorFormat::ShimadzuLcm);
    assert_eq!(reading.instrument.as_deref(), Some("Shimadzu CBM-20A"));
    assert_eq!(reading.method.solvents, vec!["0.05%磷酸水溶液", "乙腈"]);
    assert_eq!(reading.method.flow_rate_ml_min, Some(1.2));
    assert_eq!(reading.method.column_temp_c, Some(50.0));
    assert_eq!(reading.method.wavelength_nm, Some(210.0));
    assert_eq!(
        reading.method.gradient,
        vec![
            GradientStep {
                time_min: 0.0,
                b_pct: 0.0,
                c_pct: 16.0,
                d_pct: 84.0,
            },
            GradientStep {
                time_min: 20.0,
                b_pct: 0.0,
                c_pct: 80.0,
                d_pct: 20.0,
            },
        ]
    );
    assert_eq!(reading.method.gradient[0].a_pct(), 0.0);
    assert_eq!(reading.method.gradient[1].a_pct(), 0.0);
    assert!(reading.peaks.is_empty());
    assert!(reading.traces.is_empty());

    let word_doc = build_compound_file(&[("WordDocument", b"not a shimadzu method")]);
    assert_eq!(
        sequency_instruments::detect(None, &word_doc),
        VendorFormat::Ole2
    );
    assert!(matches!(
        read_lcm(&word_doc),
        Err(InstrumentError::InvalidShimadzu(_))
    ));

    let run_prop = b"<?xml version=\"1.0\"?><SampleInfo><smpl_name>TBP-1</smpl_name></SampleInfo>";
    let gumm_lcd = build_compound_file(&[
        ("File Property", run_prop),
        ("GUMM_Information/ShimadzuLC.1/GUC.1.METHOD", &mth_u16),
    ]);
    assert_eq!(
        sequency_instruments::detect(None, &gumm_lcd),
        VendorFormat::ShimadzuLcd
    );

    let prop_only = build_compound_file(&[("File Property", run_prop)]);
    assert_eq!(
        sequency_instruments::detect(None, &prop_only),
        VendorFormat::ShimadzuLcd
    );
}

#[test]
fn shimadzu_read_run_includes_parsed_method_and_all_peaks() {
    let prop = format!(
        "<?xml version=\"1.0\"?><SampleInfo><smpl_name>{}</smpl_name><inj_vol>@FtoX@40400000</inj_vol></SampleInfo>",
        stox_gb18030(" TBP-230538 ")
    );
    let method_xml = concat!(
        "<GUD>",
        "<UPD ID=\"PumpA_Flow\"><Val>0000803F</Val><Type>4</Type></UPD>",
        "<UPD ID=\"OVEN$use\"><Val>1</Val><Type>2</Type></UPD>",
        "<UPD ID=\"OVEN$ovent\"><Val>0000F041</Val><Type>4</Type></UPD>",
        "</GUD>"
    );
    let mth_u16 = encode_utf16le(method_xml);
    let mut pt = b"VER1".to_vec();
    pt.extend_from_slice(&2i32.to_le_bytes());
    pt.extend_from_slice(&[0u8; 12]);
    pt.extend_from_slice(&pack_peak_record(300_000, 1000.0, 200.0, 64));
    pt.extend_from_slice(&pack_peak_record(600_000, 9000.0, 1800.0, 64));

    let lcd = build_compound_file(&[
        ("File Property", prop.as_bytes()),
        ("GUMM_Information/ShimadzuLC.1/GUC.1.METHOD", &mth_u16),
        ("LSS Data Processing/PT-PDA.1", &pt),
    ]);

    let run = read_run(&lcd).unwrap();
    assert_eq!(run.sample_name.as_deref(), Some("TBP-230538"));
    assert_eq!(run.peaks.len(), 2);
    assert_eq!(run.method.flow_rate_ml_min, Some(1.0));
    assert_eq!(run.method.column_temp_c, Some(30.0));
    assert_eq!(run.method.injection_volume_ul, Some(3.0));
}
