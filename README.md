# sequency-instruments

This is a mirror of a folder in our internal codebase. While issues and PRs will be monitored, we won't accept direct contributions to this repository.

`sequency-instruments` parses binary analytical chemistry instrument files in pure Rust, with zero C dependencies and zero unsafe code. Python bindings are included in the `python` package and published to PyPI as `sequency-instruments`.

## Supported formats

The library extracts unified metadata, method settings, peak tables, and traces from eight file families:

- parses Shimadzu LabSolutions `.lcd` run files and `.lcm` method files stored in the OLE2 Compound File Binary Format, including `File Property` fields, `GUC.1.METHOD` mobile phases, flow, oven temperature, and gradients, `PDA.1.METHOD` wavelengths, `PT-*` peak tables, and delta-encoded chromatograms
- reads Bruker TopSpin NMR runs, including `acqus` and `procs` JCAMP-DX parameters, `BYTORDP` byte-order handling, `NC_proc` power-of-two scaling, and `1r` real spectra
- parses Agilent OpenLab CDS `.acaml` XML result sets and zipped `.dx` and `.rx` archives, including `Base/InjectionACAML` peaks and ChemStation `.CH` channels
- parses standalone Agilent ChemStation `.D` folders, `.ch` chromatograms, and `.ms` (`MSD1.MS`) total ion current traces
- parses ANDI Chromatography netCDF-3 `.cdf` files (CDF-1 and CDF-2)
- decodes Thermo Chromeleon 7 `SignHdr` `.raw` signal streams (`PtsLL2Df` blocks)
- reads Waters MassLynx `.raw` LC-MS directories (`_HEADER.TXT`, `_extern.inf`, `_FUNC*.IDX`, and 2-, 4-, 6-, and 8-byte `_FUNC*.DAT` traces)
- reads Newstyle (Qingbohua) NU3000 prep-HPLC XML datasets (`.Data`, `.FilterData`, `.icm` peaks, and `.sys` device metadata)

`detect` identifies file formats from magic bytes and filenames. `Trace::downsample_lttb` and `Trace::downsample_minmax` reduce dense curves for plotting while preserving narrow peaks.

## Rust usage

Add `sequency-instruments` to `Cargo.toml`:

```toml
[dependencies]
sequency-instruments = "0.1"
```

Read a single file with automatic format detection:

```rust
use std::fs;

let bytes = fs::read("sample.lcd")?;
let reading = sequency_instruments::read(Some("sample.lcd"), &bytes)?;
println!("technique={} peaks={} traces={}", reading.technique, reading.peaks.len(), reading.traces.len());
```

Read a multi-file bundle such as a Bruker NMR directory or an OpenLab result set:

```rust
let reading = sequency_instruments::bruker::read_dir(std::path::Path::new("10"))?;
for trace in &reading.traces {
    let preview = trace.downsample_minmax(400);
    println!("{}: {} points", preview.label, preview.points.len());
}
```

## Python usage

Install the prebuilt wheel from PyPI:

```sh
pip install sequency-instruments
```

Parse files from Python:

```python
from pathlib import Path
import sequency_instruments as si

data = Path("sample.lcd").read_bytes()
print(si.detect(data, "sample.lcd"))

reading = si.read(data, "sample.lcd")
print(reading.technique, reading.sample_name, len(reading.peaks))

for trace in reading.traces:
    preview = trace.downsample_minmax(400)
    print(trace.label, len(preview.points))
```

Format-specific helpers are also exported for OLE2 containers (`ole_layout`, `ole_place`), Shimadzu `.lcd` and `.lcm` files (`shimadzu_read`, `shimadzu_read_run`, `shimadzu_read_lcm`, `shimadzu_properties`, `shimadzu_filetime`, `shimadzu_peak_table`, `shimadzu_deltas`, `shimadzu_trace`, `shimadzu_chromatograms`), Bruker NMR files (`bruker_parse_acqus`, `bruker_read_acqus`, `bruker_read_pdata`, `bruker_read_dir`), and OpenLab archives (`openlab_parse_acaml`, `openlab_read_acaml`, `openlab_read_dx`, `openlab_read_rx`, `openlab_parse_rx_peaks`).

## Development

Run the Rust checks from the repository root:

```sh
cargo test
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
```

Check or build the Python extension from the repository root:

```sh
cargo check --manifest-path python/Cargo.toml
maturin build --manifest-path python/Cargo.toml
```

## License

Dual-licensed under the [MIT License](LICENSE-MIT) and the [Apache License, Version 2.0](LICENSE-APACHE).
