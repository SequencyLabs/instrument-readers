//! Pure-Rust binary instrument file readers for analytical chemistry formats.
//!
//! Provides unified metadata, method/gradient, peak-table, and trace extraction for:
//! - Shimadzu LabSolutions `.lcd` run files and `.lcm` method files (OLE2 Compound File Binary Format)
//! - Bruker TopSpin NMR (`acqus`, `procs`, and `1r` binary real spectra)
//! - Agilent OpenLab CDS (`.acaml` XML metadata and zipped `.dx`/`.rx` containers)
//! - Agilent ChemStation `.D` folders, standalone `.ch` chromatograms, and `.ms` mass spectra
//! - ANDI Chromatography netCDF-3 (`.cdf` CDF-1 and CDF-2 chromatograms and peak tables)
//! - Thermo Chromeleon 7 `SignHdr` `.raw` detector signal streams
//! - Waters MassLynx `.raw` LC-MS directories (`_HEADER.TXT`, `_extern.inf`, `_FUNC*.IDX`/`.DAT`)
//! - Newstyle (Qingbohua) NU3000 prep-HPLC XML datasets (`.Data`, `.FilterData`, `.icm`, `.sys`)
//! - Automatic format detection by filename and magic bytes, plus LTTB and min-max bin downsampling.

#![cfg_attr(not(test), deny(clippy::panic))]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

pub mod andi;
mod binary;
pub mod bruker;
pub mod chemstation;
pub mod chromeleon;
pub mod detect;
pub mod error;
pub mod masshunter;
pub mod mnova;
pub mod netcdf;
pub mod newstyle;
pub mod ole;
pub mod omnic;
pub mod openlab;
pub mod shimadzu;
pub mod thermo;
pub mod types;
pub mod uvtext;
pub mod waters;
mod xml;
pub mod zip_archive;

pub use detect::{detect, read, read_files};
pub use error::{InstrumentError, Result};
pub use types::{
    GradientStep, InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat, lttb,
    minmax_bins,
};
