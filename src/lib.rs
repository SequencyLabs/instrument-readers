//! Pure-Rust binary instrument file readers for analytical chemistry formats.
//!
//! Provides unified metadata, method/gradient, peak-table, and trace extraction for:
//! - Shimadzu LabSolutions `.lcd` run files and `.lcm` method files (OLE2 Compound File Binary Format)
//! - Bruker TopSpin NMR (`acqus`, `procs`, and `1r` binary real spectra)
//! - Agilent OpenLab CDS (`.acaml` XML metadata and zipped `.dx`/`.rx` containers)
//! - ANDI Chromatography netCDF-3 (`.cdf` CDF-1 and CDF-2 chromatograms and peak tables)
//! - Automatic format detection by filename and magic bytes, plus LTTB and min-max bin downsampling.

#![cfg_attr(not(test), deny(clippy::panic))]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

pub mod andi;
mod binary;
pub mod bruker;
mod chemstation;
pub mod detect;
pub mod error;
pub mod netcdf;
pub mod ole;
pub mod openlab;
pub mod shimadzu;
pub mod types;
mod xml;
pub mod zip_archive;

pub use detect::{detect, read, read_files};
pub use error::{InstrumentError, Result};
pub use types::{
    GradientStep, InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat, lttb,
    minmax_bins,
};
