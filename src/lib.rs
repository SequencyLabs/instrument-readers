//! Pure-Rust binary instrument file readers for analytical chemistry formats.
//!
//! Provides unified metadata, method/gradient, peak-table, and trace extraction for:
//! - Shimadzu LabSolutions `.lcd` (OLE2 Compound File Binary Format)
//! - Bruker TopSpin NMR (`acqus`, `procs`, and `1r` binary real spectra)
//! - Agilent OpenLab CDS (`.acaml` XML metadata and zipped `.dx`/`.rx` containers)
//! - Automatic format detection by filename and magic bytes, plus LTTB and min-max bin downsampling.

#![cfg_attr(not(test), deny(clippy::panic))]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

mod binary;
pub mod bruker;
mod chemstation;
pub mod detect;
pub mod error;
pub mod ole;
pub mod openlab;
pub mod shimadzu;
pub mod types;
pub mod zip_archive;

pub use detect::{detect, read, read_files};
pub use error::{InstrumentError, Result};
pub use types::{
    GradientStep, InstrumentMethod, InstrumentReading, Peak, Technique, Trace, VendorFormat, lttb,
    minmax_bins,
};
