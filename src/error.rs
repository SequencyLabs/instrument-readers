//! Typed error definitions for instrument file readers.

use thiserror::Error;

/// Errors returned while detecting or parsing binary and text instrument files.
#[derive(Debug, Error)]
pub enum InstrumentError {
    /// The input is not a valid OLE2 Compound File Binary Format container.
    #[error("invalid OLE2 compound file: {0}")]
    InvalidOle(&'static str),

    /// A required stream is missing from an OLE2 compound file.
    #[error("missing OLE2 stream: {0}")]
    MissingStream(String),

    /// A FAT, MiniFAT, DIFAT, or directory chain leaves table bounds or loops.
    #[error("corrupt OLE2 chain at index {0}")]
    CorruptChain(u32),

    /// A Shimadzu LabSolutions `.lcd` stream could not be decoded.
    #[error("invalid Shimadzu .lcd data: {0}")]
    InvalidShimadzu(String),

    /// A Bruker TopSpin parameter or `1r` spectrum file is malformed.
    #[error("invalid Bruker NMR data: {0}")]
    InvalidBruker(String),

    /// An Agilent OpenLab `.acaml`, `.dx`, or `.rx` file is malformed.
    #[error("invalid OpenLab data: {0}")]
    InvalidOpenLab(String),

    /// An ANDI Chromatography netCDF-3 `.cdf` file is truncated or malformed.
    #[error("invalid ANDI netCDF data: {0}")]
    InvalidAndi(String),

    /// A Newstyle (Qingbohua) NU3000 XML chromatogram or peak file is malformed.
    #[error("invalid Newstyle NU3000 data: {0}")]
    InvalidNewstyle(String),

    /// An Agilent ChemStation `.ch` or `.ms` file is truncated or malformed.
    #[error("invalid ChemStation data: {0}")]
    InvalidChemStation(String),

    /// A Thermo Chromeleon `.raw` signal file is truncated or uses an unsupported block type.
    #[error("invalid Chromeleon signal data: {0}")]
    InvalidChromeleon(String),

    /// A Waters MassLynx `.raw` dataset is truncated or malformed.
    #[error("invalid Waters MassLynx data: {0}")]
    InvalidWaters(String),

    /// A Thermo Nicolet OMNIC `.spa` spectrum is truncated or malformed.
    #[error("invalid OMNIC .spa data: {0}")]
    InvalidOmnic(String),

    /// A plain-text UV export is malformed or contains no numeric points.
    #[error("invalid UV text data: {0}")]
    InvalidUvText(String),

    /// A Thermo Finnigan `.raw` or `.sld` file is truncated or uses an unsupported version.
    #[error("invalid Thermo Finnigan data: {0}")]
    InvalidThermo(String),

    /// An Agilent MassHunter `.d` scan or calibration file is missing or malformed.
    #[error("invalid Agilent MassHunter data: {0}")]
    InvalidMassHunter(String),

    /// A MestReNova `.mnova` document is truncated or malformed.
    #[error("invalid MestReNova data: {0}")]
    InvalidMnova(String),

    /// A ZIP container (`.dx` or `.rx`) is truncated, corrupt, or uses unsupported compression.
    #[error("invalid ZIP container: {0}")]
    InvalidZip(&'static str),

    /// XML parsing failed.
    #[error("XML parse error: {0}")]
    Xml(String),

    /// The input format is unrecognized or has no binary reader in this crate.
    #[error("unsupported instrument format: {0}")]
    UnsupportedFormat(String),

    /// Filesystem I/O failed.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Alias for results returned by `sequency-instruments`.
pub type Result<T> = std::result::Result<T, InstrumentError>;
