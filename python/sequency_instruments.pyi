from collections.abc import Sequence

class GradientStep:
    time_min: float
    b_pct: float
    c_pct: float
    d_pct: float
    def a_pct(self) -> float: ...

class InstrumentMethod:
    name: str | None
    column: str | None
    solvents: list[str]
    gradient: list[GradientStep]
    flow_rate_ml_min: float | None
    column_temp_c: float | None
    injection_volume_ul: float | None
    wavelength_nm: float | None
    nucleus: str | None
    frequency_mhz: float | None
    pulse_program: str | None
    scans: int | None
    relaxation_delay_s: float | None

class Peak:
    rt_min: float | None
    area: float | None
    area_pct: float | None
    height: float | None
    compound: str | None

class Trace:
    label: str
    x_unit: str
    y_unit: str
    points: list[tuple[float, float]]
    def downsample_lttb(self, threshold: int) -> Trace: ...
    def downsample_minmax(self, bins: int) -> Trace: ...

class InstrumentReading:
    technique: str
    vendor_format: str
    acquired_at: str | None
    sample_name: str | None
    instrument: str | None
    method: InstrumentMethod
    peaks: list[Peak]
    traces: list[Trace]
    def main_peak(self) -> Peak | None: ...

class ShimadzuRun:
    injected: str | None
    sample_name: str | None
    data_file: str | None
    method_file: str | None
    instrument: str
    method: InstrumentMethod
    peaks: list[Peak]

class OleStream:
    name: str
    data: bytes
    extents: list[tuple[int, int]]

class OleLayout:
    streams: list[OleStream]
    unallocated: list[tuple[int, int]]

class AcqusParams:
    acquired_at: str | None
    nucleus: str
    pulse_program: str
    instrument: str
    solvent: str
    scans: int
    relaxation_delay_s: float
    sfo1_mhz: float | None
    bf1_mhz: float | None
    sw_h_hz: float | None
    data_path: str | None

class OpenLabInjection:
    id: str | None
    data_file: str
    sample_name: str
    sample_type: str
    vial: str
    instrument: str | None
    acquired_utc: str | None
    method: str
    processing_method: str | None
    status: str | None
    result_file: str | None
    injection_volume_ul: float | None

class AcamlResultSet:
    injections: list[OpenLabInjection]
    logins: list[str]

def detect(data: bytes, filename: str | None = None) -> str: ...
def read(data: bytes, filename: str | None = None) -> InstrumentReading: ...
def read_files(files: Sequence[tuple[str, bytes]]) -> InstrumentReading: ...
def ole_layout(data: bytes) -> OleLayout: ...
def ole_place(extents: Sequence[tuple[int, int]], start: int, end: int) -> list[tuple[int, int]]: ...
def shimadzu_read(data: bytes) -> InstrumentReading: ...
def shimadzu_read_run(data: bytes) -> ShimadzuRun: ...
def shimadzu_read_lcm(data: bytes) -> InstrumentReading: ...
def shimadzu_properties(stream: bytes) -> dict[str, str]: ...
def shimadzu_filetime(low: str | None = None, high: str | None = None) -> str | None: ...
def shimadzu_peak_table(stream: bytes) -> list[Peak]: ...
def shimadzu_deltas(sub_block: bytes) -> list[int]: ...
def shimadzu_trace(stream: bytes) -> tuple[list[float], list[int]]: ...
def shimadzu_chromatograms(data: bytes) -> list[Trace]: ...
def bruker_parse_acqus(text: str) -> AcqusParams: ...
def bruker_read_acqus(data: bytes) -> InstrumentReading: ...
def bruker_read_pdata(procs: bytes, data_1r: bytes, acqus: bytes | None = None) -> InstrumentReading: ...
def bruker_read_dir(exp_dir: str) -> InstrumentReading: ...
def openlab_parse_acaml(acaml: bytes) -> AcamlResultSet: ...
def openlab_read_acaml(acaml: bytes) -> InstrumentReading: ...
def openlab_read_dx(dx: bytes) -> InstrumentReading: ...
def openlab_read_rx(rx: bytes) -> InstrumentReading: ...
def openlab_parse_rx_peaks(rx: bytes) -> list[Peak]: ...

class AndiDataset:
    sample_name: str | None
    operator_name: str | None
    source_file_reference: str | None
    detection_method_name: str | None
    injection_date_time_stamp: str | None
    acquired_at: str | None
    sampling_interval_s: float | None
    peaks: list[Peak]
    traces: list[Trace]

class OmnicHeader:
    title: str
    y_unit: str | None
    acquired_epoch_secs: int

class OmnicSpectrum:
    wavenumbers: list[float]
    intensities: list[float]
    title: str
    y_unit: str | None
    acquired_epoch_secs: int

class UvTextTable:
    wavelengths: list[float]
    scans: list[list[float]]

class ThermoRow:
    sample_name: str
    sample_id: str
    data_file: str
    vial: str
    method: str

class ThermoFile:
    first: str
    second: str
    logins: list[str]
    rows: list[ThermoRow]
    instrument: str | None
    scans: int
    minutes: tuple[float, float] | None

class ThermoScans:
    minutes: list[float]
    tic: list[float]
    offsets: list[int]
    data_offset: int

class WatersUvChannel:
    wavelength: float
    values: list[float]

class WatersPolarityRun:
    polarity: str
    acquired: str | None
    instrument: str | None
    inlet_method: str | None
    ms_method: str | None
    rt: list[float]
    mz: list[float]
    intensity: list[float]
    uv_rt: list[float] | None
    uv: list[WatersUvChannel]

class MassHunterScan:
    minutes: float
    polarity: str
    tic: float
    offset: int
    length: int
    points: int
    calibration: list[float]
    flags: int | None

class MnovaEmbedded:
    name: str
    entry: tuple[int, int]
    start: int
    end: int
    path: tuple[int, int]
    after: tuple[int, int]
    text: bytes

class MnovaRow:
    span: tuple[int, int]
    label: str | None
    key: str | None
    formula: str | None
    value: tuple[int, int]

class MnovaSpectrum:
    dataset: str
    expno: str
    position: int
    acqus: AcqusParams

class MnovaJ2kFrame:
    at: int
    position: int | None
    nucleus: str
    low_hz: float
    mhz: float
    width_hz: float
    low_val: float
    high_val: float
    depth: int
    codestream: bytes

def andi_parse_cdf(data: bytes) -> AndiDataset: ...
def chromeleon_signal_name(data: bytes) -> str | None: ...
def chromeleon_decode_signal(data: bytes) -> tuple[list[float], list[float]]: ...
def omnic_parse_header(data: bytes) -> OmnicHeader: ...
def omnic_parse_spa(data: bytes) -> OmnicSpectrum: ...
def uvtext_parse(data: bytes) -> UvTextTable: ...
def thermo_parse_run(data: bytes) -> ThermoFile: ...
def thermo_parse_sequence(data: bytes) -> ThermoFile: ...
def thermo_ms_header_offset(data: bytes) -> int | None: ...
def thermo_decode_scan_index(entries: bytes, data_offset: int) -> ThermoScans: ...
def thermo_decode_peak_record(record: bytes) -> tuple[int, int]: ...
def thermo_decode_peak_pairs(bytes: bytes) -> tuple[list[float], list[float]]: ...
def waters_read_polarities(files: list[tuple[str, bytes]], bin_width: float = 0.1) -> dict[str, WatersPolarityRun]: ...
def masshunter_profile_scans(
    msscan_bin: bytes,
    msscan_xsd: bytes,
    msts_xml: bytes | None = None,
    masscal_bin: bytes | None = None,
    default_masscal_xml: bytes | None = None,
) -> list[MassHunterScan]: ...
def masshunter_decode_spectrum(
    segment: bytes,
    num_points: int,
    calibration: Sequence[float],
    use_flags: int | None = None,
) -> tuple[list[float], list[float]]: ...
def mnova_embedded(content: bytes) -> list[MnovaEmbedded]: ...
def mnova_rows(content: bytes) -> list[MnovaRow]: ...
def mnova_measurements(content: bytes) -> list[tuple[int, int]]: ...
def mnova_spectra(content: bytes) -> list[MnovaSpectrum]: ...
def mnova_j2k_frames(content: bytes) -> list[MnovaJ2kFrame]: ...
