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
