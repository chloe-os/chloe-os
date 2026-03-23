"""Shared data models for the chloe-heart cardiac monitoring pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


class ArrhythmiaType(StrEnum):
    """Classification of arrhythmia events."""

    NORMAL = "normal_sinus"
    SINUS_ARRHYTHMIA = "sinus_arrhythmia"  # Normal in dogs
    PVC = "premature_ventricular_complex"
    PAC = "premature_atrial_complex"
    AFIB = "atrial_fibrillation"
    VT = "ventricular_tachycardia"
    BRADYCARDIA = "bradycardia"
    TACHYCARDIA = "tachycardia"
    SECOND_DEGREE_AVB = "second_degree_av_block"
    OTHER = "other"


class RiskLevel(StrEnum):
    """Cardiac risk classification."""

    GREEN = "green"
    YELLOW = "yellow"
    RED = "red"


class DogSize(StrEnum):
    """Dog size category for breed-specific baselines."""

    GIANT = "giant"
    LARGE = "large"
    MEDIUM = "medium"
    SMALL = "small"
    TOY = "toy"


@dataclass
class BreedBaseline:
    """Breed-specific cardiac reference ranges."""

    size: DogSize
    hr_min: int  # Normal resting HR lower bound (BPM)
    hr_max: int  # Normal resting HR upper bound (BPM)
    sdnn_min: float  # Normal SDNN lower bound (ms)
    sdnn_max: float  # Normal SDNN upper bound (ms)
    dcm_risk: str  # "low", "moderate", "high"
    notes: str = ""


# Breed-specific cardiac baselines
BREED_BASELINES: dict[str, BreedBaseline] = {
    "great_dane": BreedBaseline(DogSize.GIANT, 60, 100, 80, 200, "high", "DCM with AFib common"),
    "irish_wolfhound": BreedBaseline(DogSize.GIANT, 60, 100, 80, 200, "high", "DCM + AFib"),
    "newfoundland": BreedBaseline(DogSize.GIANT, 60, 100, 80, 200, "moderate", ""),
    "doberman": BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "high", "VPCs are early DCM sign"),
    "boxer": BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "high", "ARVC common"),
    "labrador": BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "moderate", ""),
    "golden_retriever": BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "moderate", ""),
    "german_shepherd": BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "low", ""),
    "staffy": BreedBaseline(DogSize.MEDIUM, 70, 130, 50, 130, "low", ""),
    "beagle": BreedBaseline(DogSize.MEDIUM, 70, 130, 50, 130, "low", ""),
    "cocker_spaniel": BreedBaseline(DogSize.MEDIUM, 70, 130, 50, 130, "low", "DCM occasionally"),
    "bulldog": BreedBaseline(DogSize.MEDIUM, 70, 130, 50, 130, "low", ""),
    "cavalier": BreedBaseline(DogSize.SMALL, 80, 140, 40, 100, "low", "MMVD near-universal by 10"),
    "miniature_poodle": BreedBaseline(DogSize.SMALL, 80, 140, 40, 100, "low", ""),
    "dachshund": BreedBaseline(DogSize.SMALL, 80, 140, 40, 100, "low", "MMVD risk"),
    "chihuahua": BreedBaseline(DogSize.TOY, 90, 160, 30, 80, "low", ""),
    "yorkshire_terrier": BreedBaseline(DogSize.TOY, 90, 160, 30, 80, "low", ""),
    "pomeranian": BreedBaseline(DogSize.TOY, 90, 160, 30, 80, "low", ""),
}

# Default baselines by size category
SIZE_BASELINES: dict[DogSize, BreedBaseline] = {
    DogSize.GIANT: BreedBaseline(DogSize.GIANT, 60, 100, 80, 200, "moderate"),
    DogSize.LARGE: BreedBaseline(DogSize.LARGE, 60, 120, 60, 150, "moderate"),
    DogSize.MEDIUM: BreedBaseline(DogSize.MEDIUM, 70, 130, 50, 130, "low"),
    DogSize.SMALL: BreedBaseline(DogSize.SMALL, 80, 140, 40, 100, "low"),
    DogSize.TOY: BreedBaseline(DogSize.TOY, 90, 160, 30, 80, "low"),
}


def get_baseline(breed: str | None = None, size: DogSize | None = None) -> BreedBaseline:
    """Look up the cardiac baseline for a breed or size category."""
    if breed:
        key = breed.lower().replace(" ", "_").replace("-", "_")
        if key in BREED_BASELINES:
            return BREED_BASELINES[key]
    if size:
        return SIZE_BASELINES[size]
    return SIZE_BASELINES[DogSize.MEDIUM]


@dataclass
class ECGSignal:
    """Raw ECG signal data from Stage 1 ingestion."""

    signal: NDArray[np.float64]  # Voltage values
    sample_rate: float  # Samples per second (Hz)
    duration_seconds: float
    source_file: str
    source_format: str  # "csv", "edf", "wfdb"
    channel_name: str = "ECG"
    units: str = "mV"

    @property
    def num_samples(self) -> int:
        return len(self.signal)


@dataclass
class HRTimeSeries:
    """Heart rate / HRV time series from a wearable device."""

    timestamps: NDArray[np.float64]  # Unix timestamps or seconds from start
    heart_rates: NDArray[np.float64]  # BPM values
    rr_intervals: NDArray[np.float64] | None = None  # R-R intervals in ms (if available)
    source_file: str = ""
    device_name: str = ""


@dataclass
class QRSComplex:
    """A detected QRS complex (heartbeat) in the ECG."""

    r_peak_index: int  # Sample index of R-peak
    r_peak_time: float  # Time in seconds
    r_peak_amplitude: float  # Amplitude in mV
    qrs_onset: int | None = None  # Sample index
    qrs_offset: int | None = None  # Sample index
    qrs_duration_ms: float | None = None


@dataclass
class ProcessedECG:
    """Output of Stage 2: preprocessed ECG with detected beats."""

    filtered_signal: NDArray[np.float64]
    sample_rate: float
    duration_seconds: float
    qrs_complexes: list[QRSComplex]
    rr_intervals_ms: NDArray[np.float64]  # R-R intervals in milliseconds

    @property
    def num_beats(self) -> int:
        return len(self.qrs_complexes)

    @property
    def mean_hr_bpm(self) -> float:
        if len(self.rr_intervals_ms) == 0:
            return 0.0
        mean_rr_s = float(np.mean(self.rr_intervals_ms)) / 1000.0
        return 60.0 / mean_rr_s if mean_rr_s > 0 else 0.0


@dataclass
class ArrhythmiaEvent:
    """A detected arrhythmia event."""

    arrhythmia_type: ArrhythmiaType
    time_seconds: float
    duration_seconds: float | None = None
    beat_index: int | None = None
    confidence: float = 1.0  # 0-1
    details: str = ""


@dataclass
class HRVMetrics:
    """Heart rate variability metrics."""

    # Time-domain
    sdnn: float  # Standard deviation of R-R intervals (ms)
    rmssd: float  # Root mean square of successive differences (ms)
    pnn50: float  # % of intervals differing by > 50ms
    mean_rr: float  # Mean R-R interval (ms)
    mean_hr: float  # Mean heart rate (BPM)

    # Frequency-domain (optional)
    lf_power: float | None = None  # Low-frequency power (ms²)
    hf_power: float | None = None  # High-frequency power (ms²)
    lf_hf_ratio: float | None = None  # Sympathetic/parasympathetic balance
    total_power: float | None = None  # Total spectral power (ms²)


@dataclass
class CardiacAnalysis:
    """Output of Stage 3: complete cardiac analysis results."""

    arrhythmia_events: list[ArrhythmiaEvent]
    hrv_metrics: HRVMetrics
    mean_hr_bpm: float
    min_hr_bpm: float
    max_hr_bpm: float
    total_beats: int
    duration_seconds: float
    abnormal_beat_count: int

    @property
    def arrhythmia_burden_pct(self) -> float:
        """Percentage of beats that are abnormal."""
        if self.total_beats == 0:
            return 0.0
        return (self.abnormal_beat_count / self.total_beats) * 100.0

    @property
    def has_vtach(self) -> bool:
        return any(e.arrhythmia_type == ArrhythmiaType.VT for e in self.arrhythmia_events)

    @property
    def has_afib(self) -> bool:
        return any(e.arrhythmia_type == ArrhythmiaType.AFIB for e in self.arrhythmia_events)


@dataclass
class CardiacHealthScore:
    """Output of Stage 4: composite cardiac health score."""

    overall_score: float  # 0-100
    risk_level: RiskLevel

    # Component scores (0-1 each)
    rhythm_score: float
    hrv_score: float
    heart_rate_score: float
    trend_score: float
    exercise_score: float

    # Risk flags
    risk_flags: list[str] = field(default_factory=list)

    # Weights used
    weights: dict[str, float] = field(default_factory=dict)


@dataclass
class HeartConfig:
    """Configuration for a chloe-heart analysis run."""

    breed: str | None = None
    size: DogSize | None = None
    ecg_source: str | None = None  # Path to ECG file
    hr_source: str | None = None  # Path to HR/HRV data
    output_path: str = "cardiac_report.html"
    include_ai_interpretation: bool = False
    llm_provider: str = "anthropic"
    llm_model: str | None = None

    @property
    def baseline(self) -> BreedBaseline:
        return get_baseline(self.breed, self.size)
