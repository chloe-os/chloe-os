"""Tests for chloe-heart data models."""

import numpy as np
from chloe_heart.models import (
    ArrhythmiaEvent,
    ArrhythmiaType,
    CardiacAnalysis,
    DogSize,
    HRVMetrics,
    ProcessedECG,
    QRSComplex,
    RiskLevel,
    get_baseline,
)


def test_get_baseline_by_breed():
    baseline = get_baseline(breed="labrador")
    assert baseline.size == DogSize.LARGE
    assert baseline.hr_min == 60
    assert baseline.hr_max == 120


def test_get_baseline_by_size():
    baseline = get_baseline(size=DogSize.TOY)
    assert baseline.hr_min == 90
    assert baseline.hr_max == 160


def test_get_baseline_default():
    baseline = get_baseline()
    assert baseline.size == DogSize.MEDIUM


def test_get_baseline_breed_normalization():
    b1 = get_baseline(breed="golden_retriever")
    b2 = get_baseline(breed="Golden Retriever")
    assert b1.hr_min == b2.hr_min


def test_processed_ecg_mean_hr():
    qrs = [
        QRSComplex(r_peak_index=i * 500, r_peak_time=i * 0.667, r_peak_amplitude=1.0)
        for i in range(10)
    ]
    rr = np.array([667.0] * 9)  # 667ms = 90 BPM
    processed = ProcessedECG(
        filtered_signal=np.zeros(5000),
        sample_rate=500.0,
        duration_seconds=10.0,
        qrs_complexes=qrs,
        rr_intervals_ms=rr,
    )
    assert processed.num_beats == 10
    assert abs(processed.mean_hr_bpm - 90.0) < 1.0


def test_cardiac_analysis_arrhythmia_burden():
    hrv = HRVMetrics(sdnn=80.0, rmssd=40.0, pnn50=15.0, mean_rr=667.0, mean_hr=90.0)
    analysis = CardiacAnalysis(
        arrhythmia_events=[],
        hrv_metrics=hrv,
        mean_hr_bpm=90.0,
        min_hr_bpm=75.0,
        max_hr_bpm=110.0,
        total_beats=100,
        duration_seconds=60.0,
        abnormal_beat_count=5,
    )
    assert abs(analysis.arrhythmia_burden_pct - 5.0) < 0.01


def test_cardiac_analysis_vtach_detection():
    hrv = HRVMetrics(sdnn=80.0, rmssd=40.0, pnn50=15.0, mean_rr=667.0, mean_hr=90.0)
    events = [
        ArrhythmiaEvent(arrhythmia_type=ArrhythmiaType.VT, time_seconds=10.0),
    ]
    analysis = CardiacAnalysis(
        arrhythmia_events=events,
        hrv_metrics=hrv,
        mean_hr_bpm=90.0,
        min_hr_bpm=75.0,
        max_hr_bpm=110.0,
        total_beats=100,
        duration_seconds=60.0,
        abnormal_beat_count=10,
    )
    assert analysis.has_vtach is True
    assert analysis.has_afib is False


def test_risk_level_enum():
    assert RiskLevel.GREEN == "green"
    assert RiskLevel.YELLOW == "yellow"
    assert RiskLevel.RED == "red"


def test_breed_specific_notes():
    doberman = get_baseline(breed="doberman")
    assert "VPC" in doberman.notes or "DCM" in doberman.notes
    cavalier = get_baseline(breed="cavalier")
    assert "MMVD" in cavalier.notes
