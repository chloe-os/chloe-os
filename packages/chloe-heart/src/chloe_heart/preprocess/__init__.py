"""Stage 2 — ECG signal preprocessing.

Applies filtering, QRS detection, and R-R interval computation to produce
a ``ProcessedECG`` ready for cardiac analysis.
"""

from __future__ import annotations

import numpy as np

from chloe_heart.models import ECGSignal, ProcessedECG
from chloe_heart.preprocess.filters import bandpass_filter, notch_filter
from chloe_heart.preprocess.qrs import compute_rr_intervals, detect_qrs

__all__ = [
    "bandpass_filter",
    "compute_rr_intervals",
    "detect_qrs",
    "notch_filter",
    "preprocess_ecg",
]


def preprocess_ecg(ecg: ECGSignal) -> ProcessedECG:
    """Run the full Stage 2 preprocessing pipeline on a raw ECG signal.

    Steps
    -----
    1. Remove powerline interference with a 50 Hz notch filter.
    2. Apply a 0.5 -- 40 Hz Butterworth bandpass filter to remove baseline
       wander and high-frequency noise.
    3. Detect QRS complexes using the Pan-Tompkins algorithm.
    4. Compute R-R intervals from the detected beats.

    Parameters
    ----------
    ecg : ECGSignal
        Raw ECG signal from Stage 1 ingestion.

    Returns
    -------
    ProcessedECG
        Filtered signal, detected QRS complexes, and R-R intervals.
    """
    signal = ecg.signal.astype(np.float64)
    sr = ecg.sample_rate

    # Step 1: Notch filter for powerline noise (50 Hz)
    filtered = notch_filter(signal, sr, freq=50.0)

    # Step 2: Bandpass filter (0.5 - 40 Hz)
    filtered = bandpass_filter(filtered, sr, low=0.5, high=40.0)

    # Step 3: QRS detection
    qrs_complexes = detect_qrs(filtered, sr)

    # Step 4: R-R intervals
    rr_intervals = compute_rr_intervals(qrs_complexes, sr)

    return ProcessedECG(
        filtered_signal=filtered,
        sample_rate=sr,
        duration_seconds=ecg.duration_seconds,
        qrs_complexes=qrs_complexes,
        rr_intervals_ms=rr_intervals,
    )
