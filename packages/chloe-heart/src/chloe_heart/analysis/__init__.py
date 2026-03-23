"""Stage 3 — Cardiac analysis: arrhythmia detection and HRV computation.

The public entry point is :func:`analyze_cardiac`, which orchestrates
arrhythmia classification and HRV metric computation and assembles
a :class:`~chloe_heart.models.CardiacAnalysis` result.
"""

from __future__ import annotations

import numpy as np

from chloe_heart.analysis.arrhythmia import classify_arrhythmias
from chloe_heart.analysis.hrv import compute_hrv
from chloe_heart.models import (
    ArrhythmiaType,
    BreedBaseline,
    CardiacAnalysis,
    ProcessedECG,
)

__all__ = ["analyze_cardiac", "classify_arrhythmias", "compute_hrv"]


def analyze_cardiac(
    processed: ProcessedECG,
    baseline: BreedBaseline,
) -> CardiacAnalysis:
    """Run the full Stage-3 cardiac analysis pipeline.

    1. Classify arrhythmia events from the processed ECG.
    2. Compute HRV metrics from R-R intervals.
    3. Derive heart-rate statistics and assemble the result.

    Parameters
    ----------
    processed:
        Stage-2 output — filtered ECG signal with detected QRS complexes
        and R-R intervals.
    baseline:
        Breed-specific cardiac reference ranges.

    Returns
    -------
    CardiacAnalysis
        Complete analysis results including arrhythmia events, HRV metrics,
        and beat-level statistics.
    """
    # --- Arrhythmia classification ---
    arrhythmia_events = classify_arrhythmias(processed, baseline)

    # --- HRV metrics ---
    rr = processed.rr_intervals_ms
    if len(rr) >= 2:
        hrv_metrics = compute_hrv(rr)
    else:
        # Not enough beats for meaningful HRV — return zeroed metrics.
        from chloe_heart.models import HRVMetrics

        hrv_metrics = HRVMetrics(
            sdnn=0.0,
            rmssd=0.0,
            pnn50=0.0,
            mean_rr=0.0,
            mean_hr=0.0,
        )

    # --- Heart-rate statistics ---
    if len(rr) > 0:
        inst_hr = 60_000.0 / rr  # instantaneous HR from each R-R interval
        mean_hr_bpm = float(np.mean(inst_hr))
        min_hr_bpm = float(np.min(inst_hr))
        max_hr_bpm = float(np.max(inst_hr))
    else:
        mean_hr_bpm = 0.0
        min_hr_bpm = 0.0
        max_hr_bpm = 0.0

    # --- Abnormal beat count ---
    # Count beats flagged as PVC or part of a VT run.
    abnormal_types = {
        ArrhythmiaType.PVC,
        ArrhythmiaType.VT,
        ArrhythmiaType.AFIB,
    }
    abnormal_beat_indices: set[int] = set()
    for event in arrhythmia_events:
        if event.arrhythmia_type in abnormal_types and event.beat_index is not None:
            abnormal_beat_indices.add(event.beat_index)
            # For VT runs, estimate the number of beats involved.
            if (
                event.arrhythmia_type == ArrhythmiaType.VT
                and event.duration_seconds is not None
                and event.duration_seconds > 0
            ):
                # Parse the run length from the details if possible,
                # otherwise estimate from duration and mean RR.
                detail_beats = _parse_run_length(event.details)
                if detail_beats is not None:
                    for offset in range(detail_beats):
                        abnormal_beat_indices.add(event.beat_index + offset)

    abnormal_beat_count = len(abnormal_beat_indices)

    return CardiacAnalysis(
        arrhythmia_events=arrhythmia_events,
        hrv_metrics=hrv_metrics,
        mean_hr_bpm=mean_hr_bpm,
        min_hr_bpm=min_hr_bpm,
        max_hr_bpm=max_hr_bpm,
        total_beats=processed.num_beats,
        duration_seconds=processed.duration_seconds,
        abnormal_beat_count=abnormal_beat_count,
    )


def _parse_run_length(details: str) -> int | None:
    """Extract the number of consecutive PVCs from VT detail text."""
    # Expected format: "5 consecutive PVCs at 200 BPM"
    parts = details.split()
    if len(parts) >= 1:
        try:
            return int(parts[0])
        except ValueError:
            return None
    return None
