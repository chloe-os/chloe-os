"""Rule-based arrhythmia classification for canine ECG.

Dogs have physiological sinus arrhythmia (respiratory-linked R-R variation)
that is perfectly normal.  This module separates benign sinus arrhythmia from
pathological rhythms such as PVCs, AFib, and VT using heuristic rules tuned
for veterinary cardiology.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from chloe_heart.models import (
    ArrhythmiaEvent,
    ArrhythmiaType,
    BreedBaseline,
    ProcessedECG,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def classify_arrhythmias(
    processed: ProcessedECG,
    baseline: BreedBaseline,
) -> list[ArrhythmiaEvent]:
    """Classify arrhythmia events from a processed ECG recording.

    Detection rules applied (in priority order):
    1. Ventricular tachycardia (VT) — 3+ consecutive PVCs at HR > 180 BPM.
    2. Premature ventricular complexes (PVC) — short R-R followed by
       compensatory pause.
    3. Atrial fibrillation (AFib) — irregularly irregular R-R intervals
       over a sliding window.
    4. Sinus arrhythmia — gradual R-R variation correlated with the
       respiratory cycle.  *Normal* in dogs.
    5. Sustained bradycardia / tachycardia — mean HR outside breed range
       for >30 s.

    Parameters
    ----------
    processed:
        Stage-2 output containing the filtered signal, detected QRS
        complexes, and R-R intervals.
    baseline:
        Breed-specific cardiac reference ranges.

    Returns
    -------
    list[ArrhythmiaEvent]
        Chronologically ordered list of detected events.
    """
    rr = processed.rr_intervals_ms
    if len(rr) < 3:
        return []

    events: list[ArrhythmiaEvent] = []

    # Beat times (seconds) for each R-R interval.  The i-th R-R interval
    # starts at the time of the i-th QRS complex.
    beat_times = np.array([q.r_peak_time for q in processed.qrs_complexes], dtype=np.float64)

    # Running average of R-R intervals (centred, window=7).
    running_avg = _running_mean(rr, window=7)

    # --- 1 & 2. PVC detection (also feeds VT) -------------------------
    pvc_indices = _detect_pvcs(rr, running_avg, processed)
    pvc_set = set(pvc_indices)

    for idx in pvc_indices:
        events.append(
            ArrhythmiaEvent(
                arrhythmia_type=ArrhythmiaType.PVC,
                time_seconds=float(beat_times[idx]),
                beat_index=idx,
                confidence=0.8,
                details="Short R-R with compensatory pause",
            )
        )

    # --- 1. VT: 3+ consecutive PVCs at HR > 180 BPM -------------------
    vt_events = _detect_vt(pvc_indices, rr, beat_times)
    events.extend(vt_events)

    # --- 3. AFib: irregularly irregular R-R over sliding window --------
    afib_events = _detect_afib(rr, beat_times)
    events.extend(afib_events)

    # --- 4. Sinus arrhythmia (normal in dogs) --------------------------
    sinus_events = _detect_sinus_arrhythmia(rr, beat_times, pvc_set)
    events.extend(sinus_events)

    # --- 5. Sustained bradycardia / tachycardia ------------------------
    rate_events = _detect_sustained_rate_abnormalities(
        rr, beat_times, baseline, processed.duration_seconds
    )
    events.extend(rate_events)

    # Sort chronologically.
    events.sort(key=lambda e: e.time_seconds)
    return events


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _running_mean(arr: NDArray[np.float64], window: int = 7) -> NDArray[np.float64]:
    """Compute a centred running mean, padding edges with the global mean."""
    n = len(arr)
    if n == 0:
        return arr.copy()
    if n < window:
        return np.full(n, np.mean(arr))
    kernel = np.ones(window) / window
    padded = np.pad(arr, (window // 2, window // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")[:n]


def _detect_pvcs(
    rr: NDArray[np.float64],
    running_avg: NDArray[np.float64],
    processed: ProcessedECG,
) -> list[int]:
    """Detect PVC beats by short R-R + compensatory pause pattern.

    A PVC is flagged when:
    - The R-R interval is < 80 % of the running average, **and**
    - The *next* R-R interval is > 120 % of the running average
      (compensatory pause).
    - Optionally, QRS duration is wider than normal (> 80 ms).
    """
    pvc_indices: list[int] = []
    n = len(rr)
    for i in range(n - 1):
        short_rr = rr[i] < 0.80 * running_avg[i]
        long_next = rr[i + 1] > 1.20 * running_avg[i + 1]

        if short_rr and long_next:
            # Optional QRS-width check.
            qrs = processed.qrs_complexes[i]
            if qrs.qrs_duration_ms is not None and qrs.qrs_duration_ms <= 60:
                # Narrow QRS — less likely a true PVC, skip.
                continue
            pvc_indices.append(i)

    return pvc_indices


def _detect_vt(
    pvc_indices: list[int],
    rr: NDArray[np.float64],
    beat_times: NDArray[np.float64],
) -> list[ArrhythmiaEvent]:
    """Detect ventricular tachycardia: 3+ consecutive PVCs at HR > 180 BPM."""
    events: list[ArrhythmiaEvent] = []
    if len(pvc_indices) < 3:
        return events

    # Find runs of consecutive PVC indices.
    run_start = 0
    for i in range(1, len(pvc_indices) + 1):
        # End of a consecutive run?
        if i == len(pvc_indices) or pvc_indices[i] != pvc_indices[i - 1] + 1:
            run_len = i - run_start
            if run_len >= 3:
                run_idx = pvc_indices[run_start:i]
                # Compute instantaneous HR over the run.
                run_rr = rr[run_idx]
                mean_rr_s = float(np.mean(run_rr)) / 1000.0
                inst_hr = 60.0 / mean_rr_s if mean_rr_s > 0 else 0.0

                if inst_hr > 180.0:
                    start_time = float(beat_times[run_idx[0]])
                    end_time = float(beat_times[run_idx[-1]])
                    events.append(
                        ArrhythmiaEvent(
                            arrhythmia_type=ArrhythmiaType.VT,
                            time_seconds=start_time,
                            duration_seconds=end_time - start_time,
                            beat_index=run_idx[0],
                            confidence=0.85,
                            details=(f"{run_len} consecutive PVCs at {inst_hr:.0f} BPM"),
                        )
                    )
            run_start = i

    return events


def _detect_afib(
    rr: NDArray[np.float64],
    beat_times: NDArray[np.float64],
    window: int = 30,
    cv_threshold: float = 0.15,
) -> list[ArrhythmiaEvent]:
    """Detect atrial fibrillation via coefficient of variation of R-R.

    AFib is flagged when:
    - CV of R-R in a sliding window exceeds *cv_threshold* (0.15), **and**
    - Variation is *irregularly* irregular (no clear periodic pattern).

    The second criterion distinguishes AFib from sinus arrhythmia, which
    has a gradual, quasi-periodic modulation.
    """
    events: list[ArrhythmiaEvent] = []
    n = len(rr)
    if n < window:
        return events

    i = 0
    while i <= n - window:
        segment = rr[i : i + window]
        mean_rr = float(np.mean(segment))
        if mean_rr == 0:
            i += 1
            continue
        cv = float(np.std(segment)) / mean_rr

        if cv > cv_threshold and _is_irregularly_irregular(segment):
            start_time = float(beat_times[i])
            end_time = float(beat_times[min(i + window, len(beat_times) - 1)])
            events.append(
                ArrhythmiaEvent(
                    arrhythmia_type=ArrhythmiaType.AFIB,
                    time_seconds=start_time,
                    duration_seconds=end_time - start_time,
                    beat_index=i,
                    confidence=0.7,
                    details=f"R-R CV={cv:.3f} over {window}-beat window",
                )
            )
            # Advance past this window to avoid duplicate detections.
            i += window
        else:
            i += 1

    return events


def _is_irregularly_irregular(rr_segment: NDArray[np.float64]) -> bool:
    """Return True if R-R variation lacks a periodic (respiratory) pattern.

    Sinus arrhythmia produces a smooth, quasi-sinusoidal modulation.
    AFib produces chaotic variation.  We distinguish them by checking
    whether the autocorrelation of successive differences has a clear
    peak — if it does, the variation is periodic (sinus arrhythmia).
    """
    diffs = np.diff(rr_segment)
    if len(diffs) < 6:
        return True  # Too short to judge; default to irregular.

    # Normalise.
    diffs = diffs - np.mean(diffs)
    std = float(np.std(diffs))
    if std < 1e-9:
        return False  # No variation at all.
    diffs = diffs / std

    # Autocorrelation (unbiased).
    n = len(diffs)
    acf = np.correlate(diffs, diffs, mode="full")[n - 1 :]
    acf = acf / acf[0]  # Normalise so lag-0 = 1.

    # A respiratory pattern would show a secondary peak at lag ~3-8 beats.
    secondary = acf[3:9] if len(acf) > 8 else acf[1:]
    if len(secondary) == 0:
        return True
    max_secondary = float(np.max(secondary))

    # If the best secondary peak is weak, the rhythm is irregularly irregular.
    return max_secondary < 0.3


def _detect_sinus_arrhythmia(
    rr: NDArray[np.float64],
    beat_times: NDArray[np.float64],
    pvc_set: set[int],
    window: int = 20,
    variation_threshold: float = 0.05,
) -> list[ArrhythmiaEvent]:
    """Detect sinus arrhythmia (normal respiratory-linked R-R variation).

    Sinus arrhythmia is characterised by *gradual*, periodic R-R changes
    — the heart speeds up during inspiration and slows during expiration.

    We flag it when:
    - R-R coefficient of variation exceeds a minimum threshold (>5 %),
    - Successive R-R differences are *small* (gradual), and
    - The segment is not already classified as PVC or AFib.
    """
    events: list[ArrhythmiaEvent] = []
    n = len(rr)
    if n < window:
        return events

    i = 0
    while i <= n - window:
        segment = rr[i : i + window]
        mean_rr = float(np.mean(segment))
        if mean_rr == 0:
            i += 1
            continue

        cv = float(np.std(segment)) / mean_rr

        # Must have meaningful variation.
        if cv < variation_threshold:
            i += window
            continue

        # Check that no beats in this window are PVCs.
        window_indices = set(range(i, i + window))
        if window_indices & pvc_set:
            i += 1
            continue

        # Gradual change: median absolute successive difference should be
        # small relative to mean R-R.
        abs_diffs = np.abs(np.diff(segment))
        median_diff = float(np.median(abs_diffs))
        gradual = median_diff < 0.10 * mean_rr

        if gradual and not _is_irregularly_irregular(segment):
            start_time = float(beat_times[i])
            end_time = float(beat_times[min(i + window - 1, len(beat_times) - 1)])
            events.append(
                ArrhythmiaEvent(
                    arrhythmia_type=ArrhythmiaType.SINUS_ARRHYTHMIA,
                    time_seconds=start_time,
                    duration_seconds=end_time - start_time,
                    beat_index=i,
                    confidence=0.9,
                    details="Gradual R-R variation — normal in dogs",
                )
            )
            i += window
        else:
            i += 1

    return events


def _detect_sustained_rate_abnormalities(
    rr: NDArray[np.float64],
    beat_times: NDArray[np.float64],
    baseline: BreedBaseline,
    total_duration: float,
    sustain_seconds: float = 30.0,
) -> list[ArrhythmiaEvent]:
    """Detect sustained bradycardia or tachycardia.

    The mean HR must remain outside the breed-specific resting range for
    at least *sustain_seconds* continuously.
    """
    events: list[ArrhythmiaEvent] = []
    n = len(rr)
    if n < 2:
        return events

    # Convert R-R intervals to instantaneous HR (BPM).
    inst_hr = 60_000.0 / rr  # rr is in ms

    # Walk through beats, tracking sustained periods.
    segment_start: int | None = None
    current_type: ArrhythmiaType | None = None

    for i in range(n):
        hr_val = float(inst_hr[i])

        if hr_val < baseline.hr_min:
            beat_type = ArrhythmiaType.BRADYCARDIA
        elif hr_val > baseline.hr_max:
            beat_type = ArrhythmiaType.TACHYCARDIA
        else:
            beat_type = None

        if beat_type is not None and beat_type == current_type:
            # Continue the current segment — check duration.
            elapsed = float(beat_times[i] - beat_times[segment_start])  # type: ignore[index]
            if elapsed >= sustain_seconds:
                # Emit the event once per sustained segment.
                events.append(
                    ArrhythmiaEvent(
                        arrhythmia_type=beat_type,
                        time_seconds=float(beat_times[segment_start]),  # type: ignore[index]
                        duration_seconds=elapsed,
                        beat_index=segment_start,
                        confidence=0.9,
                        details=(
                            f"Sustained {beat_type.value} for "
                            f"{elapsed:.1f}s (breed range "
                            f"{baseline.hr_min}–{baseline.hr_max} BPM)"
                        ),
                    )
                )
                # Reset so we don't keep emitting for the same segment.
                segment_start = i
        elif beat_type is not None:
            # New segment.
            segment_start = i
            current_type = beat_type
        else:
            # Normal HR — reset.
            segment_start = None
            current_type = None

    return events
