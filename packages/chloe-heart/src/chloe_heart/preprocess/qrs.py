"""QRS complex detection using the Pan-Tompkins algorithm.

Reference
---------
Pan, J. & Tompkins, W.J. (1985). "A real-time QRS detection algorithm."
IEEE Transactions on Biomedical Engineering, BME-32(3), 230-236.

The implementation follows the classical five-stage signal-processing chain
described in the paper, with adaptive dual-threshold logic for R-peak
identification.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.signal import butter, filtfilt

from chloe_heart.models import QRSComplex

if TYPE_CHECKING:
    from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _pan_tompkins_bandpass(
    signal: NDArray[np.float64],
    sample_rate: float,
) -> NDArray[np.float64]:
    """Stage 1 — Bandpass filter (5-15 Hz).

    This narrow band retains the steep slopes of the QRS complex while
    suppressing P/T waves, baseline wander, and high-frequency noise.
    """
    nyquist = sample_rate / 2.0
    low = 5.0 / nyquist
    high = min(15.0 / nyquist, 1.0 - 1e-6)
    b, a = butter(2, [low, high], btype="band")
    return filtfilt(b, a, signal).astype(np.float64)


def _differentiate(signal: NDArray[np.float64]) -> NDArray[np.float64]:
    """Stage 2 — Five-point derivative.

    Approximates the slope of the QRS complex.  The five-point formula from
    the original paper is:
        y[n] = (1/8T)(-x[n-2] - 2x[n-1] + 2x[n+1] + x[n+2])
    We use ``np.gradient`` which gives a comparable central-difference result.
    """
    return np.gradient(signal).astype(np.float64)


def _square(signal: NDArray[np.float64]) -> NDArray[np.float64]:
    """Stage 3 — Point-wise squaring.

    Makes all values positive and amplifies the high-slope QRS features
    relative to the lower-amplitude P and T waves.
    """
    return (signal**2).astype(np.float64)


def _moving_window_integration(
    signal: NDArray[np.float64],
    sample_rate: float,
) -> NDArray[np.float64]:
    """Stage 4 — Moving-window integration.

    Window width of ~150 ms matches the typical QRS duration.  The integrator
    smooths the squared-derivative signal into a single pulse per heartbeat.
    """
    window_size = max(1, int(round(0.150 * sample_rate)))
    kernel = np.ones(window_size, dtype=np.float64) / window_size
    integrated = np.convolve(signal, kernel, mode="same")
    return integrated.astype(np.float64)


def _adaptive_threshold(
    integrated: NDArray[np.float64],
    sample_rate: float,
) -> list[int]:
    """Stage 5 — Adaptive dual-threshold peak detection.

    Maintains running estimates of the *signal peak* level (SPKI) and the
    *noise peak* level (NPKI), from which two thresholds are derived:

        THRESHOLD_I  = NPKI + 0.25 * (SPKI - NPKI)
        THRESHOLD_II = 0.5 * THRESHOLD_I

    A refractory period of 200 ms prevents double-detections within the
    same heartbeat.  If no peak is found in a 1.66x expected R-R interval,
    the algorithm performs a *search-back* at the lower threshold.
    """
    # Minimum refractory period between successive R-peaks (200 ms)
    refractory_samples = int(round(0.200 * sample_rate))

    # Initialise peak-level estimates from the first second of data
    init_window = min(len(integrated), int(sample_rate))
    spki = float(np.max(integrated[:init_window]))  # signal peak estimate
    npki = float(0.5 * np.mean(integrated[:init_window]))  # noise peak estimate

    threshold_i = npki + 0.25 * (spki - npki)

    # Find all local maxima in the integrated signal
    peaks: list[int] = []
    for i in range(1, len(integrated) - 1):
        if integrated[i] > integrated[i - 1] and integrated[i] >= integrated[i + 1]:
            peaks.append(i)

    if not peaks:
        return []

    # Walk through candidate peaks and apply adaptive thresholding
    r_peaks: list[int] = []
    rr_recent: list[int] = []  # recent R-R intervals in samples

    for peak_idx in peaks:
        peak_val = float(integrated[peak_idx])

        # Enforce refractory period
        if r_peaks and (peak_idx - r_peaks[-1]) < refractory_samples:
            continue

        if peak_val > threshold_i:
            # --- Accepted as a signal peak ---
            # Search-back check: if the R-R interval is too long (> 1.66x
            # the running average), scan the gap at the lower threshold.
            if r_peaks and rr_recent:
                rr_mean = float(np.mean(rr_recent[-8:]))
                rr_limit = 1.66 * rr_mean
                gap = peak_idx - r_peaks[-1]

                if gap > rr_limit:
                    threshold_ii = 0.5 * threshold_i
                    # Search back through the gap for a missed peak
                    search_start = r_peaks[-1] + refractory_samples
                    search_end = peak_idx
                    best_sb_idx: int | None = None
                    best_sb_val: float = 0.0

                    for sb_idx in range(search_start, search_end):
                        if sb_idx < 1 or sb_idx >= len(integrated) - 1:
                            continue
                        val = float(integrated[sb_idx])
                        if (
                            val > threshold_ii
                            and val > integrated[sb_idx - 1]
                            and val >= integrated[sb_idx + 1]
                            and val > best_sb_val
                        ):
                            best_sb_val = val
                            best_sb_idx = sb_idx

                    if best_sb_idx is not None:
                        r_peaks.append(best_sb_idx)
                        if len(r_peaks) >= 2:
                            rr_recent.append(r_peaks[-1] - r_peaks[-2])
                        # Update signal-peak estimate with the back-filled
                        # peak (weighted toward previous estimate)
                        spki = 0.75 * spki + 0.25 * best_sb_val

            r_peaks.append(peak_idx)
            if len(r_peaks) >= 2:
                rr_recent.append(r_peaks[-1] - r_peaks[-2])
            spki = 0.875 * spki + 0.125 * peak_val
        else:
            # --- Classified as noise ---
            npki = 0.875 * npki + 0.125 * peak_val

        # Recompute threshold from updated estimates
        threshold_i = npki + 0.25 * (spki - npki)

    return sorted(r_peaks)


def _refine_r_peaks(
    r_peak_indices: list[int],
    original_signal: NDArray[np.float64],
    sample_rate: float,
) -> list[int]:
    """Refine R-peak locations on the *original* (unfiltered) signal.

    The Pan-Tompkins chain locates R-peaks in the integrated waveform.
    This step searches a small window around each detected index in the
    original signal to find the true maximum amplitude.
    """
    half_window = max(1, int(round(0.075 * sample_rate)))  # +/- 75 ms
    refined: list[int] = []

    for idx in r_peak_indices:
        lo = max(0, idx - half_window)
        hi = min(len(original_signal), idx + half_window + 1)
        local_max = int(lo + np.argmax(np.abs(original_signal[lo:hi])))
        refined.append(local_max)

    return refined


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def detect_qrs(
    signal: NDArray[np.float64],
    sample_rate: float,
) -> list[QRSComplex]:
    """Detect QRS complexes (R-peaks) in an ECG signal.

    Implements the Pan-Tompkins algorithm:
      1. Bandpass filter  (5-15 Hz)
      2. Differentiate    (five-point derivative)
      3. Square           (amplify QRS slopes)
      4. Moving-window integration (150 ms window)
      5. Adaptive dual-threshold peak detection

    After detection in the processed signal, R-peak locations are refined on
    the original signal to obtain accurate amplitudes and timing.

    Parameters
    ----------
    signal : NDArray
        ECG voltage values (should already be bandpass-filtered for general
        noise removal, but the algorithm applies its own narrow bandpass
        internally).
    sample_rate : float
        Sampling frequency in Hz.

    Returns
    -------
    list[QRSComplex]
        Detected QRS complexes sorted by time.
    """
    if len(signal) < int(sample_rate * 0.5):
        # Signal too short for meaningful detection
        return []

    # --- Pan-Tompkins processing chain ---
    bp = _pan_tompkins_bandpass(signal, sample_rate)
    diff = _differentiate(bp)
    sq = _square(diff)
    mwi = _moving_window_integration(sq, sample_rate)

    # --- Adaptive threshold R-peak detection ---
    raw_peaks = _adaptive_threshold(mwi, sample_rate)

    if not raw_peaks:
        return []

    # --- Refine on original signal ---
    refined_peaks = _refine_r_peaks(raw_peaks, signal, sample_rate)

    # --- Build QRSComplex objects ---
    # Estimate QRS boundaries as +/- ~50 ms from the R-peak
    qrs_half = max(1, int(round(0.050 * sample_rate)))

    complexes: list[QRSComplex] = []
    for idx in refined_peaks:
        onset = max(0, idx - qrs_half)
        offset = min(len(signal) - 1, idx + qrs_half)
        duration_ms = float((offset - onset) / sample_rate * 1000.0)

        complexes.append(
            QRSComplex(
                r_peak_index=idx,
                r_peak_time=float(idx / sample_rate),
                r_peak_amplitude=float(signal[idx]),
                qrs_onset=onset,
                qrs_offset=offset,
                qrs_duration_ms=duration_ms,
            )
        )

    # Deduplicate: if two QRS complexes share the same R-peak index
    # (possible after refinement), keep only the first.
    seen: set[int] = set()
    unique: list[QRSComplex] = []
    for qrs in complexes:
        if qrs.r_peak_index not in seen:
            seen.add(qrs.r_peak_index)
            unique.append(qrs)

    return sorted(unique, key=lambda q: q.r_peak_index)


def compute_rr_intervals(
    qrs_complexes: list[QRSComplex],
    sample_rate: float,
) -> NDArray[np.float64]:
    """Compute R-R intervals in milliseconds from detected QRS complexes.

    Parameters
    ----------
    qrs_complexes : list[QRSComplex]
        Sorted list of detected QRS complexes.
    sample_rate : float
        Sampling frequency in Hz.

    Returns
    -------
    NDArray[np.float64]
        Array of R-R intervals in milliseconds.  Length is
        ``len(qrs_complexes) - 1``.
    """
    if len(qrs_complexes) < 2:
        return np.array([], dtype=np.float64)

    indices = np.array([q.r_peak_index for q in qrs_complexes], dtype=np.float64)
    rr_samples = np.diff(indices)
    rr_ms = rr_samples / sample_rate * 1000.0
    return rr_ms.astype(np.float64)
