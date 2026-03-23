"""Heart rate variability (HRV) metric computation.

Computes both time-domain and frequency-domain HRV metrics from a series
of R-R intervals.  Frequency-domain analysis uses Welch's method and is
attempted when SciPy is available and the recording is long enough.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from chloe_heart.models import HRVMetrics

if TYPE_CHECKING:
    from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_hrv(rr_intervals_ms: NDArray[np.float64]) -> HRVMetrics:
    """Compute HRV metrics from an array of R-R intervals in milliseconds.

    Parameters
    ----------
    rr_intervals_ms:
        1-D array of successive R-R intervals (ms).  Must contain at least
        two values.

    Returns
    -------
    HRVMetrics
        Time-domain metrics are always populated.  Frequency-domain fields
        are filled when the recording is long enough (>= ~60 s of data) and
        SciPy is installed; otherwise they remain ``None``.

    Raises
    ------
    ValueError
        If fewer than two R-R intervals are provided.
    """
    rr = np.asarray(rr_intervals_ms, dtype=np.float64)
    if len(rr) < 2:
        raise ValueError(f"At least 2 R-R intervals are required; got {len(rr)}")

    # ---- Time-domain metrics ----
    sdnn = float(np.std(rr, ddof=0))
    successive_diffs = np.diff(rr)
    rmssd = float(np.sqrt(np.mean(successive_diffs**2)))
    pnn50 = float(np.sum(np.abs(successive_diffs) > 50.0) / len(successive_diffs) * 100.0)
    mean_rr = float(np.mean(rr))
    mean_hr = 60_000.0 / mean_rr if mean_rr > 0 else 0.0

    # ---- Frequency-domain metrics (optional) ----
    lf_power, hf_power, lf_hf_ratio, total_power = _frequency_domain(rr)

    return HRVMetrics(
        sdnn=sdnn,
        rmssd=rmssd,
        pnn50=pnn50,
        mean_rr=mean_rr,
        mean_hr=mean_hr,
        lf_power=lf_power,
        hf_power=hf_power,
        lf_hf_ratio=lf_hf_ratio,
        total_power=total_power,
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# Frequency band limits (Hz).
_LF_LOW = 0.04
_LF_HIGH = 0.15
_HF_LOW = 0.15
_HF_HIGH = 0.40


def _frequency_domain(
    rr_ms: NDArray[np.float64],
) -> tuple[float | None, float | None, float | None, float | None]:
    """Compute LF and HF power via Welch's periodogram.

    R-R intervals are first interpolated onto a uniform time grid (4 Hz)
    so that standard spectral estimation can be applied.

    Returns ``(None, None, None, None)`` when:
    - SciPy is not installed.
    - The recording is too short (< 60 s).
    """
    try:
        from scipy import interpolate as scipy_interp  # noqa: WPS433
        from scipy import signal as scipy_signal  # noqa: WPS433
    except ImportError:
        return None, None, None, None

    # Cumulative time axis in seconds.
    rr_s = rr_ms / 1000.0
    t_rr = np.cumsum(rr_s)
    t_rr = t_rr - t_rr[0]  # start at 0

    total_duration = float(t_rr[-1])
    if total_duration < 60.0:
        return None, None, None, None

    # Interpolate onto a uniform 4 Hz grid.
    fs = 4.0  # Hz
    t_uniform = np.arange(0.0, total_duration, 1.0 / fs)
    interp_fn = scipy_interp.interp1d(t_rr, rr_ms, kind="cubic", fill_value="extrapolate")
    rr_uniform = interp_fn(t_uniform)

    # Remove mean (detrend).
    rr_uniform = rr_uniform - np.mean(rr_uniform)

    # Welch PSD.
    nperseg = min(256, len(rr_uniform))
    freqs, psd = scipy_signal.welch(
        rr_uniform,
        fs=fs,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        scaling="density",
    )

    # Band powers (trapezoidal integration).
    lf_mask = (freqs >= _LF_LOW) & (freqs < _LF_HIGH)
    hf_mask = (freqs >= _HF_LOW) & (freqs < _HF_HIGH)

    lf_power = float(np.trapz(psd[lf_mask], freqs[lf_mask])) if np.any(lf_mask) else 0.0
    hf_power = float(np.trapz(psd[hf_mask], freqs[hf_mask])) if np.any(hf_mask) else 0.0
    total_power = float(np.trapz(psd, freqs))

    lf_hf_ratio = lf_power / hf_power if hf_power > 0 else None

    return lf_power, hf_power, lf_hf_ratio, total_power
