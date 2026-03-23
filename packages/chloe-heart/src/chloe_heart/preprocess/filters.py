"""Signal filtering utilities for ECG preprocessing.

Provides Butterworth bandpass and IIR notch filters implemented with
scipy.signal for cleaning raw canine ECG recordings.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.signal import butter, filtfilt, iirnotch

if TYPE_CHECKING:
    from numpy.typing import NDArray


def bandpass_filter(
    signal: NDArray[np.float64],
    sample_rate: float,
    low: float = 0.5,
    high: float = 40.0,
    order: int = 4,
) -> NDArray[np.float64]:
    """Apply a zero-phase Butterworth bandpass filter.

    The default passband of 0.5 -- 40 Hz is appropriate for canine ECG:
    it removes baseline wander (< 0.5 Hz) and high-frequency noise / EMG
    artefacts (> 40 Hz) while preserving the clinically relevant QRS, P, and
    T wave morphology.

    Parameters
    ----------
    signal : NDArray
        Raw ECG voltage samples.
    sample_rate : float
        Sampling frequency in Hz.
    low : float
        Lower cutoff frequency in Hz.
    high : float
        Upper cutoff frequency in Hz.
    order : int
        Filter order.  Higher orders give a steeper roll-off but may
        introduce ringing on short signals.

    Returns
    -------
    NDArray
        Filtered signal (same length as input).
    """
    nyquist = sample_rate / 2.0
    low_norm = low / nyquist
    high_norm = high / nyquist

    # Clamp to valid Nyquist range (0, 1) exclusive
    low_norm = max(low_norm, 1e-6)
    high_norm = min(high_norm, 1.0 - 1e-6)

    b, a = butter(order, [low_norm, high_norm], btype="band")
    return filtfilt(b, a, signal).astype(np.float64)


def notch_filter(
    signal: NDArray[np.float64],
    sample_rate: float,
    freq: float = 50.0,
    quality: float = 30.0,
) -> NDArray[np.float64]:
    """Apply a zero-phase IIR notch filter to remove powerline interference.

    Parameters
    ----------
    signal : NDArray
        ECG voltage samples (may already be bandpass-filtered).
    sample_rate : float
        Sampling frequency in Hz.
    freq : float
        Centre frequency of the notch in Hz (50 Hz for most of the world,
        60 Hz in North America / parts of Asia).
    quality : float
        Quality factor *Q* — controls notch width.  Higher values give a
        narrower notch.  30 is a reasonable default for powerline rejection.

    Returns
    -------
    NDArray
        Filtered signal (same length as input).
    """
    nyquist = sample_rate / 2.0

    # If the notch frequency is at or above Nyquist, the filter is
    # inapplicable — return the signal unchanged rather than erroring.
    if freq >= nyquist:
        return signal.copy()

    b, a = iirnotch(freq, quality, fs=sample_rate)
    return filtfilt(b, a, signal).astype(np.float64)
