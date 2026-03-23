"""Stage 1 — ECG signal ingestion.

Provides format-aware loaders that return a unified ``ECGSignal`` dataclass
regardless of the source file format.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from chloe_heart.ingest.csv_reader import load_csv_ecg

if TYPE_CHECKING:
    from chloe_heart.models import ECGSignal

__all__ = ["load_csv_ecg", "load_ecg"]

# Registry mapping file extensions to their loader functions.
_LOADERS: dict[str, object] = {
    ".csv": load_csv_ecg,
}


def load_ecg(
    file_path: str,
    sample_rate: float | None = None,
) -> ECGSignal:
    """Load an ECG recording, auto-detecting the format from the file extension.

    Parameters
    ----------
    file_path : str
        Path to the ECG file.  Supported extensions: ``.csv``.
    sample_rate : float or None
        Sampling frequency in Hz.  Passed through to the format-specific
        loader.  If ``None``, the loader will try to infer it from the data.

    Returns
    -------
    ECGSignal
        Parsed signal ready for preprocessing.

    Raises
    ------
    ValueError
        If the file extension is not recognised.
    """
    import os

    ext = os.path.splitext(file_path)[1].lower()

    if ext not in _LOADERS:
        supported = ", ".join(sorted(_LOADERS.keys()))
        raise ValueError(f"Unsupported ECG file format '{ext}'. Supported formats: {supported}")

    loader = _LOADERS[ext]
    return loader(file_path, sample_rate=sample_rate)  # type: ignore[operator]
