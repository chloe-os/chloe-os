"""CSV ECG file reader for chloe-heart Stage 1 ingestion.

Parses simple two-column CSV files where the first column is time (seconds)
and the second column is voltage (mV). Supports files with or without a
header row.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from chloe_heart.models import ECGSignal

if TYPE_CHECKING:
    from numpy.typing import NDArray


def _is_header_row(row: list[str]) -> bool:
    """Detect whether a CSV row is a header (non-numeric first cell)."""
    if not row:
        return False
    try:
        float(row[0])
        return False
    except ValueError:
        return True


def _infer_sample_rate(time_values: NDArray[np.float64]) -> float:
    """Infer sample rate from the time column by computing the median time step.

    Uses the median rather than the mean to be robust against occasional
    timing jitter or missing samples.

    Raises
    ------
    ValueError
        If fewer than two time points are provided or the computed sample rate
        is non-positive.
    """
    if len(time_values) < 2:
        raise ValueError("Cannot infer sample rate: need at least two time points.")

    dt = np.median(np.diff(time_values))

    if dt <= 0:
        raise ValueError(
            f"Cannot infer sample rate: median time step is {dt:.6e} s "
            "(must be positive). Check that the time column is monotonically "
            "increasing."
        )

    return float(1.0 / dt)


def load_csv_ecg(
    file_path: str,
    sample_rate: float | None = None,
) -> ECGSignal:
    """Load an ECG recording from a two-column CSV file.

    Parameters
    ----------
    file_path : str
        Path to the CSV file. Expected columns: time (seconds), voltage (mV).
    sample_rate : float or None
        Sampling frequency in Hz. If ``None``, the rate is inferred from the
        time column using the median inter-sample interval.

    Returns
    -------
    ECGSignal
        Parsed signal dataclass ready for downstream processing.

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    ValueError
        If the file is empty, has fewer than two columns, or the sample rate
        cannot be determined.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"ECG file not found: {file_path}")

    time_values: list[float] = []
    voltage_values: list[float] = []

    with open(path, newline="") as fh:
        reader = csv.reader(fh)

        for row in reader:
            # Skip blank lines
            if not row or all(cell.strip() == "" for cell in row):
                continue

            # Auto-detect and skip a header row
            if _is_header_row(row):
                continue

            if len(row) < 2:
                raise ValueError(
                    f"Expected at least 2 columns (time, voltage), got {len(row)} in row: {row}"
                )

            time_values.append(float(row[0].strip()))
            voltage_values.append(float(row[1].strip()))

    if not voltage_values:
        raise ValueError(f"No data rows found in {file_path}")

    time_arr = np.asarray(time_values, dtype=np.float64)
    voltage_arr = np.asarray(voltage_values, dtype=np.float64)

    # Determine sample rate
    if sample_rate is None:
        sample_rate = _infer_sample_rate(time_arr)

    duration = float(time_arr[-1] - time_arr[0])
    # Ensure duration accounts for the last sample interval
    if len(time_arr) > 1:
        duration += float(np.median(np.diff(time_arr)))

    return ECGSignal(
        signal=voltage_arr,
        sample_rate=sample_rate,
        duration_seconds=duration,
        source_file=os.path.basename(file_path),
        source_format="csv",
        channel_name="ECG",
        units="mV",
    )
