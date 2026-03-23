# Sample Canine ECG Data

This directory contains **synthetic example data** for testing Chloe Heart.

**This is NOT a real ECG recording from a real dog.** The signal is computationally generated to demonstrate the pipeline's functionality.

## Files

- `ecg.csv` — Synthetic ECG signal (30 seconds, 500 Hz, ~45 beats at 90 BPM)

## Signal Characteristics

- Normal sinus rhythm with respiratory sinus arrhythmia (normal in dogs)
- One premature ventricular complex (PVC) at approximately beat 15
- Heart rate ~90 BPM (typical for a medium-sized dog)
- P waves, QRS complexes, and T waves present
- Low-level background noise added for realism

## Usage

```bash
chloe heart analyze --ecg examples/sample_canine_ecg/ecg.csv --breed staffy --output cardiac_report.html
```
