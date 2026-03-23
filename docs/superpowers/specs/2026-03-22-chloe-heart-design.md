# chloe-heart — Canine Cardiac Monitoring Analysis

## Context

Cardiac disease is the second leading cause of death in dogs after cancer. Many conditions (DCM, arrhythmias, early CHF) go undetected until a crisis because continuous cardiac monitoring and expert interpretation are inaccessible to most pet owners. Meanwhile, pet wearables and portable ECG devices are becoming more affordable — the gap is in making the data actionable.

chloe-heart processes ECG recordings and heart rate data from monitors and wearables, detects arrhythmias, computes cardiac health metrics, and generates a plain-language report a pet owner can share with their veterinarian.

## Architecture

New package `packages/chloe-heart/` in the existing chloe-os monorepo, following the same patterns as chloe-core.

### Module Structure

```
packages/chloe-heart/
├── src/chloe_heart/
│   ├── __init__.py
│   ├── models.py              # Shared dataclasses
│   ├── ingest/                # Stage 1: signal ingestion
│   │   ├── __init__.py
│   │   ├── edf.py             # EDF format parser
│   │   ├── csv_reader.py      # CSV (time, voltage) parser
│   │   ├── wfdb_reader.py     # PhysioNet WFDB format
│   │   └── wearable.py        # HR/HRV from consumer devices
│   ├── preprocess/            # Stage 2: signal preprocessing
│   │   ├── __init__.py
│   │   ├── filters.py         # Bandpass, notch filtering
│   │   └── qrs.py             # QRS detection, R-R extraction
│   ├── analysis/              # Stage 3: cardiac analysis
│   │   ├── __init__.py
│   │   ├── arrhythmia.py      # Arrhythmia classification
│   │   ├── hrv.py             # HRV metric computation
│   │   └── trends.py          # Long-term trend analysis
│   ├── scoring/               # Stage 4: health scoring
│   │   ├── __init__.py
│   │   └── scorer.py          # Composite cardiac health score
│   ├── interpretation/        # Stage 5: AI interpretation
│   │   ├── __init__.py
│   │   └── interpreter.py     # LLM plain-language explanation
│   ├── report/                # Stage 6: report generation
│   │   ├── __init__.py
│   │   ├── generator.py       # HTML report builder
│   │   └── template.html      # Jinja2 template with ECG plots
│   └── community/             # Anonymized data sharing
│       ├── __init__.py
│       └── sharing.py
├── tests/
└── pyproject.toml
```

## Pipeline — 6 Stages

### Stage 1: Signal Ingestion (`chloe_heart.ingest`)

**Input:** ECG file or HR/HRV export
**Output:** `ECGSignal` or `HRTimeSeries` dataclass

Supported formats:
- **EDF** (European Data Format) — standard for clinical Holter monitors. Parsed via `pyedflib`.
- **CSV** — time + voltage columns. Common export from AliveCor Vet and research devices.
- **WFDB** (PhysioNet) — research standard. Parsed via `wfdb` library.
- **Wearable JSON/CSV** — processed HR/HRV data from PetPace, Fi, Whistle. Schema detection with fallback to column mapping.

Auto-detects format from file extension and header inspection.

### Stage 2: Preprocessing (`chloe_heart.preprocess`)

**Input:** `ECGSignal`
**Output:** `ProcessedECG` with detected beats and R-R intervals

- **Bandpass filter:** 0.5–40 Hz (canine ECG frequency range) via scipy Butterworth filter
- **Notch filter:** 50/60 Hz powerline noise removal
- **Baseline wander correction:** High-pass filter or wavelet detrending
- **QRS detection:** Pan-Tompkins algorithm adapted for canine ECG (different QRS morphology)
- **R-R interval extraction:** Time between consecutive R-peaks
- **Beat segmentation:** Extract individual heartbeat waveforms for classification

For HR/HRV data (no raw ECG), this stage is bypassed — data goes directly to analysis.

### Stage 3: Analysis (`chloe_heart.analysis`)

**Input:** `ProcessedECG` or `HRTimeSeries`
**Output:** `CardiacAnalysis` with arrhythmia events, HRV metrics, trends

#### Arrhythmia Detection (rule-based + ML)

Rule-based classification (MVP):
- **Premature Ventricular Complexes (PVCs):** R-R interval < 80% of running average, wide QRS
- **Atrial Fibrillation (AFib):** Irregularly irregular R-R intervals, absent P-waves
- **Ventricular Tachycardia (VT):** ≥3 consecutive PVCs, HR > 180 BPM
- **Bradycardia:** HR < breed-specific lower bound
- **Tachycardia:** HR > breed-specific upper bound
- **Sinus Arrhythmia:** NORMAL in dogs — R-R varies with respiration. Must NOT flag this.

ML model (v0.2):
- Fine-tune a pre-trained 1D-CNN ECG model on canine data
- Start with transfer learning from human ECG models (MIT-BIH trained)
- Community-contributed labeled data improves the model over time

#### HRV Metrics

Time-domain:
- **SDNN** — Standard deviation of R-R intervals (overall HRV)
- **RMSSD** — Root mean square of successive R-R differences (parasympathetic)
- **pNN50** — % of R-R intervals differing by > 50ms

Frequency-domain:
- **LF/HF ratio** — Sympathetic/parasympathetic balance
- **Total power** — Overall autonomic activity

#### Trend Analysis

For multi-session data:
- Resting heart rate trends over days/weeks
- HRV trends (declining SDNN may indicate worsening cardiac function)
- Exercise recovery rate
- Arrhythmia burden over time (% of abnormal beats per session)

### Stage 4: Health Scoring (`chloe_heart.scoring`)

**Input:** `CardiacAnalysis`
**Output:** `CardiacHealthScore`

Composite score (0-100) based on:
- **Rhythm regularity** (weight 0.30): Arrhythmia burden, type severity
- **HRV health** (weight 0.25): SDNN and RMSSD relative to breed norms
- **Heart rate** (weight 0.20): Resting HR relative to breed-specific range
- **Trend trajectory** (weight 0.15): Improving, stable, or declining
- **Exercise tolerance** (weight 0.10): Recovery rate after activity

Risk flags:
- RED: VT detected, high arrhythmia burden (>10%), rapid HRV decline
- YELLOW: Frequent PVCs, HRV below breed norm, elevated resting HR
- GREEN: Normal rhythm, healthy HRV, stable trends

### Stage 5: AI Interpretation (`chloe_heart.interpretation`)

**Input:** `CardiacHealthScore` + `CardiacAnalysis`
**Output:** `CardiacInterpretation`

Same pattern as chloe-core:
- Serialize findings into structured prompt
- LLM produces plain-language explanation of cardiac status
- Per-finding explanations ("We detected 47 premature beats in 24 hours — here's what that means...")
- Trend interpretation ("Your dog's resting heart rate has been gradually increasing over the past 2 weeks")
- Urgency assessment and suggested next steps
- Clearly labeled as AI-generated, not medical advice
- Skippable — works without API key

### Stage 6: Report Generation (`chloe_heart.report`)

**Input:** All pipeline outputs
**Output:** Self-contained HTML report

Sections:
- **Executive summary** — cardiac health score, risk level, key findings
- **ECG Waveform Viewer** — inline SVG plots of ECG with annotated beats (normal in green, PVCs in red, etc.)
- **Arrhythmia Event Log** — timeline of detected events with timestamps
- **HRV Dashboard** — SDNN, RMSSD, pNN50 visualized against breed norms
- **Heart Rate Trends** — resting HR over time with trend line
- **Health Score Breakdown** — spider/radar chart of scoring components
- **Vet-Ready Summary** — formatted for sharing with a cardiologist
- **Disclaimers** — not a substitute for veterinary cardiology

Uses inline SVG for all visualizations (no external JS libraries needed).

## Breed-Specific Baselines

| Size | Example Breeds | Normal HR (BPM) | Normal SDNN (ms) | DCM Risk |
|------|---------------|-----------------|-------------------|----------|
| Giant | Great Dane, Irish Wolfhound, Newfoundland | 60–100 | 80–200 | High |
| Large | Labrador, Golden Retriever, Boxer, Doberman | 60–120 | 60–150 | Moderate–High |
| Medium | Staffy, Beagle, Cocker Spaniel, Bulldog | 70–130 | 50–130 | Low–Moderate |
| Small | Cavalier KCS, Miniature Poodle, Dachshund | 80–140 | 40–100 | Low (MMVD risk) |
| Toy | Chihuahua, Yorkshire Terrier, Pomeranian | 90–160 | 30–80 | Low |

Key breed-specific conditions:
- **Doberman:** High DCM prevalence, VPCs are early indicator
- **Boxer:** Arrhythmogenic right ventricular cardiomyopathy (ARVC)
- **Cavalier King Charles:** Nearly universal MMVD by age 10
- **Great Dane:** DCM with atrial fibrillation

## CLI Integration

Adds `heart` subcommand group to existing `chloe` CLI:

```bash
# Guided mode
chloe heart guided

# ECG analysis
chloe heart analyze --ecg recording.edf --breed labrador --output report.html

# HR/HRV wearable data
chloe heart analyze --hr-data wearable_export.csv --breed cavalier --output report.html

# Multi-session trend analysis
chloe heart trends --data-dir ./recordings/ --breed boxer --output trends_report.html
```

## Community Data Sharing

Same opt-in model as chloe-core:

**Shared (anonymized):**
- Breed and size category
- Arrhythmia event counts and types
- HRV metrics summary
- Health score
- Device type used

**Never shared:**
- Raw ECG signals
- Owner/pet identifying info
- Location data

Community data enables:
- Breed-specific cardiac baseline refinement
- Arrhythmia prevalence studies
- Training data for ML models (with separate consent)

## Dependencies

```toml
[project]
dependencies = [
    "numpy>=1.24.0",
    "scipy>=1.11.0",
    "jinja2>=3.1.0",
    "pydantic>=2.0.0",
]

[project.optional-dependencies]
edf = ["pyedflib>=0.1.35"]
wfdb = ["wfdb>=4.1.0"]
ml = ["torch>=2.0.0", "transformers>=4.40.0"]
interpretation = ["anthropic>=0.40.0"]
```

## Ethical Considerations

- Same framework as chloe-core: this is a screening tool, not a diagnostic
- ECG interpretation requires veterinary cardiology expertise for definitive diagnosis
- False negatives are dangerous — disclaimers must emphasize that a "normal" result does not guarantee absence of disease
- False positives cause anxiety — explain that flagged events need professional confirmation
- Sinus arrhythmia handling is critical — flagging normal canine respiratory sinus arrhythmia as abnormal would undermine trust

## MVP Scope

**v0.1:**
- Stages 1-2 (ECG ingestion from CSV + preprocessing)
- Stage 3 (rule-based arrhythmia detection + HRV metrics)
- Stage 4 (health scoring)
- Stage 6 (HTML report with ECG plots)
- Breed-specific baselines
- CLI integration (`chloe heart analyze`, `chloe heart guided`)
- Disclaimers and ethics section

**v0.2:**
- EDF and WFDB format support
- Wearable HR/HRV data ingestion
- ML arrhythmia classification model
- AI interpretation (Stage 5)
- Multi-session trend analysis
- Community data sharing

## Verification

1. **Unit tests:** `uv run pytest packages/chloe-heart/tests/`
2. **Integration test with sample ECG:** `chloe heart analyze --ecg examples/sample_canine_ecg/ecg.csv --breed labrador --output test_cardiac_report.html`
3. **Report inspection:** Open HTML report, verify ECG plots render and arrhythmia events are annotated
4. **Breed baseline validation:** Test that breed-specific thresholds produce correct risk flags
