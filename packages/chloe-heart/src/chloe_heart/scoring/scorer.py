"""Stage 4 — Cardiac health scoring.

Computes a composite 0-100 cardiac health score from the Stage-3 analysis
results and breed-specific baselines, then assigns a risk level
(GREEN / YELLOW / RED).
"""

from __future__ import annotations

from chloe_heart.models import (
    ArrhythmiaType,
    BreedBaseline,
    CardiacAnalysis,
    CardiacHealthScore,
    RiskLevel,
)

# Component weights (must sum to 1.0).
_WEIGHT_RHYTHM = 0.30
_WEIGHT_HRV = 0.25
_WEIGHT_HR = 0.20
_WEIGHT_TREND = 0.15
_WEIGHT_EXERCISE = 0.10

_WEIGHTS = {
    "rhythm": _WEIGHT_RHYTHM,
    "hrv": _WEIGHT_HRV,
    "heart_rate": _WEIGHT_HR,
    "trend": _WEIGHT_TREND,
    "exercise": _WEIGHT_EXERCISE,
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def score_cardiac_health(
    analysis: CardiacAnalysis,
    baseline: BreedBaseline,
    *,
    trend_score_override: float | None = None,
    exercise_score_override: float | None = None,
) -> CardiacHealthScore:
    """Compute a composite cardiac health score.

    Each component yields a sub-score in the range ``[0, 1]``.  The weighted
    combination is scaled to ``[0, 100]``.

    Parameters
    ----------
    analysis:
        Stage-3 cardiac analysis output.
    baseline:
        Breed-specific cardiac reference ranges.
    trend_score_override:
        Optional override for the trend sub-score (0-1).  Defaults to 0.7
        for single-session data (no trend history).
    exercise_score_override:
        Optional override for the exercise sub-score (0-1).  Defaults to
        0.7 for MVP (no exercise data).

    Returns
    -------
    CardiacHealthScore
        Composite score, component sub-scores, risk level, and risk flags.
    """
    rhythm = _rhythm_score(analysis)
    hrv = _hrv_score(analysis, baseline)
    hr = _heart_rate_score(analysis, baseline)
    trend = trend_score_override if trend_score_override is not None else 0.7
    exercise = exercise_score_override if exercise_score_override is not None else 0.7

    # Clamp all sub-scores to [0, 1].
    rhythm = _clamp(rhythm)
    hrv = _clamp(hrv)
    hr = _clamp(hr)
    trend = _clamp(trend)
    exercise = _clamp(exercise)

    overall = (
        rhythm * _WEIGHT_RHYTHM
        + hrv * _WEIGHT_HRV
        + hr * _WEIGHT_HR
        + trend * _WEIGHT_TREND
        + exercise * _WEIGHT_EXERCISE
    ) * 100.0

    risk_flags = _compute_risk_flags(analysis, baseline, overall)
    risk_level = _determine_risk_level(analysis, risk_flags, overall)

    return CardiacHealthScore(
        overall_score=round(overall, 1),
        risk_level=risk_level,
        rhythm_score=round(rhythm, 3),
        hrv_score=round(hrv, 3),
        heart_rate_score=round(hr, 3),
        trend_score=round(trend, 3),
        exercise_score=round(exercise, 3),
        risk_flags=risk_flags,
        weights=_WEIGHTS,
    )


# ---------------------------------------------------------------------------
# Component scoring
# ---------------------------------------------------------------------------


def _rhythm_score(analysis: CardiacAnalysis) -> float:
    """Score rhythm regularity based on arrhythmia burden.

    - 0 % arrhythmia burden -> 1.0
    - >= 10 % arrhythmia burden -> 0.0
    - Linear interpolation in between.
    - Presence of VT or AFib caps the score at 0.3.
    """
    burden_pct = analysis.arrhythmia_burden_pct
    score = max(0.0, 1.0 - burden_pct / 10.0)

    # Hard cap for dangerous rhythms.
    if analysis.has_vtach or analysis.has_afib:
        score = min(score, 0.3)

    return score


def _hrv_score(analysis: CardiacAnalysis, baseline: BreedBaseline) -> float:
    """Score HRV (SDNN) relative to breed baseline range.

    - SDNN within ``[sdnn_min, sdnn_max]`` -> 1.0
    - SDNN below sdnn_min -> proportional reduction (linear to 0 at half
      the minimum).
    - SDNN above sdnn_max -> slight reduction (generally not concerning;
      cap at 0.9 if very high).
    """
    sdnn = analysis.hrv_metrics.sdnn
    sdnn_min = baseline.sdnn_min
    sdnn_max = baseline.sdnn_max

    if sdnn_min <= sdnn <= sdnn_max:
        return 1.0

    if sdnn < sdnn_min:
        # Linear decrease.  Reaches 0 when SDNN is at half the minimum.
        floor = sdnn_min * 0.5
        if floor <= 0:
            return 0.0
        return max(0.0, (sdnn - floor) / (sdnn_min - floor))

    # Above range — mild reduction.
    overshoot = sdnn - sdnn_max
    return max(0.5, 1.0 - overshoot / (sdnn_max * 2))


def _heart_rate_score(analysis: CardiacAnalysis, baseline: BreedBaseline) -> float:
    """Score resting heart rate relative to breed range.

    - Mean HR within ``[hr_min, hr_max]`` -> 1.0
    - Outside range -> proportional reduction (linear).
    """
    mean_hr = analysis.mean_hr_bpm
    hr_min = float(baseline.hr_min)
    hr_max = float(baseline.hr_max)

    if hr_min <= mean_hr <= hr_max:
        return 1.0

    if mean_hr < hr_min:
        # Deviation below.  0 at 50 % of minimum.
        floor = hr_min * 0.5
        if floor <= 0:
            return 0.0
        return max(0.0, (mean_hr - floor) / (hr_min - floor))

    # Deviation above.
    ceiling = hr_max * 1.5
    return max(0.0, (ceiling - mean_hr) / (ceiling - hr_max))


# ---------------------------------------------------------------------------
# Risk flags & risk level
# ---------------------------------------------------------------------------


def _compute_risk_flags(
    analysis: CardiacAnalysis,
    baseline: BreedBaseline,
    overall_score: float,
) -> list[str]:
    """Build a list of human-readable risk flag strings."""
    flags: list[str] = []

    if analysis.has_vtach:
        flags.append("Ventricular tachycardia detected")

    burden = analysis.arrhythmia_burden_pct
    if burden > 10.0:
        flags.append(f"High arrhythmia burden ({burden:.1f}%)")

    if overall_score < 40.0:
        flags.append(f"Low overall cardiac score ({overall_score:.1f})")

    # Count PVCs.
    pvc_count = sum(
        1 for e in analysis.arrhythmia_events if e.arrhythmia_type == ArrhythmiaType.PVC
    )
    if pvc_count > 50:
        flags.append(f"Frequent PVCs ({pvc_count})")

    if analysis.has_afib:
        flags.append("Atrial fibrillation detected")

    if analysis.hrv_metrics.sdnn < baseline.sdnn_min:
        flags.append(
            f"HRV below breed norm (SDNN {analysis.hrv_metrics.sdnn:.1f} ms "
            f"vs min {baseline.sdnn_min:.1f} ms)"
        )

    return flags


def _determine_risk_level(
    analysis: CardiacAnalysis,
    risk_flags: list[str],
    overall_score: float,
) -> RiskLevel:
    """Assign a risk level based on flags and overall score.

    RED triggers:
    - VT detected
    - Arrhythmia burden > 10 %
    - Overall score < 40

    YELLOW triggers:
    - Frequent PVCs (> 50)
    - AFib detected
    - HRV below breed norm
    - Overall score 40-70

    GREEN:
    - Score >= 70 with no RED/YELLOW triggers.
    """
    # RED conditions.
    if analysis.has_vtach:
        return RiskLevel.RED
    if analysis.arrhythmia_burden_pct > 10.0:
        return RiskLevel.RED
    if overall_score < 40.0:
        return RiskLevel.RED

    # YELLOW conditions.
    pvc_count = sum(
        1 for e in analysis.arrhythmia_events if e.arrhythmia_type == ArrhythmiaType.PVC
    )
    if pvc_count > 50:
        return RiskLevel.YELLOW
    if analysis.has_afib:
        return RiskLevel.YELLOW

    # Check for HRV flag.
    if any("HRV below breed norm" in f for f in risk_flags):
        return RiskLevel.YELLOW

    if overall_score < 70.0:
        return RiskLevel.YELLOW

    return RiskLevel.GREEN


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    """Clamp a value to the range ``[lo, hi]``."""
    return max(lo, min(hi, value))
