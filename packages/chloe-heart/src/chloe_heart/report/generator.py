"""HTML report generation for the Chloe Heart cardiac monitoring pipeline.

Stage 6 of the cardiac analysis pipeline. Renders a self-contained HTML report
summarising ECG analysis, arrhythmia detection, HRV metrics, and the composite
cardiac health score. The output is a single file with all CSS inlined — no
external dependencies required to view it.
"""

from __future__ import annotations

import logging
import math
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from jinja2 import Environment, FileSystemLoader

from chloe_heart.models import (  # noqa: TCH001
    ArrhythmiaType,
    CardiacAnalysis,
    CardiacHealthScore,
    HeartConfig,
    ProcessedECG,
    RiskLevel,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Template directory (sibling to this module file)
# ---------------------------------------------------------------------------
_TEMPLATE_DIR = Path(__file__).resolve().parent
_TEMPLATE_NAME = "template.html"


# ---------------------------------------------------------------------------
# ECG waveform SVG generation
# ---------------------------------------------------------------------------


def _generate_ecg_svg(
    processed: ProcessedECG,
    arrhythmia_events: list,
    max_seconds: float = 10.0,
) -> str:
    """Generate an inline SVG of the ECG waveform with annotated R-peaks.

    Parameters
    ----------
    processed:
        Preprocessed ECG data with filtered signal and detected beats.
    arrhythmia_events:
        List of ArrhythmiaEvent objects for marking abnormal beats.
    max_seconds:
        Maximum duration to render (defaults to 10 seconds).

    Returns
    -------
    str
        Complete SVG markup as a string for inline embedding in HTML.
    """
    sample_rate = processed.sample_rate
    duration = min(processed.duration_seconds, max_seconds)
    num_samples = int(duration * sample_rate)
    signal = processed.filtered_signal[:num_samples]

    if len(signal) == 0:
        return (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 300">'
            '<text x="500" y="150" text-anchor="middle" fill="#999" '
            'font-size="14">No ECG data available</text></svg>'
        )

    # SVG dimensions
    svg_width = 1000
    svg_height = 300
    margin_left = 50
    margin_right = 20
    margin_top = 30
    margin_bottom = 40
    plot_width = svg_width - margin_left - margin_right
    plot_height = svg_height - margin_top - margin_bottom

    # Scale signal to plot area
    sig_min = float(np.min(signal))
    sig_max = float(np.max(signal))
    sig_range = sig_max - sig_min if sig_max != sig_min else 1.0

    def x_pos(sample_idx: int) -> float:
        return margin_left + (sample_idx / num_samples) * plot_width

    def y_pos(value: float) -> float:
        # Invert Y axis (SVG Y increases downward)
        normalized = (value - sig_min) / sig_range
        return margin_top + plot_height * (1.0 - normalized)

    # Build abnormal beat times set for quick lookup
    abnormal_times: set[float] = set()
    for evt in arrhythmia_events:
        if evt.arrhythmia_type not in (ArrhythmiaType.NORMAL, ArrhythmiaType.SINUS_ARRHYTHMIA):
            abnormal_times.add(round(evt.time_seconds, 4))

    # Start building SVG
    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_width} {svg_height}" '
        f'style="width:100%;height:auto;background:#fff;'
        f'border:1px solid #e5e7eb;border-radius:8px;">'
    )

    # ECG paper grid background
    # Major grid lines (every 1 second horizontally, every 0.5 mV vertically)
    parts.append('<g class="ecg-grid">')

    # Vertical grid lines (time) — minor every 0.2s, major every 1.0s
    time_step_minor = 0.2
    t = 0.0
    while t <= duration:
        px = margin_left + (t / duration) * plot_width
        is_major = abs(t - round(t)) < 0.01
        if is_major:
            parts.append(
                f'<line x1="{px:.1f}" y1="{margin_top}" '
                f'x2="{px:.1f}" y2="{margin_top + plot_height}" '
                f'stroke="#f4c2c2" stroke-width="0.8" />'
            )
            # Time label
            parts.append(
                f'<text x="{px:.1f}" y="{svg_height - 8}" text-anchor="middle" '
                f'fill="#999" font-size="10" font-family="sans-serif">{t:.0f}s</text>'
            )
        else:
            parts.append(
                f'<line x1="{px:.1f}" y1="{margin_top}" '
                f'x2="{px:.1f}" y2="{margin_top + plot_height}" '
                f'stroke="#fce4e4" stroke-width="0.4" />'
            )
        t += time_step_minor

    # Horizontal grid lines (amplitude)
    num_h_major = 5
    for i in range(num_h_major + 1):
        py = margin_top + (i / num_h_major) * plot_height
        parts.append(
            f'<line x1="{margin_left}" y1="{py:.1f}" x2="{margin_left + plot_width}" y2="{py:.1f}" '
            f'stroke="#f4c2c2" stroke-width="0.8" />'
        )
        # Amplitude label
        amp_val = sig_max - (i / num_h_major) * sig_range
        parts.append(
            f'<text x="{margin_left - 6}" y="{py + 3:.1f}" text-anchor="end" '
            f'fill="#999" font-size="9" font-family="sans-serif">{amp_val:.1f}</text>'
        )

    # Minor horizontal grid
    num_h_minor = num_h_major * 5
    for i in range(num_h_minor + 1):
        if i % 5 == 0:
            continue
        py = margin_top + (i / num_h_minor) * plot_height
        parts.append(
            f'<line x1="{margin_left}" y1="{py:.1f}" x2="{margin_left + plot_width}" y2="{py:.1f}" '
            f'stroke="#fce4e4" stroke-width="0.3" />'
        )

    parts.append("</g>")

    # Plot border
    parts.append(
        f'<rect x="{margin_left}" y="{margin_top}" width="{plot_width}" height="{plot_height}" '
        f'fill="none" stroke="#e5e7eb" stroke-width="1" />'
    )

    # ECG waveform path — downsample if too many points
    max_points = 2000
    step = max(1, num_samples // max_points)
    indices = range(0, num_samples, step)

    path_parts: list[str] = []
    for i, idx in enumerate(indices):
        px = x_pos(idx)
        py = y_pos(float(signal[idx]))
        if i == 0:
            path_parts.append(f"M{px:.1f},{py:.1f}")
        else:
            path_parts.append(f"L{px:.1f},{py:.1f}")

    parts.append(
        f'<path d="{" ".join(path_parts)}" fill="none" stroke="#c0392b" '
        f'stroke-width="1.5" stroke-linejoin="round" stroke-linecap="round" />'
    )

    # R-peak markers
    parts.append('<g class="r-peaks">')
    for qrs in processed.qrs_complexes:
        if qrs.r_peak_time > duration:
            break
        px = x_pos(qrs.r_peak_index) if qrs.r_peak_index < num_samples else x_pos(num_samples - 1)
        if qrs.r_peak_index >= num_samples:
            continue
        py = y_pos(qrs.r_peak_amplitude)

        # Determine if this beat is abnormal
        is_abnormal = round(qrs.r_peak_time, 4) in abnormal_times
        # Also check proximity for floating-point tolerance
        if not is_abnormal:
            for at in abnormal_times:
                if abs(qrs.r_peak_time - at) < 0.05:
                    is_abnormal = True
                    break

        color = "#e74c3c" if is_abnormal else "#27ae60"
        radius = 5 if is_abnormal else 4
        parts.append(
            f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{radius}" '
            f'fill="{color}" stroke="#fff" stroke-width="1.5" />'
        )
    parts.append("</g>")

    # Legend
    legend_y = svg_height - 8
    parts.append(
        f'<circle cx="{margin_left + plot_width - 180}" cy="{legend_y - 3}" r="4" fill="#27ae60" />'
        f'<text x="{margin_left + plot_width - 172}" y="{legend_y}" fill="#666" '
        f'font-size="10" font-family="sans-serif">Normal beat</text>'
    )
    parts.append(
        f'<circle cx="{margin_left + plot_width - 90}" cy="{legend_y - 3}" r="4" fill="#e74c3c" />'
        f'<text x="{margin_left + plot_width - 82}" y="{legend_y}" fill="#666" '
        f'font-size="10" font-family="sans-serif">Abnormal</text>'
    )

    # Y-axis label
    parts.append(
        f'<text x="14" y="{margin_top + plot_height / 2}" text-anchor="middle" '
        f'fill="#999" font-size="10" font-family="sans-serif" '
        f'transform="rotate(-90, 14, {margin_top + plot_height / 2})">Amplitude (mV)</text>'
    )

    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Helpers — prepare template context
# ---------------------------------------------------------------------------


def _risk_color(risk_level: RiskLevel) -> str:
    """Return a CSS colour for the risk level."""
    return {
        RiskLevel.GREEN: "#27ae60",
        RiskLevel.YELLOW: "#f39c12",
        RiskLevel.RED: "#e74c3c",
    }.get(risk_level, "#999")


def _risk_label(risk_level: RiskLevel) -> str:
    """Return a human-readable label for the risk level."""
    return {
        RiskLevel.GREEN: "Low Risk",
        RiskLevel.YELLOW: "Moderate Risk",
        RiskLevel.RED: "High Risk",
    }.get(risk_level, "Unknown")


def _arrhythmia_label(arrhythmia_type: ArrhythmiaType) -> str:
    """Return a human-readable label for the arrhythmia type."""
    return {
        ArrhythmiaType.NORMAL: "Normal Sinus Rhythm",
        ArrhythmiaType.SINUS_ARRHYTHMIA: "Sinus Arrhythmia",
        ArrhythmiaType.PVC: "Premature Ventricular Complex",
        ArrhythmiaType.PAC: "Premature Atrial Complex",
        ArrhythmiaType.AFIB: "Atrial Fibrillation",
        ArrhythmiaType.VT: "Ventricular Tachycardia",
        ArrhythmiaType.BRADYCARDIA: "Bradycardia",
        ArrhythmiaType.TACHYCARDIA: "Tachycardia",
        ArrhythmiaType.SECOND_DEGREE_AVB: "2nd Degree AV Block",
        ArrhythmiaType.OTHER: "Other",
    }.get(arrhythmia_type, str(arrhythmia_type))


def _arrhythmia_severity(arrhythmia_type: ArrhythmiaType) -> str:
    """Return a severity class for color-coding arrhythmia events."""
    high = {ArrhythmiaType.VT, ArrhythmiaType.AFIB, ArrhythmiaType.SECOND_DEGREE_AVB}
    moderate = {
        ArrhythmiaType.PVC,
        ArrhythmiaType.PAC,
        ArrhythmiaType.BRADYCARDIA,
        ArrhythmiaType.TACHYCARDIA,
    }
    if arrhythmia_type in high:
        return "high"
    if arrhythmia_type in moderate:
        return "moderate"
    return "low"


def _hrv_comparison(value: float, low: float, high: float) -> str:
    """Compare a metric value against a reference range."""
    if value < low:
        return "below"
    if value > high:
        return "above"
    return "within"


def _generate_score_gauge_svg(score: float, risk_level: RiskLevel) -> str:
    """Generate an inline SVG circular gauge for the cardiac health score.

    Creates a donut-arc gauge that fills proportional to the 0-100 score,
    coloured by the risk level.
    """
    color = _risk_color(risk_level)
    size = 200
    cx = size / 2
    cy = size / 2
    radius = 80
    stroke_width = 16

    # Arc spans from -225 degrees to +45 degrees (270 degree sweep)
    start_angle = -225
    end_angle = 45
    total_sweep = 270

    # Calculate the filled portion
    fill_fraction = max(0.0, min(1.0, score / 100.0))
    fill_sweep = total_sweep * fill_fraction

    def polar_to_cart(angle_deg: float) -> tuple[float, float]:
        rad = math.radians(angle_deg)
        return cx + radius * math.cos(rad), cy + radius * math.sin(rad)

    # Background arc (full 270 degrees)
    bg_start = polar_to_cart(start_angle)
    bg_end = polar_to_cart(end_angle)

    # Filled arc
    fill_end_angle = start_angle + fill_sweep
    fill_end = polar_to_cart(fill_end_angle)

    # SVG arc flags
    bg_large_arc = 1 if total_sweep > 180 else 0
    fill_large_arc = 1 if fill_sweep > 180 else 0

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" '
        f'style="width:200px;height:200px;">'
    )

    # Background track
    parts.append(
        f'<path d="M{bg_start[0]:.1f},{bg_start[1]:.1f} '
        f'A{radius},{radius} 0 {bg_large_arc},1 {bg_end[0]:.1f},{bg_end[1]:.1f}" '
        f'fill="none" stroke="#e5e7eb" stroke-width="{stroke_width}" stroke-linecap="round" />'
    )

    # Filled arc
    if fill_sweep > 0.5:
        parts.append(
            f'<path d="M{bg_start[0]:.1f},{bg_start[1]:.1f} '
            f'A{radius},{radius} 0 {fill_large_arc},1 {fill_end[0]:.1f},{fill_end[1]:.1f}" '
            f'fill="none" stroke="{color}" stroke-width="{stroke_width}" stroke-linecap="round" />'
        )

    # Score text in center
    parts.append(
        f'<text x="{cx}" y="{cy - 8}" text-anchor="middle" '
        f'fill="{color}" font-size="42" font-weight="700" font-family="sans-serif">'
        f"{score:.0f}</text>"
    )
    parts.append(
        f'<text x="{cx}" y="{cy + 16}" text-anchor="middle" '
        f'fill="#6b7280" font-size="13" font-family="sans-serif">out of 100</text>'
    )

    # Risk label below
    label = _risk_label(risk_level)
    parts.append(
        f'<text x="{cx}" y="{cy + 38}" text-anchor="middle" '
        f'fill="{color}" font-size="14" font-weight="600" font-family="sans-serif">'
        f"{label}</text>"
    )

    parts.append("</svg>")
    return "\n".join(parts)


def _build_template_context(
    processed: ProcessedECG,
    analysis: CardiacAnalysis,
    score: CardiacHealthScore,
    config: HeartConfig,
) -> dict:
    """Assemble the full Jinja2 template context dictionary."""
    baseline = config.baseline

    # ECG SVG
    ecg_svg = _generate_ecg_svg(processed, analysis.arrhythmia_events)

    # Score gauge SVG
    score_gauge_svg = _generate_score_gauge_svg(score.overall_score, score.risk_level)

    # Arrhythmia events table
    arrhythmia_rows = []
    for evt in analysis.arrhythmia_events:
        if evt.arrhythmia_type in (ArrhythmiaType.NORMAL, ArrhythmiaType.SINUS_ARRHYTHMIA):
            continue
        arrhythmia_rows.append(
            {
                "time": f"{evt.time_seconds:.2f}",
                "type": _arrhythmia_label(evt.arrhythmia_type),
                "confidence": f"{evt.confidence * 100:.0f}%",
                "details": evt.details or "--",
                "severity": _arrhythmia_severity(evt.arrhythmia_type),
            }
        )

    # HRV metrics with breed comparisons
    hrv = analysis.hrv_metrics
    hrv_cards = [
        {
            "label": "SDNN",
            "value": f"{hrv.sdnn:.1f}",
            "unit": "ms",
            "comparison": _hrv_comparison(hrv.sdnn, baseline.sdnn_min, baseline.sdnn_max),
            "ref_range": f"{baseline.sdnn_min:.0f}–{baseline.sdnn_max:.0f} ms",
        },
        {
            "label": "RMSSD",
            "value": f"{hrv.rmssd:.1f}",
            "unit": "ms",
            "comparison": "within",  # No specific breed range for RMSSD
            "ref_range": "Breed-specific",
        },
        {
            "label": "pNN50",
            "value": f"{hrv.pnn50:.1f}",
            "unit": "%",
            "comparison": "within",
            "ref_range": "Breed-specific",
        },
        {
            "label": "Mean HR",
            "value": f"{hrv.mean_hr:.0f}",
            "unit": "BPM",
            "comparison": _hrv_comparison(hrv.mean_hr, baseline.hr_min, baseline.hr_max),
            "ref_range": f"{baseline.hr_min}–{baseline.hr_max} BPM",
        },
    ]

    # Score breakdown (5 components)
    score_components = [
        {
            "name": "Rhythm Regularity",
            "score": score.rhythm_score,
            "pct": score.rhythm_score * 100,
            "weight": score.weights.get("rhythm", 0.0),
        },
        {
            "name": "Heart Rate Variability",
            "score": score.hrv_score,
            "pct": score.hrv_score * 100,
            "weight": score.weights.get("hrv", 0.0),
        },
        {
            "name": "Heart Rate",
            "score": score.heart_rate_score,
            "pct": score.heart_rate_score * 100,
            "weight": score.weights.get("heart_rate", 0.0),
        },
        {
            "name": "Trend Stability",
            "score": score.trend_score,
            "pct": score.trend_score * 100,
            "weight": score.weights.get("trend", 0.0),
        },
        {
            "name": "Exercise Tolerance",
            "score": score.exercise_score,
            "pct": score.exercise_score * 100,
            "weight": score.weights.get("exercise", 0.0),
        },
    ]

    # Risk flags
    risk_flags = []
    for flag_text in score.risk_flags:
        # Determine flag level from the text heuristically
        flag_lower = flag_text.lower()
        if any(
            kw in flag_lower
            for kw in (
                "critical",
                "vtach",
                "ventricular tachycardia",
                "afib",
                "atrial fibrillation",
            )
        ):
            level = "red"
        elif any(kw in flag_lower for kw in ("elevated", "pvc", "abnormal", "low hrv", "high")):
            level = "yellow"
        else:
            level = "yellow"
        risk_flags.append({"text": flag_text, "level": level})

    # Rhythm assessment
    abnormal_types = [
        _arrhythmia_label(e.arrhythmia_type)
        for e in analysis.arrhythmia_events
        if e.arrhythmia_type not in (ArrhythmiaType.NORMAL, ArrhythmiaType.SINUS_ARRHYTHMIA)
    ]
    unique_arrhythmias = sorted(set(abnormal_types))
    if not unique_arrhythmias:
        rhythm_assessment = "Normal sinus rhythm"
    else:
        rhythm_assessment = ", ".join(unique_arrhythmias)

    return {
        # Metadata
        "generated_at": datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
        "chloe_version": "0.1.0",
        # Config
        "breed": config.breed or "Not specified",
        "size": config.size.value.title() if config.size else "Not specified",
        # Cardiac health score
        "overall_score": f"{score.overall_score:.0f}",
        "risk_level": score.risk_level.value,
        "risk_label": _risk_label(score.risk_level),
        "risk_color": _risk_color(score.risk_level),
        "score_gauge_svg": score_gauge_svg,
        # Executive summary
        "mean_hr": f"{analysis.mean_hr_bpm:.0f}",
        "min_hr": f"{analysis.min_hr_bpm:.0f}",
        "max_hr": f"{analysis.max_hr_bpm:.0f}",
        "total_beats": analysis.total_beats,
        "duration_seconds": f"{analysis.duration_seconds:.1f}",
        "duration_minutes": f"{analysis.duration_seconds / 60:.1f}",
        "abnormal_beat_count": analysis.abnormal_beat_count,
        "arrhythmia_burden_pct": f"{analysis.arrhythmia_burden_pct:.1f}",
        "rhythm_assessment": rhythm_assessment,
        # ECG waveform
        "ecg_svg": ecg_svg,
        # Arrhythmia events
        "arrhythmia_rows": arrhythmia_rows,
        "has_arrhythmias": len(arrhythmia_rows) > 0,
        # HRV dashboard
        "hrv_cards": hrv_cards,
        "hrv": {
            "lf_power": f"{hrv.lf_power:.1f}" if hrv.lf_power is not None else None,
            "hf_power": f"{hrv.hf_power:.1f}" if hrv.hf_power is not None else None,
            "lf_hf_ratio": f"{hrv.lf_hf_ratio:.2f}" if hrv.lf_hf_ratio is not None else None,
            "total_power": f"{hrv.total_power:.1f}" if hrv.total_power is not None else None,
        },
        # Score breakdown
        "score_components": score_components,
        # Risk flags
        "risk_flags": risk_flags,
        "has_risk_flags": len(risk_flags) > 0,
        # Breed baseline
        "baseline": {
            "hr_min": baseline.hr_min,
            "hr_max": baseline.hr_max,
            "sdnn_min": baseline.sdnn_min,
            "sdnn_max": baseline.sdnn_max,
            "dcm_risk": baseline.dcm_risk,
            "notes": baseline.notes,
            "size": baseline.size.value.title(),
        },
        # Unique arrhythmia types found
        "unique_arrhythmias": unique_arrhythmias,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_cardiac_report(
    processed: ProcessedECG,
    analysis: CardiacAnalysis,
    score: CardiacHealthScore,
    config: HeartConfig,
    output_path: str,
) -> str:
    """Render the HTML cardiac report and write it to *output_path*.

    Parameters
    ----------
    processed:
        Preprocessed ECG data with filtered signal and detected beats
        (Stage 2 output).
    analysis:
        Cardiac analysis results including arrhythmias and HRV metrics
        (Stage 3 output).
    score:
        Composite cardiac health score and component breakdown
        (Stage 4 output).
    config:
        Pipeline configuration including breed/size and output settings.
    output_path:
        Filesystem path where the HTML report will be written.

    Returns
    -------
    str
        The *output_path* that was written to (for convenience in chaining).
    """
    logger.info("Generating cardiac HTML report -> %s", output_path)

    context = _build_template_context(processed, analysis, score, config)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=True,
    )
    template = env.get_template(_TEMPLATE_NAME)
    html = template.render(**context)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")

    logger.info(
        "Cardiac report written (score=%s, risk=%s, %d KB)",
        context["overall_score"],
        context["risk_level"],
        len(html) // 1024,
    )
    return str(out)
