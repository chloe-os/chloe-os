"""HTML report generation for the Chloe OS neoantigen pipeline.

Stage 7 of the neoantigen identification pipeline. Renders a self-contained
HTML report summarising variant calls, MHC binding predictions, and ranked
neoantigen candidates. The output is a single file with all CSS inlined —
no external dependencies required to view it.
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from chloe_core.models import (
    AnnotatedVariantSet,
    NeoantigenCandidate,
    PipelineConfig,
    PredictionResults,
    RankedNeoantigens,
    VariantSet,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Template directory (sibling to this module file)
# ---------------------------------------------------------------------------
_TEMPLATE_DIR = Path(__file__).resolve().parent
_TEMPLATE_NAME = "template.html"


# ---------------------------------------------------------------------------
# Helpers — prepare template context
# ---------------------------------------------------------------------------


def _binding_category_label(candidate: NeoantigenCandidate) -> str:
    """Return a human-readable binding category for a candidate."""
    best = candidate.variant_prediction.best_mutant_binding
    if best is None:
        return "N/A"
    return best.binding_category.replace("_", " ").title()


def _format_ic50(value: float | None) -> str:
    """Format an IC50 value for display."""
    if value is None:
        return "N/A"
    if value < 1:
        return f"{value:.2f}"
    if value < 100:
        return f"{value:.1f}"
    return f"{value:.0f}"


def _build_candidate_details(candidate: NeoantigenCandidate) -> dict:
    """Flatten a NeoantigenCandidate into a plain dict for template use."""
    av = candidate.variant_prediction.annotated_variant
    best_mut = candidate.variant_prediction.best_mutant_binding
    best_wt = candidate.variant_prediction.best_wildtype_binding

    return {
        "rank": candidate.rank,
        "gene": candidate.gene or "Unknown",
        "protein_change": candidate.protein_change or "N/A",
        "best_peptide": candidate.best_peptide or "N/A",
        "best_allele": candidate.best_allele or "N/A",
        "best_ic50": candidate.best_ic50,
        "best_ic50_fmt": _format_ic50(candidate.best_ic50),
        "binding_category": _binding_category_label(candidate),
        "composite_score": round(candidate.composite_score, 4),
        "binding_score": round(candidate.binding_score, 4),
        "agretopicity_score": round(candidate.agretopicity_score, 4),
        "vaf_score": round(candidate.vaf_score, 4),
        "consequence_score": round(candidate.consequence_score, 4),
        "expression_score": (
            round(candidate.expression_score, 4)
            if candidate.expression_score is not None
            else None
        ),
        "variant_id": av.variant.variant_id,
        "chrom": av.variant.chrom,
        "pos": av.variant.pos,
        "ref": av.variant.ref,
        "alt": av.variant.alt,
        "consequence": av.consequence.value,
        "impact": av.impact or "N/A",
        "transcript_id": av.transcript_id or "N/A",
        "wildtype_peptide": av.wildtype_peptide or "N/A",
        "mutant_peptide": av.mutant_peptide or "N/A",
        "agretopicity": (
            round(candidate.variant_prediction.agretopicity, 4)
            if candidate.variant_prediction.agretopicity is not None
            else None
        ),
        "best_wt_ic50": best_wt.ic50 if best_wt else None,
        "best_wt_ic50_fmt": _format_ic50(best_wt.ic50) if best_wt else "N/A",
        "best_wt_peptide": best_wt.peptide if best_wt else "N/A",
        "allele_frequency": av.variant.allele_frequency,
        "read_depth": av.variant.read_depth,
        "all_mutant_predictions": [
            {
                "peptide": p.peptide,
                "allele": p.allele,
                "ic50": p.ic50,
                "ic50_fmt": _format_ic50(p.ic50),
                "percentile_rank": (
                    round(p.percentile_rank, 2)
                    if p.percentile_rank is not None
                    else None
                ),
                "binding_category": p.binding_category.replace("_", " ").title(),
                "is_strong_binder": p.is_strong_binder,
                "is_binder": p.is_binder,
            }
            for p in candidate.variant_prediction.mutant_predictions
        ],
    }


def _build_template_context(
    ranked: RankedNeoantigens,
    variant_set: VariantSet,
    annotated_set: AnnotatedVariantSet,
    predictions: PredictionResults,
    config: PipelineConfig,
) -> dict:
    """Assemble the full Jinja2 template context dictionary."""
    candidates = [_build_candidate_details(c) for c in ranked.candidates]

    # Count binding categories across all candidates
    strong_count = sum(1 for c in candidates if c["binding_category"] == "Strong")
    weak_count = sum(1 for c in candidates if c["binding_category"] == "Weak")

    return {
        # Metadata
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "chloe_version": "0.1.0",
        # Executive summary numbers
        "total_variants_raw": variant_set.total_variants_raw,
        "total_variants_filtered": variant_set.total_variants_filtered,
        "somatic_count": variant_set.somatic_count,
        "total_annotated": annotated_set.total_annotated,
        "total_with_protein_change": annotated_set.total_with_protein_change,
        "total_binders": predictions.total_binders,
        "num_candidates": len(ranked.candidates),
        "strong_count": strong_count,
        "weak_count": weak_count,
        # Pipeline config
        "species": config.species.value.title(),
        "assembly": config.assembly,
        "breed": config.breed or "Not specified",
        "predictor": config.predictor,
        "peptide_lengths": ", ".join(str(n) for n in config.peptide_lengths),
        "alleles_used": predictions.alleles_used,
        "binding_threshold": config.binding_threshold,
        "strong_binding_threshold": config.strong_binding_threshold,
        "min_quality": config.min_quality,
        "min_read_depth": config.min_read_depth,
        "min_allele_frequency": config.min_allele_frequency,
        "vep_version": annotated_set.vep_version or "N/A",
        "predictor_version": predictions.predictor_version or "N/A",
        "scoring_weights": ranked.scoring_weights,
        # Candidates
        "candidates": candidates,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_report(
    ranked: RankedNeoantigens,
    variant_set: VariantSet,
    annotated_set: AnnotatedVariantSet,
    predictions: PredictionResults,
    config: PipelineConfig,
    output_path: str,
) -> str:
    """Render the HTML neoantigen report and write it to *output_path*.

    Parameters
    ----------
    ranked:
        Scored and ranked neoantigen candidates (Stage 4 output).
    variant_set:
        Parsed somatic variants (Stage 1 output).
    annotated_set:
        Annotated variants with gene/protein impact (Stage 2 output).
    predictions:
        MHC binding predictions (Stage 3 output).
    config:
        Pipeline configuration used for this run.
    output_path:
        Filesystem path where the HTML report will be written.

    Returns
    -------
    str
        The *output_path* that was written to (for convenience in chaining).
    """
    logger.info("Generating HTML report → %s", output_path)

    context = _build_template_context(
        ranked, variant_set, annotated_set, predictions, config
    )

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
        "Report written (%d candidates, %d KB)",
        len(ranked.candidates),
        len(html) // 1024,
    )
    return str(out)
