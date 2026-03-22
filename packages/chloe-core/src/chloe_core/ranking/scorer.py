"""Neoantigen candidate scoring and ranking (Stage 4).

Scores each variant's neoantigen potential using a weighted composite of four
criteria and returns a ranked list of candidates for vaccine design.

Scoring rationale
-----------------
The composite score combines four orthogonal signals that together predict
whether a mutant peptide will provoke a tumour-specific T-cell response:

1. **Binding affinity** — a peptide that cannot bind MHC class I cannot be
   presented to T cells.  Lower IC50 is better; we normalise against a
   5 000 nM ceiling so strong binders (< 50 nM) score close to 1.0.

2. **Agretopicity** — the mutant peptide should bind *better* than the
   corresponding wildtype peptide (agretopicity < 1).  This ensures the
   immune system sees the neo-epitope as foreign.  We convert the ratio to
   a 0-1 score via ``1 - agretopicity``.

3. **Variant allele frequency (VAF)** — a higher VAF implies the mutation
   is clonal rather than sub-clonal, so a vaccine targeting it will cover
   more tumour cells.  Score is linearly scaled so a VAF of 0.5 already
   reaches the maximum.

4. **Consequence type** — frameshifts produce completely novel peptide
   stretches and score highest; missense mutations are the most common
   actionable class; synonymous and other silent variants score lowest.

Default weights (binding=0.4, agretopicity=0.25, vaf=0.2, consequence=0.15)
reflect relative importance and can be overridden per run.
"""

from __future__ import annotations

import logging
from typing import Any

from chloe_core.models import (
    NeoantigenCandidate,
    PredictionResults,
    RankedNeoantigens,
    VariantConsequence,
    VariantPrediction,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default scoring weights
# ---------------------------------------------------------------------------

DEFAULT_WEIGHTS: dict[str, float] = {
    "binding": 0.40,
    "agretopicity": 0.25,
    "vaf": 0.20,
    "consequence": 0.15,
}

# ---------------------------------------------------------------------------
# Consequence-type priority map
# ---------------------------------------------------------------------------

_CONSEQUENCE_SCORES: dict[VariantConsequence, float] = {
    VariantConsequence.FRAMESHIFT: 1.0,
    VariantConsequence.MISSENSE: 0.8,
    VariantConsequence.NONSENSE: 0.6,
    VariantConsequence.INFRAME_INSERTION: 0.5,
    VariantConsequence.INFRAME_DELETION: 0.5,
}

_DEFAULT_CONSEQUENCE_SCORE: float = 0.2


# ---------------------------------------------------------------------------
# Individual scoring functions (each returns a value in [0, 1])
# ---------------------------------------------------------------------------


def _score_binding(vp: VariantPrediction, threshold: float = 5000.0) -> float:
    """Score the mutant peptide's MHC binding affinity.

    The score is computed as ``max(0, 1 - ic50 / threshold)`` so that an IC50
    of 0 nM yields 1.0 and an IC50 at or above *threshold* yields 0.0.  With
    the default 5 000 nM ceiling, a strong binder at 50 nM scores ~0.99.

    Parameters
    ----------
    vp:
        Variant prediction containing mutant binding data.
    threshold:
        IC50 value (nM) at which the score becomes 0.  Predictions with
        IC50 >= threshold are considered non-immunogenic.

    Returns
    -------
    float
        Binding score in [0, 1].
    """
    best = vp.best_mutant_binding
    if best is None:
        return 0.0
    return max(0.0, 1.0 - best.ic50 / threshold)


def _score_agretopicity(vp: VariantPrediction) -> float:
    """Score differential binding between the mutant and wildtype peptides.

    Agretopicity is the ratio ``mutant_IC50 / wildtype_IC50``.  A value below
    1.0 means the mutant binds MHC better than the wildtype — exactly what we
    want for tumour-specific presentation.  The score is
    ``clamp(1 - agretopicity, 0, 1)``.

    Returns
    -------
    float
        Agretopicity score in [0, 1].  Returns 0.0 when the ratio is
        unavailable (e.g. no wildtype predictions).
    """
    agretopicity = vp.agretopicity
    if agretopicity is None:
        return 0.0
    return max(0.0, min(1.0, 1.0 - agretopicity))


def _score_vaf(vp: VariantPrediction) -> float:
    """Score the variant allele frequency (clonality proxy).

    Higher VAF suggests the variant is present in a larger fraction of tumour
    cells (clonal rather than sub-clonal).  The score saturates at a VAF of
    0.5 (``min(1, vaf * 2)``), since heterozygous somatic mutations in a pure
    tumour rarely exceed 0.5.

    Returns
    -------
    float
        VAF score in [0, 1].  Returns 0.0 when VAF is unknown.
    """
    vaf = vp.annotated_variant.variant.allele_frequency
    if vaf is None:
        return 0.0
    return min(1.0, vaf * 2.0)


def _score_consequence(vp: VariantPrediction) -> float:
    """Score the variant by its functional consequence on the protein.

    Frameshifts produce entirely novel peptide sequences (score 1.0), making
    them the most immunogenic class.  Missense mutations (0.8) are the
    commonest actionable category.  Synonymous and unclassified variants
    receive a low default score (0.2).

    Returns
    -------
    float
        Consequence score in [0, 1].
    """
    consequence = vp.annotated_variant.consequence
    return _CONSEQUENCE_SCORES.get(consequence, _DEFAULT_CONSEQUENCE_SCORE)


# ---------------------------------------------------------------------------
# Composite scoring and ranking
# ---------------------------------------------------------------------------


def _compute_composite_score(
    vp: VariantPrediction,
    weights: dict[str, float],
    binding_threshold: float = 5000.0,
) -> dict[str, Any]:
    """Compute all individual scores and the weighted composite.

    Returns a dict with keys ``binding_score``, ``agretopicity_score``,
    ``vaf_score``, ``consequence_score``, and ``composite_score``.
    """
    binding = _score_binding(vp, threshold=binding_threshold)
    agretopicity = _score_agretopicity(vp)
    vaf = _score_vaf(vp)
    consequence = _score_consequence(vp)

    composite = (
        weights["binding"] * binding
        + weights["agretopicity"] * agretopicity
        + weights["vaf"] * vaf
        + weights["consequence"] * consequence
    )

    return {
        "binding_score": binding,
        "agretopicity_score": agretopicity,
        "vaf_score": vaf,
        "consequence_score": consequence,
        "composite_score": composite,
    }


def rank_neoantigens(
    predictions: PredictionResults,
    weights: dict[str, float] | None = None,
    top_n: int = 20,
) -> RankedNeoantigens:
    """Score and rank neoantigen candidates from Stage 3 binding predictions.

    This is the main entry point for Stage 4 of the pipeline.  Each variant
    prediction is scored across four criteria (binding affinity, agretopicity,
    VAF, and consequence type), combined into a single composite score using
    the supplied weights, and then ranked in descending order.  The top *n*
    candidates are returned.

    Parameters
    ----------
    predictions:
        Output from Stage 3 — per-variant MHC binding predictions for both
        mutant and wildtype peptides.
    weights:
        Scoring-dimension weights.  Keys must be ``"binding"``,
        ``"agretopicity"``, ``"vaf"``, and ``"consequence"``.  If ``None``,
        the default weights are used (0.40, 0.25, 0.20, 0.15).
    top_n:
        Maximum number of candidates to return after ranking.

    Returns
    -------
    RankedNeoantigens
        The ranked candidate list together with metadata about how many
        variants were evaluated and how many produced binders.
    """
    if weights is None:
        weights = dict(DEFAULT_WEIGHTS)

    # Validate weight keys
    expected_keys = {"binding", "agretopicity", "vaf", "consequence"}
    if set(weights.keys()) != expected_keys:
        raise ValueError(f"weights must contain exactly {expected_keys}; got {set(weights.keys())}")

    total_evaluated = len(predictions.variant_predictions)
    total_binders = predictions.total_binders

    logger.info(
        "Scoring %d variant predictions (%d binders)",
        total_evaluated,
        total_binders,
    )

    # Score every variant prediction
    scored: list[tuple[float, dict[str, Any], VariantPrediction]] = []
    for vp in predictions.variant_predictions:
        scores = _compute_composite_score(vp, weights)
        scored.append((scores["composite_score"], scores, vp))

    # Sort descending by composite score
    scored.sort(key=lambda item: item[0], reverse=True)

    # Build ranked candidate objects (top_n only)
    candidates: list[NeoantigenCandidate] = []
    for rank_idx, (composite, scores, vp) in enumerate(scored[:top_n], start=1):
        candidate = NeoantigenCandidate(
            variant_prediction=vp,
            composite_score=composite,
            rank=rank_idx,
            binding_score=scores["binding_score"],
            agretopicity_score=scores["agretopicity_score"],
            vaf_score=scores["vaf_score"],
            consequence_score=scores["consequence_score"],
        )
        candidates.append(candidate)

    logger.info(
        "Ranked %d candidates (top %d of %d evaluated)",
        len(candidates),
        top_n,
        total_evaluated,
    )

    return RankedNeoantigens(
        candidates=candidates,
        scoring_weights=weights,
        total_variants_evaluated=total_evaluated,
        total_binders=total_binders,
    )
