"""Abstract MHC binding prediction interface and shared utilities.

Defines the pluggable predictor contract and the main entry point for
running Stage 3 of the pipeline: predicting which mutant peptides bind
to canine MHC (DLA) molecules.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod

from chloe_core.models import (
    AnnotatedVariant,
    AnnotatedVariantSet,
    BindingPrediction,
    PredictionResults,
    VariantPrediction,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Abstract predictor interface
# ---------------------------------------------------------------------------


class MHCPredictor(ABC):
    """Backend-agnostic interface for MHC-peptide binding prediction.

    Concrete implementations wrap a specific prediction engine
    (e.g. MHCflurry, NetMHCpan) behind this common API so that the
    rest of the pipeline is decoupled from the tool choice.
    """

    @abstractmethod
    def predict(
        self, peptides: list[str], alleles: list[str]
    ) -> list[BindingPrediction]:
        """Predict binding affinities for every peptide-allele combination.

        Parameters
        ----------
        peptides:
            Amino-acid sequences to evaluate (variable lengths accepted).
        alleles:
            MHC allele names, e.g. ``["DLA-88*001:01"]``.

        Returns
        -------
        list[BindingPrediction]
            One entry per peptide-allele pair.
        """

    @abstractmethod
    def name(self) -> str:
        """Short human-readable name of this predictor backend."""

    @abstractmethod
    def version(self) -> str | None:
        """Version string of the underlying prediction engine, or ``None``."""


# ---------------------------------------------------------------------------
# Peptide generation
# ---------------------------------------------------------------------------


def generate_peptides(sequence: str, lengths: list[int]) -> list[str]:
    """Slide a window of each requested length across *sequence*.

    Parameters
    ----------
    sequence:
        Full amino-acid sequence (e.g. a 25-mer window around a mutation).
    lengths:
        K-mer sizes to generate, e.g. ``[8, 9, 10, 11]``.

    Returns
    -------
    list[str]
        All unique peptides, preserving insertion order.  Peptides
        shorter than the smallest requested length are silently skipped.
    """
    seen: set[str] = set()
    peptides: list[str] = []
    for k in lengths:
        if k > len(sequence):
            continue
        for start in range(len(sequence) - k + 1):
            pep = sequence[start : start + k]
            if pep not in seen:
                seen.add(pep)
                peptides.append(pep)
    return peptides


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------


def _predict_for_variant(
    variant: AnnotatedVariant,
    predictor: MHCPredictor,
    alleles: list[str],
    peptide_lengths: list[int],
) -> VariantPrediction | None:
    """Run predictions for a single annotated variant.

    Returns ``None`` when the variant has no usable peptide sequences.
    """
    if not variant.has_protein_change:
        logger.debug(
            "Skipping variant %s — no protein change",
            variant.variant.variant_id,
        )
        return None

    # These are guaranteed non-None by has_protein_change.
    assert variant.mutant_peptide is not None
    assert variant.wildtype_peptide is not None

    mutant_peptides = generate_peptides(variant.mutant_peptide, peptide_lengths)
    wildtype_peptides = generate_peptides(variant.wildtype_peptide, peptide_lengths)

    if not mutant_peptides:
        logger.debug(
            "No peptides generated for variant %s (sequence too short?)",
            variant.variant.variant_id,
        )
        return None

    mutant_predictions = predictor.predict(mutant_peptides, alleles)
    wildtype_predictions = (
        predictor.predict(wildtype_peptides, alleles) if wildtype_peptides else []
    )

    return VariantPrediction(
        annotated_variant=variant,
        mutant_predictions=mutant_predictions,
        wildtype_predictions=wildtype_predictions,
    )


def predict_variants(
    annotated_set: AnnotatedVariantSet,
    predictor: MHCPredictor,
    alleles: list[str],
    peptide_lengths: list[int] | None = None,
) -> PredictionResults:
    """Run MHC binding predictions for every protein-changing variant.

    This is the main entry point for Stage 3 of the pipeline.  It:

    1. Filters the annotated variant set to those with protein changes.
    2. Generates overlapping k-mer peptides from mutant and wildtype
       flanking sequences.
    3. Predicts binding affinity for each peptide against every supplied
       DLA allele using the given *predictor* backend.

    Parameters
    ----------
    annotated_set:
        Output of Stage 2 (annotation).
    predictor:
        Any concrete ``MHCPredictor`` implementation.
    alleles:
        DLA allele names to predict against.
    peptide_lengths:
        K-mer lengths for peptide generation.  Defaults to
        ``[8, 9, 10, 11]`` if not supplied.

    Returns
    -------
    PredictionResults
        Aggregated per-variant predictions with metadata.
    """
    if peptide_lengths is None:
        peptide_lengths = [8, 9, 10, 11]

    protein_changing = annotated_set.protein_changing_variants
    logger.info(
        "Running MHC predictions: %d protein-changing variants, "
        "%d alleles, peptide lengths %s",
        len(protein_changing),
        len(alleles),
        peptide_lengths,
    )

    variant_predictions: list[VariantPrediction] = []
    for idx, variant in enumerate(protein_changing, 1):
        logger.debug(
            "Processing variant %d/%d: %s",
            idx,
            len(protein_changing),
            variant.variant.variant_id,
        )
        result = _predict_for_variant(variant, predictor, alleles, peptide_lengths)
        if result is not None:
            variant_predictions.append(result)

    total_binders = sum(
        1
        for vp in variant_predictions
        if vp.best_mutant_binding and vp.best_mutant_binding.is_binder
    )
    logger.info(
        "Prediction complete: %d variants evaluated, %d with at least one binder",
        len(variant_predictions),
        total_binders,
    )

    return PredictionResults(
        variant_predictions=variant_predictions,
        alleles_used=alleles,
        predictor_name=predictor.name(),
        predictor_version=predictor.version(),
        peptide_lengths=peptide_lengths,
    )
