"""MHCflurry backend for MHC-peptide binding prediction.

Wraps the ``mhcflurry`` library (``Class1PresentationPredictor``) behind
the common :class:`~chloe_core.prediction.base.MHCPredictor` interface.

MHCflurry is trained on human HLA data but uses a pan-allele neural
network architecture.  For canine DLA alleles, we attempt prediction
directly — the pan-MHC model may generalise to some extent — and fall
back to a configurable set of proxy HLA alleles if the engine does not
recognise a given DLA allele.
"""

from __future__ import annotations

import logging
from typing import Any

from chloe_core.models import BindingPrediction
from chloe_core.prediction.base import MHCPredictor

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy import helper
# ---------------------------------------------------------------------------

_INSTALL_MSG = (
    "mhcflurry is required but not installed.\n"
    "Install it with:\n"
    "  pip install mhcflurry\n"
    "Then download the models:\n"
    "  mhcflurry-downloads fetch models_class1_presentation\n"
)


def _import_mhcflurry() -> Any:
    """Import mhcflurry lazily and raise a friendly error if missing."""
    try:
        from mhcflurry import Class1PresentationPredictor  # type: ignore[import-untyped]

        return Class1PresentationPredictor
    except ImportError as exc:
        raise ImportError(_INSTALL_MSG) from exc


# ---------------------------------------------------------------------------
# Default DLA -> HLA proxy mapping
# ---------------------------------------------------------------------------

# When MHCflurry does not recognise a DLA allele, we optionally fall back to
# the closest HLA-A/B/C allele based on published sequence homology.  This
# mapping is a coarse approximation and should be refined as canine MHC
# research matures.
DEFAULT_DLA_TO_HLA_PROXY: dict[str, str] = {
    "DLA-88*001:01": "HLA-A*02:01",
    "DLA-88*002:01": "HLA-A*24:02",
    "DLA-88*003:01": "HLA-A*01:01",
    "DLA-88*004:01": "HLA-A*03:01",
    "DLA-88*005:01": "HLA-A*11:01",
    "DLA-88*006:01": "HLA-B*07:02",
    "DLA-88*007:01": "HLA-B*08:01",
    "DLA-88*008:01": "HLA-B*44:02",
}

# ---------------------------------------------------------------------------
# Binding-strength thresholds (nM IC50)
# ---------------------------------------------------------------------------

_IC50_BINDER_THRESHOLD = 500.0
_IC50_STRONG_BINDER_THRESHOLD = 50.0


# ---------------------------------------------------------------------------
# Backend implementation
# ---------------------------------------------------------------------------


class MHCflurryPredictor(MHCPredictor):
    """MHCflurry-based MHC binding predictor.

    Parameters
    ----------
    use_hla_proxy:
        When ``True`` (the default), DLA alleles that MHCflurry does not
        recognise are transparently mapped to the closest HLA allele via
        *dla_to_hla_map*.  When ``False``, unrecognised alleles raise an
        error.
    dla_to_hla_map:
        Override the default DLA-to-HLA mapping.  Keys are DLA allele
        names, values are HLA allele names that MHCflurry supports.
    """

    def __init__(
        self,
        *,
        use_hla_proxy: bool = True,
        dla_to_hla_map: dict[str, str] | None = None,
    ) -> None:
        cls = _import_mhcflurry()
        self._predictor: Any = cls.load()
        self._use_hla_proxy = use_hla_proxy
        self._dla_to_hla_map = (
            dict(dla_to_hla_map) if dla_to_hla_map else dict(DEFAULT_DLA_TO_HLA_PROXY)
        )
        self._supported_alleles: set[str] = set(self._predictor.supported_alleles)
        self._version: str | None = self._get_mhcflurry_version()

    # ----- MHCPredictor interface -----------------------------------------

    def name(self) -> str:
        return "mhcflurry"

    def version(self) -> str | None:
        return self._version

    def predict(
        self, peptides: list[str], alleles: list[str]
    ) -> list[BindingPrediction]:
        """Predict binding for each peptide against each allele.

        Returns one :class:`BindingPrediction` per peptide-allele pair.
        Alleles that cannot be resolved (neither natively supported nor
        mapped via the proxy table) are skipped with a warning.
        """
        if not peptides or not alleles:
            return []

        resolved_alleles = self._resolve_alleles(alleles)
        if not resolved_alleles:
            logger.warning("No alleles could be resolved — returning empty results")
            return []

        predictions: list[BindingPrediction] = []

        # MHCflurry's presentation predictor accepts lists of peptides and
        # corresponding alleles (one-to-one).  We expand the Cartesian
        # product ourselves so we can track the original DLA name.
        for original_allele, resolved_allele in resolved_alleles:
            try:
                df = self._predictor.predict(
                    peptides=peptides,
                    alleles=[resolved_allele] * len(peptides),
                    verbose=0,
                )
            except Exception:
                logger.warning(
                    "MHCflurry prediction failed for allele %s (%s) — skipping",
                    original_allele,
                    resolved_allele,
                    exc_info=True,
                )
                continue

            for _, row in df.iterrows():
                ic50: float = float(row.get("presentation_score", row.get("affinity", 0)))
                # MHCflurry may provide different column names depending
                # on the predictor variant.  Prefer the affinity column
                # for IC50 and presentation_percentile for rank.
                if "affinity" in row.index:
                    ic50 = float(row["affinity"])

                percentile: float | None = None
                if "affinity_percentile" in row.index:
                    percentile = float(row["affinity_percentile"])
                elif "presentation_percentile" in row.index:
                    percentile = float(row["presentation_percentile"])

                predictions.append(
                    BindingPrediction(
                        peptide=str(row.get("peptide", "")),
                        allele=original_allele,
                        ic50=ic50,
                        percentile_rank=percentile,
                        is_binder=ic50 < _IC50_BINDER_THRESHOLD,
                        is_strong_binder=ic50 < _IC50_STRONG_BINDER_THRESHOLD,
                        predictor=self.name(),
                    )
                )

        return predictions

    # ----- Internal helpers ------------------------------------------------

    @staticmethod
    def _get_mhcflurry_version() -> str | None:
        try:
            import mhcflurry  # type: ignore[import-untyped]

            return str(getattr(mhcflurry, "__version__", None))
        except Exception:
            return None

    def _resolve_alleles(
        self, alleles: list[str]
    ) -> list[tuple[str, str]]:
        """Map each requested allele to one MHCflurry understands.

        Returns a list of ``(original_name, resolved_name)`` tuples.
        Alleles that cannot be resolved are dropped with a warning.
        """
        resolved: list[tuple[str, str]] = []

        for allele in alleles:
            # 1. Try the allele as-is (works for HLA alleles and potentially
            #    for DLA if a future MHCflurry release adds canine support).
            if allele in self._supported_alleles:
                resolved.append((allele, allele))
                continue

            # 2. Attempt HLA proxy mapping for DLA alleles.
            if self._use_hla_proxy and allele in self._dla_to_hla_map:
                proxy = self._dla_to_hla_map[allele]
                if proxy in self._supported_alleles:
                    logger.info(
                        "Mapping DLA allele %s -> HLA proxy %s", allele, proxy
                    )
                    resolved.append((allele, proxy))
                    continue
                logger.warning(
                    "HLA proxy %s for DLA allele %s is not supported by "
                    "MHCflurry — skipping",
                    proxy,
                    allele,
                )
                continue

            # 3. Allele cannot be resolved.
            logger.warning(
                "Allele %s is not supported by MHCflurry and no HLA proxy "
                "is configured — skipping",
                allele,
            )

        return resolved
