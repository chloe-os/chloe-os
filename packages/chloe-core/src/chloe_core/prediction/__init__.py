"""MHC binding prediction with pluggable backends."""

from chloe_core.prediction.base import MHCPredictor, generate_peptides, predict_variants
from chloe_core.prediction.mhcflurry_backend import MHCflurryPredictor

__all__ = [
    "MHCPredictor",
    "MHCflurryPredictor",
    "generate_peptides",
    "predict_variants",
]
