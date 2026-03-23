"""Stage 4 — Cardiac health scoring.

Re-exports :func:`score_cardiac_health` from the scorer module.
"""

from chloe_heart.scoring.scorer import score_cardiac_health

__all__ = ["score_cardiac_health"]
