"""VCF variant parsing and filtering for the Chloe OS pipeline.

Stage 1 of the neoantigen identification pipeline. Parses VCF files
containing tumor mutations, applies quality filters, and identifies
somatic variants by comparing tumor and normal samples.

Supports two modes:
  - Tumor-only: A single VCF from a tumor sample, assumed pre-filtered.
  - Paired: Tumor + matched normal VCFs; somatic variants are those
    present in the tumor but absent from the normal sample.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from cyvcf2 import VCF

from chloe_core.models import PipelineConfig, Species, Variant, VariantSet

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default filter thresholds (used when no PipelineConfig is supplied)
# ---------------------------------------------------------------------------
DEFAULT_MIN_QUALITY: float = 20.0
DEFAULT_MIN_DEPTH: int = 10
DEFAULT_MIN_AF: float = 0.05


# ---------------------------------------------------------------------------
# VCF Parsing
# ---------------------------------------------------------------------------


def parse_vcf(path: str) -> list[Variant]:
    """Parse a VCF file and return a list of :class:`Variant` objects.

    Each VCF record may contain multiple ALT alleles; this function
    creates one ``Variant`` per ALT allele so downstream code never
    needs to handle multi-allelic records.

    Parameters
    ----------
    path:
        Filesystem path to the VCF (plain-text or bgzipped).

    Returns
    -------
    list[Variant]
        Parsed variants.  An empty list is returned when the VCF
        contains no records.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file cannot be opened as a valid VCF.
    """
    vcf_path = Path(path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    try:
        vcf_reader = VCF(str(vcf_path))
    except Exception as exc:
        raise ValueError(f"Failed to open VCF file {path}: {exc}") from exc

    variants: list[Variant] = []

    for record in vcf_reader:
        quality = record.QUAL if record.QUAL is not None else 0.0

        # Determine filter status — cyvcf2 returns None when PASS
        filter_status = "PASS" if record.FILTER is None else str(record.FILTER)

        # Extract read depth (DP) — try INFO first, then FORMAT
        read_depth = _extract_depth(record)

        # Extract allele frequency (AF) — try INFO first, then FORMAT
        allele_frequency = _extract_af(record)

        # Collect INFO fields into a plain dict for downstream stages
        info = _extract_info(record)

        # Handle multi-allelic sites: one Variant per ALT allele
        alt_alleles: list[str] = record.ALT if record.ALT else []
        if not alt_alleles:
            # No ALT allele — skip (monomorphic reference site)
            continue

        for alt in alt_alleles:
            variants.append(
                Variant(
                    chrom=str(record.CHROM),
                    pos=int(record.POS),
                    ref=str(record.REF),
                    alt=str(alt),
                    quality=float(quality),
                    read_depth=read_depth,
                    allele_frequency=allele_frequency,
                    filter_status=filter_status,
                    info=info,
                )
            )

    vcf_reader.close()
    logger.info("Parsed %d variant(s) from %s", len(variants), path)
    return variants


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------


def filter_variants(
    variants: list[Variant],
    min_quality: float = DEFAULT_MIN_QUALITY,
    min_depth: int = DEFAULT_MIN_DEPTH,
    min_af: float = DEFAULT_MIN_AF,
) -> list[Variant]:
    """Apply quality-based filters to a list of variants.

    Variants are **kept** when they meet *all* of the following criteria:

    * ``quality >= min_quality``
    * ``read_depth >= min_depth`` (skipped if read_depth is unknown)
    * ``allele_frequency >= min_af`` (skipped if AF is unknown)
    * ``filter_status`` is ``"PASS"`` or ``"."``

    Parameters
    ----------
    variants:
        Input variant list (not modified in place).
    min_quality:
        Minimum QUAL score.
    min_depth:
        Minimum read depth (DP).
    min_af:
        Minimum variant allele frequency.

    Returns
    -------
    list[Variant]
        Variants that passed all filters.
    """
    filtered: list[Variant] = []

    for v in variants:
        # Quality gate
        if v.quality < min_quality:
            continue

        # Read depth gate (None = unknown, so we let it pass)
        if v.read_depth is not None and v.read_depth < min_depth:
            continue

        # Allele frequency gate
        if v.allele_frequency is not None and v.allele_frequency < min_af:
            continue

        # Filter status gate — accept PASS and "." (unfiltered)
        if v.filter_status not in ("PASS", "."):
            continue

        filtered.append(v)

    logger.info(
        "Filtered %d -> %d variant(s) (QUAL>=%.1f, DP>=%d, AF>=%.3f)",
        len(variants),
        len(filtered),
        min_quality,
        min_depth,
        min_af,
    )
    return filtered


# ---------------------------------------------------------------------------
# Somatic identification (paired mode)
# ---------------------------------------------------------------------------


def identify_somatic(
    tumor_variants: list[Variant],
    normal_variants: list[Variant],
) -> list[Variant]:
    """Identify somatic variants present in the tumor but absent from normal.

    Comparison is based on the genomic coordinate key
    ``(chrom, pos, ref, alt)``.  A tumor variant is considered somatic
    if no normal variant shares the same key.

    Parameters
    ----------
    tumor_variants:
        Variants parsed (and optionally filtered) from the tumor VCF.
    normal_variants:
        Variants parsed from the matched normal VCF.

    Returns
    -------
    list[Variant]
        The subset of *tumor_variants* not found in *normal_variants*.
    """
    normal_keys: set[tuple[str, int, str, str]] = {
        (v.chrom, v.pos, v.ref, v.alt) for v in normal_variants
    }

    somatic = [v for v in tumor_variants if (v.chrom, v.pos, v.ref, v.alt) not in normal_keys]

    logger.info(
        "Somatic filtering: %d tumor, %d normal -> %d somatic variant(s)",
        len(tumor_variants),
        len(normal_variants),
        len(somatic),
    )
    return somatic


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def load_variants(
    tumor_vcf: str,
    normal_vcf: str | None = None,
    config: PipelineConfig | None = None,
) -> VariantSet:
    """Parse, filter, and (optionally) intersect tumor/normal VCFs.

    This is the primary public API for Stage 1.

    Parameters
    ----------
    tumor_vcf:
        Path to the tumor VCF file.
    normal_vcf:
        Path to the matched-normal VCF file.  When provided, only
        somatic variants (tumor-specific) are retained.
    config:
        Pipeline configuration.  Filter thresholds and metadata fields
        (species, assembly) are read from this object when supplied;
        otherwise sensible defaults are used.

    Returns
    -------
    VariantSet
        A :class:`VariantSet` containing the final list of somatic
        variants alongside summary statistics.
    """
    # Resolve filter thresholds
    min_quality = config.min_quality if config else DEFAULT_MIN_QUALITY
    min_depth = config.min_read_depth if config else DEFAULT_MIN_DEPTH
    min_af = config.min_allele_frequency if config else DEFAULT_MIN_AF
    species = config.species if config else Species.CANINE
    assembly = config.assembly if config else "CanFam_GSD"

    # Step 1 — Parse tumor VCF
    logger.info("Loading tumor VCF: %s", tumor_vcf)
    tumor_variants = parse_vcf(tumor_vcf)
    total_raw = len(tumor_variants)

    # Step 2 — Quality filtering
    tumor_variants = filter_variants(
        tumor_variants,
        min_quality=min_quality,
        min_depth=min_depth,
        min_af=min_af,
    )

    # Step 3 — Somatic identification (paired mode)
    if normal_vcf is not None:
        logger.info("Loading normal VCF: %s", normal_vcf)
        normal_variants = parse_vcf(normal_vcf)
        tumor_variants = identify_somatic(tumor_variants, normal_variants)

    total_filtered = len(tumor_variants)

    logger.info("Stage 1 complete: %d raw -> %d final variant(s)", total_raw, total_filtered)

    return VariantSet(
        variants=tumor_variants,
        total_variants_raw=total_raw,
        total_variants_filtered=total_filtered,
        tumor_vcf_path=tumor_vcf,
        normal_vcf_path=normal_vcf,
        species=species,
        assembly=assembly,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_depth(record: Any) -> int | None:
    """Extract read depth from a cyvcf2 Variant record.

    Tries the INFO/DP field first (caller-level depth), then falls
    back to the first sample's FORMAT/DP.
    """
    # INFO-level DP
    try:
        dp = record.INFO.get("DP")
        if dp is not None:
            return int(dp)
    except (AttributeError, TypeError):
        pass

    # FORMAT-level DP (per-sample; take first sample)
    try:
        fmt_dp = record.format("DP")
        if fmt_dp is not None and len(fmt_dp) > 0:
            val = int(fmt_dp[0][0])
            if val >= 0:
                return val
    except (KeyError, TypeError, IndexError, ValueError):
        pass

    return None


def _extract_af(record: Any) -> float | None:
    """Extract allele frequency from a cyvcf2 Variant record.

    Tries INFO/AF first, then FORMAT/AF.
    """
    # INFO-level AF
    try:
        af = record.INFO.get("AF")
        if af is not None:
            # AF can be a tuple for multi-allelic; take the first value
            if isinstance(af, (list, tuple)):
                return float(af[0])
            return float(af)
    except (AttributeError, TypeError, ValueError):
        pass

    # FORMAT-level AF
    try:
        fmt_af = record.format("AF")
        if fmt_af is not None and len(fmt_af) > 0:
            val = float(fmt_af[0][0])
            if 0.0 <= val <= 1.0:
                return val
    except (KeyError, TypeError, IndexError, ValueError):
        pass

    return None


def _extract_info(record: Any) -> dict[str, Any]:
    """Build a plain dict from the INFO fields of a cyvcf2 record.

    Only includes fields that are easily serializable (str, int, float,
    bool, list of those).  Skips fields that raise on access.
    """
    info: dict[str, Any] = {}
    try:
        for key, value in record.INFO:
            info[key] = value
    except (TypeError, StopIteration):
        pass
    return info
