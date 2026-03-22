"""Plain-language explanation helpers for non-technical pet owners.

Each function returns a short, jargon-free explanation of a pipeline
concept so that someone with no biology background can follow along.
"""

from __future__ import annotations


def explain_vcf() -> str:
    """Explain what a VCF file is."""
    return (
        "A VCF file is a standard format that stores the results of DNA sequencing. "
        "It lists every spot where your pet's tumor DNA differs from normal DNA. "
        "Think of it as a detailed list of all the genetic typos found in the tumor."
    )


def explain_variants(count: int) -> str:
    """Explain what variants are and how many were found."""
    return (
        f"We found {count:,} genetic changes (called variants) in the tumor. "
        "Most of these are harmless, but some may cause the tumor to produce "
        "unusual proteins that the immune system can learn to recognise and attack."
    )


def explain_annotation() -> str:
    """Explain variant annotation and VEP."""
    return (
        "We now look up each genetic change in a database to understand which gene "
        "it affects and whether it actually changes a protein. This step is like "
        "checking a dictionary to see if a spelling change alters the meaning of a word."
    )


def explain_mhc_binding() -> str:
    """Explain MHC / DLA binding prediction."""
    return (
        "The immune system uses special molecules called DLA (like a security badge "
        "scanner) to check protein fragments. We predict which tumor protein fragments "
        "fit into your pet's DLA molecules -- the ones that fit well are the best "
        "vaccine targets because the immune system will notice them."
    )


def explain_neoantigens(count: int) -> str:
    """Explain neoantigen candidates."""
    if count == 0:
        return (
            "No strong neoantigen candidates were identified in this run. "
            "This does not necessarily mean a vaccine is impossible -- different "
            "prediction settings or additional sequencing data may reveal targets."
        )
    return (
        f"We identified {count:,} neoantigen candidates -- tumor-specific protein "
        "fragments that your pet's immune system is likely to recognise. These are "
        "ranked by how strongly they bind to DLA molecules and how different they "
        "are from normal proteins."
    )


def explain_next_steps() -> str:
    """Explain what happens after the pipeline finishes."""
    return (
        "The report lists the top vaccine targets found in your pet's tumor. "
        "A veterinary oncologist can review these candidates and discuss whether "
        "a personalised vaccine is appropriate. This analysis is a research tool "
        "and is not a substitute for professional veterinary advice."
    )
