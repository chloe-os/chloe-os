"""Shared data models for the Chloe OS pipeline.

These dataclasses define the inputs and outputs of each pipeline stage,
providing a typed contract between stages.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class VariantConsequence(str, Enum):
    """Classification of a variant's effect on protein sequence."""

    MISSENSE = "missense_variant"
    FRAMESHIFT = "frameshift_variant"
    NONSENSE = "stop_gained"
    INFRAME_INSERTION = "inframe_insertion"
    INFRAME_DELETION = "inframe_deletion"
    SPLICE_SITE = "splice_donor_variant"
    SYNONYMOUS = "synonymous_variant"
    OTHER = "other"


class Species(str, Enum):
    """Supported species."""

    CANINE = "canine"


@dataclass(frozen=True)
class GenomeReference:
    """Reference genome configuration."""

    species: Species
    assembly: str  # e.g. "CanFam_GSD", "CanFam3.1"
    source: str  # URL or local path


# Default reference genomes
CANINE_REFERENCES: dict[str, GenomeReference] = {
    "CanFam_GSD": GenomeReference(
        species=Species.CANINE,
        assembly="CanFam_GSD",
        source="https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011100685.1/",
    ),
    "CanFam3.1": GenomeReference(
        species=Species.CANINE,
        assembly="CanFam3.1",
        source="https://www.ncbi.nlm.nih.gov/assembly/317138",
    ),
}


@dataclass
class Variant:
    """A single somatic variant from VCF parsing."""

    chrom: str
    pos: int
    ref: str
    alt: str
    quality: float
    read_depth: int | None = None
    allele_frequency: float | None = None
    filter_status: str = "PASS"
    info: dict[str, Any] = field(default_factory=dict)

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


@dataclass
class VariantSet:
    """Output of Stage 1: parsed and filtered somatic variants.

    Contains all variants that passed quality filters from the tumor VCF,
    along with summary statistics.
    """

    variants: list[Variant]
    total_variants_raw: int
    total_variants_filtered: int
    tumor_vcf_path: str
    normal_vcf_path: str | None = None
    species: Species = Species.CANINE
    assembly: str = "CanFam_GSD"

    @property
    def somatic_count(self) -> int:
        return len(self.variants)


@dataclass
class AnnotatedVariant:
    """A variant annotated with gene and protein impact information."""

    variant: Variant
    gene_symbol: str | None = None
    gene_id: str | None = None
    transcript_id: str | None = None
    consequence: VariantConsequence = VariantConsequence.OTHER
    protein_change: str | None = None  # e.g. "V600E"
    codon_change: str | None = None  # e.g. "GTG/GAG"
    wildtype_peptide: str | None = None  # surrounding wildtype peptide sequence
    mutant_peptide: str | None = None  # surrounding mutant peptide sequence
    impact: str | None = None  # HIGH, MODERATE, LOW, MODIFIER

    @property
    def has_protein_change(self) -> bool:
        return self.mutant_peptide is not None and self.wildtype_peptide is not None


@dataclass
class AnnotatedVariantSet:
    """Output of Stage 2: variants annotated with gene/protein impact.

    Each variant is enriched with the affected gene, protein change,
    and surrounding peptide sequences for MHC prediction.
    """

    annotated_variants: list[AnnotatedVariant]
    total_annotated: int
    total_with_protein_change: int
    vep_version: str | None = None
    assembly: str = "CanFam_GSD"

    @property
    def protein_changing_variants(self) -> list[AnnotatedVariant]:
        return [v for v in self.annotated_variants if v.has_protein_change]


@dataclass
class BindingPrediction:
    """MHC binding prediction for a single peptide-allele pair."""

    peptide: str
    allele: str  # e.g. "DLA-88*001:01"
    ic50: float  # Binding affinity in nM (lower = stronger binding)
    percentile_rank: float | None = None  # Percentile rank (lower = stronger)
    is_binder: bool = False  # ic50 < 500nM
    is_strong_binder: bool = False  # ic50 < 50nM
    predictor: str = ""  # e.g. "mhcflurry", "netmhcpan"

    @property
    def binding_category(self) -> str:
        if self.is_strong_binder:
            return "strong"
        if self.is_binder:
            return "weak"
        return "non-binder"


@dataclass
class VariantPrediction:
    """MHC binding predictions for all peptides derived from a single variant."""

    annotated_variant: AnnotatedVariant
    mutant_predictions: list[BindingPrediction]
    wildtype_predictions: list[BindingPrediction]

    @property
    def best_mutant_binding(self) -> BindingPrediction | None:
        if not self.mutant_predictions:
            return None
        return min(self.mutant_predictions, key=lambda p: p.ic50)

    @property
    def best_wildtype_binding(self) -> BindingPrediction | None:
        if not self.wildtype_predictions:
            return None
        return min(self.wildtype_predictions, key=lambda p: p.ic50)

    @property
    def agretopicity(self) -> float | None:
        """Ratio of mutant to wildtype binding. Values < 1 mean mutant binds better."""
        best_mut = self.best_mutant_binding
        best_wt = self.best_wildtype_binding
        if best_mut is None or best_wt is None or best_wt.ic50 == 0:
            return None
        return best_mut.ic50 / best_wt.ic50


@dataclass
class PredictionResults:
    """Output of Stage 3: MHC binding predictions for all variants.

    Contains per-variant binding predictions for both mutant and wildtype peptides.
    """

    variant_predictions: list[VariantPrediction]
    alleles_used: list[str]
    predictor_name: str
    predictor_version: str | None = None
    peptide_lengths: list[int] = field(default_factory=lambda: [8, 9, 10, 11])

    @property
    def total_binders(self) -> int:
        return sum(
            1
            for vp in self.variant_predictions
            if vp.best_mutant_binding and vp.best_mutant_binding.is_binder
        )


@dataclass
class NeoantigenCandidate:
    """A scored and ranked neoantigen candidate."""

    variant_prediction: VariantPrediction
    composite_score: float
    rank: int

    # Individual scoring components
    binding_score: float = 0.0
    agretopicity_score: float = 0.0
    vaf_score: float = 0.0
    consequence_score: float = 0.0
    expression_score: float | None = None

    @property
    def gene(self) -> str | None:
        return self.variant_prediction.annotated_variant.gene_symbol

    @property
    def protein_change(self) -> str | None:
        return self.variant_prediction.annotated_variant.protein_change

    @property
    def best_peptide(self) -> str | None:
        best = self.variant_prediction.best_mutant_binding
        return best.peptide if best else None

    @property
    def best_allele(self) -> str | None:
        best = self.variant_prediction.best_mutant_binding
        return best.allele if best else None

    @property
    def best_ic50(self) -> float | None:
        best = self.variant_prediction.best_mutant_binding
        return best.ic50 if best else None


@dataclass
class RankedNeoantigens:
    """Output of Stage 4: scored and ranked neoantigen candidates.

    Candidates are ordered by composite score (highest first).
    """

    candidates: list[NeoantigenCandidate]
    scoring_weights: dict[str, float] = field(default_factory=dict)
    total_variants_evaluated: int = 0
    total_binders: int = 0

    @property
    def top_candidates(self) -> list[NeoantigenCandidate]:
        """Return candidates with composite_score > 0."""
        return [c for c in self.candidates if c.composite_score > 0]


@dataclass
class PipelineConfig:
    """Configuration for the full pipeline run."""

    species: Species = Species.CANINE
    assembly: str = "CanFam_GSD"
    breed: str | None = None

    # Variant filtering
    min_quality: float = 20.0
    min_read_depth: int = 10
    min_allele_frequency: float = 0.05

    # MHC prediction
    predictor: str = "mhcflurry"  # "mhcflurry" or "netmhcpan"
    alleles: list[str] | None = None  # None = use breed defaults
    peptide_lengths: list[int] = field(default_factory=lambda: [8, 9, 10, 11])

    # Neoantigen ranking
    binding_threshold: float = 500.0  # IC50 in nM
    strong_binding_threshold: float = 50.0
    top_n: int = 20

    # Report
    output_format: str = "html"
    include_ai_interpretation: bool = False
    include_construct_design: bool = False

    # AI interpretation
    llm_provider: str = "anthropic"  # "anthropic", "openai", "ollama"
    llm_model: str | None = None


# Common DLA alleles by breed group
# DLA-88 is the most polymorphic MHC class I gene in dogs
BREED_DLA_ALLELES: dict[str, list[str]] = {
    "default": [
        "DLA-88*001:01",
        "DLA-88*002:01",
        "DLA-88*003:01",
        "DLA-88*004:01",
        "DLA-88*005:01",
        "DLA-88*006:01",
        "DLA-88*007:01",
        "DLA-88*008:01",
    ],
    "staffy": [
        "DLA-88*001:01",
        "DLA-88*002:01",
        "DLA-88*004:01",
        "DLA-88*008:01",
    ],
    "golden_retriever": [
        "DLA-88*001:01",
        "DLA-88*003:01",
        "DLA-88*005:01",
        "DLA-88*007:01",
    ],
    "labrador": [
        "DLA-88*001:01",
        "DLA-88*002:01",
        "DLA-88*003:01",
        "DLA-88*006:01",
    ],
    "german_shepherd": [
        "DLA-88*001:01",
        "DLA-88*002:01",
        "DLA-88*005:01",
        "DLA-88*007:01",
    ],
    "boxer": [
        "DLA-88*001:01",
        "DLA-88*003:01",
        "DLA-88*004:01",
        "DLA-88*006:01",
    ],
}
