"""Tests for shared data models."""

from chloe_core.models import (
    BREED_DLA_ALLELES,
    AnnotatedVariant,
    AnnotatedVariantSet,
    BindingPrediction,
    NeoantigenCandidate,
    PipelineConfig,
    PredictionResults,
    RankedNeoantigens,
    Species,
    Variant,
    VariantConsequence,
    VariantPrediction,
    VariantSet,
)


def test_variant_id():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    assert v.variant_id == "chr1:100:A>T"


def test_variant_set_somatic_count():
    variants = [
        Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0),
        Variant(chrom="chr2", pos=200, ref="G", alt="C", quality=60.0),
    ]
    vs = VariantSet(
        variants=variants,
        total_variants_raw=10,
        total_variants_filtered=2,
        tumor_vcf_path="tumor.vcf",
    )
    assert vs.somatic_count == 2


def test_annotated_variant_has_protein_change():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    av = AnnotatedVariant(
        variant=v,
        gene_symbol="BRAF",
        consequence=VariantConsequence.MISSENSE,
        wildtype_peptide="VVVVVVVVV",
        mutant_peptide="VVVEVVVVV",
    )
    assert av.has_protein_change is True

    av_no_change = AnnotatedVariant(variant=v)
    assert av_no_change.has_protein_change is False


def test_annotated_variant_set_protein_changing():
    v1 = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    v2 = Variant(chrom="chr2", pos=200, ref="G", alt="C", quality=60.0)
    avs = AnnotatedVariantSet(
        annotated_variants=[
            AnnotatedVariant(
                variant=v1,
                wildtype_peptide="AAA",
                mutant_peptide="AEA",
            ),
            AnnotatedVariant(variant=v2),
        ],
        total_annotated=2,
        total_with_protein_change=1,
    )
    assert len(avs.protein_changing_variants) == 1


def test_binding_prediction_categories():
    strong = BindingPrediction(
        peptide="SIINFEKL", allele="DLA-88*001:01", ic50=10.0,
        is_binder=True, is_strong_binder=True,
    )
    assert strong.binding_category == "strong"

    weak = BindingPrediction(
        peptide="SIINFEKL", allele="DLA-88*001:01", ic50=200.0,
        is_binder=True, is_strong_binder=False,
    )
    assert weak.binding_category == "weak"

    non = BindingPrediction(
        peptide="SIINFEKL", allele="DLA-88*001:01", ic50=5000.0,
        is_binder=False, is_strong_binder=False,
    )
    assert non.binding_category == "non-binder"


def test_variant_prediction_agretopicity():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    av = AnnotatedVariant(variant=v)
    vp = VariantPrediction(
        annotated_variant=av,
        mutant_predictions=[
            BindingPrediction(peptide="MUT", allele="DLA", ic50=50.0, is_binder=True),
        ],
        wildtype_predictions=[
            BindingPrediction(peptide="WT", allele="DLA", ic50=500.0, is_binder=True),
        ],
    )
    assert vp.agretopicity == 0.1  # 50/500 = mutant binds 10x better


def test_prediction_results_total_binders():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    av = AnnotatedVariant(variant=v)
    pr = PredictionResults(
        variant_predictions=[
            VariantPrediction(
                annotated_variant=av,
                mutant_predictions=[
                    BindingPrediction(peptide="P1", allele="DLA", ic50=100.0, is_binder=True),
                ],
                wildtype_predictions=[],
            ),
            VariantPrediction(
                annotated_variant=av,
                mutant_predictions=[
                    BindingPrediction(peptide="P2", allele="DLA", ic50=5000.0, is_binder=False),
                ],
                wildtype_predictions=[],
            ),
        ],
        alleles_used=["DLA-88*001:01"],
        predictor_name="test",
    )
    assert pr.total_binders == 1


def test_neoantigen_candidate_properties():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    av = AnnotatedVariant(
        variant=v,
        gene_symbol="TP53",
        protein_change="R175H",
    )
    vp = VariantPrediction(
        annotated_variant=av,
        mutant_predictions=[
            BindingPrediction(peptide="PEPTIDE", allele="DLA-88*001:01", ic50=25.0, is_binder=True),
        ],
        wildtype_predictions=[],
    )
    nc = NeoantigenCandidate(variant_prediction=vp, composite_score=0.85, rank=1)
    assert nc.gene == "TP53"
    assert nc.protein_change == "R175H"
    assert nc.best_peptide == "PEPTIDE"
    assert nc.best_allele == "DLA-88*001:01"
    assert nc.best_ic50 == 25.0


def test_ranked_neoantigens_top_candidates():
    v = Variant(chrom="chr1", pos=100, ref="A", alt="T", quality=50.0)
    av = AnnotatedVariant(variant=v)
    vp = VariantPrediction(annotated_variant=av, mutant_predictions=[], wildtype_predictions=[])

    rn = RankedNeoantigens(
        candidates=[
            NeoantigenCandidate(variant_prediction=vp, composite_score=0.8, rank=1),
            NeoantigenCandidate(variant_prediction=vp, composite_score=0.5, rank=2),
            NeoantigenCandidate(variant_prediction=vp, composite_score=0.0, rank=3),
        ]
    )
    assert len(rn.top_candidates) == 2  # score > 0 only


def test_pipeline_config_defaults():
    config = PipelineConfig()
    assert config.species == Species.CANINE
    assert config.assembly == "CanFam_GSD"
    assert config.predictor == "mhcflurry"
    assert config.binding_threshold == 500.0
    assert config.top_n == 20


def test_breed_dla_alleles():
    assert "default" in BREED_DLA_ALLELES
    assert "staffy" in BREED_DLA_ALLELES
    assert len(BREED_DLA_ALLELES["default"]) == 8
    assert all(a.startswith("DLA-88*") for a in BREED_DLA_ALLELES["default"])
