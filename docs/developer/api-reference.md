# API Reference

Chloe OS can be used as a Python library for custom pipelines and integrations.

## Quick Start

```python
from chloe_core.variants import load_variants
from chloe_core.annotation import annotate_variants
from chloe_core.prediction import predict_variants, MHCflurryPredictor
from chloe_core.ranking import rank_neoantigens
from chloe_core.report import generate_report
from chloe_core.models import PipelineConfig

# Configure
config = PipelineConfig(species="canine", breed="staffy")

# Run pipeline stages
variant_set = load_variants("tumor.vcf", normal_vcf="normal.vcf", config=config)
annotated = annotate_variants(variant_set)
predictor = MHCflurryPredictor()
predictions = predict_variants(annotated, predictor, alleles=config.alleles)
ranked = rank_neoantigens(predictions, top_n=20)

# Generate report
generate_report(ranked, variant_set, annotated, predictions, config, "report.html")
```

## Core Modules

### `chloe_core.variants`

```python
load_variants(
    tumor_vcf: str,
    normal_vcf: str | None = None,
    config: PipelineConfig | None = None,
) -> VariantSet
```

Parse and filter somatic variants from VCF files.

### `chloe_core.annotation`

```python
annotate_variants(
    variant_set: VariantSet,
    use_docker: bool = True,
    vep_cache_dir: str | None = None,
) -> AnnotatedVariantSet
```

Annotate variants with gene/protein impact via Ensembl VEP.

### `chloe_core.prediction`

```python
predict_variants(
    annotated_set: AnnotatedVariantSet,
    predictor: MHCPredictor,
    alleles: list[str],
    peptide_lengths: list[int] = [8, 9, 10, 11],
) -> PredictionResults
```

Predict MHC binding for all protein-changing variants.

### `chloe_core.ranking`

```python
rank_neoantigens(
    predictions: PredictionResults,
    weights: dict[str, float] | None = None,
    top_n: int = 20,
) -> RankedNeoantigens
```

Score and rank neoantigen candidates.

### `chloe_core.report`

```python
generate_report(
    ranked: RankedNeoantigens,
    variant_set: VariantSet,
    annotated_set: AnnotatedVariantSet,
    predictions: PredictionResults,
    config: PipelineConfig,
    output_path: str,
) -> str
```

Generate a self-contained HTML report.

## Data Models

All data models are defined in `chloe_core.models`. See the source code for complete field documentation.

Key types: `Variant`, `VariantSet`, `AnnotatedVariant`, `AnnotatedVariantSet`, `BindingPrediction`, `VariantPrediction`, `PredictionResults`, `NeoantigenCandidate`, `RankedNeoantigens`, `PipelineConfig`.
