# Architecture

Chloe OS is a modular Python pipeline organized as a uv workspace monorepo.

## Package Structure

```
chloe-os/
├── packages/
│   ├── chloe-core/     # Pipeline library (pip install chloe-core)
│   └── chloe-cli/      # CLI application (pip install chloe-cli)
├── docs/               # MkDocs documentation
├── docker/             # Docker setup for external tools
└── examples/           # Example data
```

### chloe-core

The core library implements the pipeline as a series of independent stages. Each stage:

- Takes a well-defined input dataclass
- Produces a well-defined output dataclass
- Can be run independently
- Has no side effects beyond its output

#### Pipeline Stages

```
VariantSet ──→ AnnotatedVariantSet ──→ PredictionResults ──→ RankedNeoantigens ──→ HTML Report
  Stage 1          Stage 2                 Stage 3              Stage 4             Stage 7
 (variants)      (annotation)           (prediction)          (ranking)           (report)
```

| Stage | Module | Input | Output |
|-------|--------|-------|--------|
| 1 | `chloe_core.variants` | VCF file(s) | `VariantSet` |
| 2 | `chloe_core.annotation` | `VariantSet` | `AnnotatedVariantSet` |
| 3 | `chloe_core.prediction` | `AnnotatedVariantSet` | `PredictionResults` |
| 4 | `chloe_core.ranking` | `PredictionResults` | `RankedNeoantigens` |
| 7 | `chloe_core.report` | All outputs | HTML file |

Stages 5 (AI interpretation) and 6 (mRNA construct design) are planned for v0.2.

#### Data Models

All inter-stage data flows through typed dataclasses defined in `chloe_core.models`. This makes the pipeline:

- **Testable** — Mock any stage input/output
- **Composable** — Use individual stages as a library
- **Type-safe** — Mypy catches interface mismatches

### chloe-cli

The CLI wraps chloe-core with two modes:

- **Guided mode** — Interactive step-by-step walkthrough with plain-language explanations
- **Direct mode** — Standard CLI flags for power users

Built with Typer (CLI framework) and Rich (terminal UI).

## External Tools

The pipeline integrates several external bioinformatics tools:

- **Ensembl VEP** — Variant annotation (via Docker)
- **MHCflurry** — MHC binding prediction (via pip)
- **NetMHCpan** — MHC binding prediction (optional, requires separate license)

The prediction module uses a pluggable backend pattern — implement the `MHCPredictor` abstract class to add new predictors.

## Technology Choices

| Choice | Rationale |
|--------|-----------|
| Python 3.11+ | Bioinformatics ecosystem, AI/ML libraries |
| uv workspaces | Fast dependency management, monorepo support |
| Pydantic/dataclasses | Type-safe data models |
| Typer + Rich | Modern CLI with beautiful output |
| Jinja2 | HTML report templating |
| Docker | Reproducible environments for complex tools |
| MkDocs Material | Documentation with great UX |
