# Chloe OS

**Open Source Intelligence for Lifelong Pet Health**

Chloe OS is an open-source pipeline that helps identify personalized cancer vaccine candidates for pets. It takes tumor DNA sequencing data and produces an actionable neoantigen report — something a pet owner can take to a veterinary oncologist or research lab.

> **Important:** Chloe OS is a research tool. It does not provide veterinary medical advice. Always consult a qualified veterinarian before making treatment decisions. See [ETHICS.md](ETHICS.md) for our complete ethical framework.

## What It Does

```
Tumor VCF File → Variant Parsing → Annotation → MHC Binding Prediction → Neoantigen Ranking → Report
```

1. **Parses** somatic mutations from tumor sequencing data (VCF files)
2. **Annotates** variants with gene and protein impact using Ensembl VEP
3. **Predicts** which mutant peptides bind to canine MHC (DLA) molecules
4. **Ranks** neoantigen candidates by predicted immunogenicity
5. **Generates** an HTML report with plain-language explanations

## Quick Start

### Prerequisites

- Python 3.11+
- [uv](https://docs.astral.sh/uv/) (recommended) or pip
- Docker (for VEP annotation)

### Installation

```bash
# Clone the repository
git clone https://github.com/chloe-os/chloe-os.git
cd chloe-os

# Install with uv
uv sync

# Or install with pip
pip install -e packages/chloe-core -e packages/chloe-cli
```

### Usage

**Guided mode** (recommended for first-time users):
```bash
chloe guided
```

**Direct mode** (for developers and researchers):
```bash
chloe run --vcf tumor.vcf --normal normal.vcf \
    --species canine --breed staffy \
    --output report.html
```

**As a Python library:**
```python
from chloe_core import Pipeline

pipeline = Pipeline(species="canine")
results = pipeline.run(tumor_vcf="tumor.vcf")
results.to_html("report.html")
```

## Project Structure

```
chloe-os/
├── packages/
│   ├── chloe-core/     # Pipeline library
│   └── chloe-cli/      # Command-line interface
├── docs/               # Educational documentation
├── docker/             # Docker setup for external tools
├── examples/           # Example data and walkthroughs
└── ETHICS.md           # Ethical framework
```

## The Science

Chloe OS is inspired by [recent cases](https://theconversation.com/) where pet owners have used AI tools and DNA sequencing to explore personalized cancer vaccines for their animals. The pipeline integrates established bioinformatics tools:

| Tool | Purpose | License |
|------|---------|---------|
| [Ensembl VEP](https://ensembl.org/vep) | Variant annotation | Apache 2.0 |
| [MHCflurry](https://github.com/openvax/mhcflurry) | MHC binding prediction | Apache 2.0 |
| [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) | MHC binding prediction (canine DLA) | Free academic |
| [ESMFold](https://github.com/facebookresearch/esm) | Protein structure prediction | MIT |

## Currently Supported

- **Species:** Dogs (canine genome — CanFam_GSD / CanFam3.1)
- **Cancer types:** Any solid tumor with available sequencing data
- **Input:** VCF files from tumor DNA sequencing

## Roadmap

- **v0.1** — Core pipeline (variant parsing → neoantigen ranking → report)
- **v0.2** — AI interpretation, mRNA construct design, community data sharing
- **v0.3** — Web interface, cat support, treatment outcome tracking

## Contributing

We welcome contributions from bioinformaticians, veterinary researchers, software engineers, and anyone passionate about pet health. See [docs/developer/contributing.md](docs/developer/contributing.md) for guidelines.

## Disclaimer

Chloe OS is an experimental research tool. It has not been validated in clinical trials. Neoantigen predictions may contain false positives and false negatives. mRNA construct designs require expert review. **Always involve a veterinary professional.**

## License

Apache License 2.0 — see [LICENSE](LICENSE) for details.
