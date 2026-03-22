# Getting Started

This guide walks you through installing Chloe OS and running your first analysis.

## Prerequisites

- **Python 3.11 or newer** — [Download Python](https://www.python.org/downloads/)
- **Docker** — Required for variant annotation. [Install Docker](https://docs.docker.com/get-docker/)
- **A VCF file** — This is the file your sequencing lab gives you containing your pet's tumor mutations

!!! tip "Don't have a VCF file yet?"
    See [Getting Sequencing Done](getting-sequencing.md) for how to get your pet's tumor sequenced. You can also try Chloe OS with our example data to see how it works.

## Installation

### Option 1: Using uv (recommended)

```bash
git clone https://github.com/chloe-os/chloe-os.git
cd chloe-os
uv sync
```

### Option 2: Using pip

```bash
git clone https://github.com/chloe-os/chloe-os.git
cd chloe-os
pip install -e packages/chloe-core -e packages/chloe-cli
```

## Your First Analysis

### Guided Mode (recommended)

The guided mode walks you through every step, explaining what's happening in plain language:

```bash
chloe guided
```

It will ask you for:

1. The path to your VCF file
2. Your dog's breed (optional — helps select the right immune markers)
3. Where to save the report

### Direct Mode

If you prefer a single command:

```bash
chloe run --vcf /path/to/tumor.vcf --breed staffy --output report.html
```

### Try with Example Data

```bash
chloe run --vcf examples/sample_canine_tumor/tumor.vcf --output example_report.html
```

## What Happens Next

After the analysis completes, you'll have an HTML report you can open in any web browser. See [Understanding Your Results](understanding-results.md) for what to do with it.

!!! warning "Important"
    **Do not act on these results without consulting a veterinarian.** The report identifies *candidates* for further investigation — it does not prescribe treatment.
