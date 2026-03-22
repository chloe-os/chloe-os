# Sample Canine Tumor Data

This directory contains **synthetic example data** for testing Chloe OS.

**These are NOT real variants from a real dog tumor.** The variant positions, alleles, and annotations are illustrative only, designed to demonstrate the pipeline's functionality.

## Files

- `tumor.vcf` — Synthetic VCF file with 21 example variants (19 PASS + 2 LowQual)

## Usage

```bash
# Run the pipeline with example data
chloe run --vcf examples/sample_canine_tumor/tumor.vcf --output example_report.html

# Or use guided mode
chloe guided
# When prompted for a VCF file, enter: examples/sample_canine_tumor/tumor.vcf
```

## About the Data

The example VCF contains:
- 19 high-quality somatic variants (PASS filter)
- 2 low-quality variants (LowQual filter) that will be filtered out
- Variant allele frequencies ranging from 0.04 to 0.45
- Read depths ranging from 8 to 120
- Distributed across multiple canine chromosomes
