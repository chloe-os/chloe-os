# Chloe OS

**Open Source Intelligence for Lifelong Pet Health**

Chloe OS is an open-source pipeline that helps identify personalized cancer vaccine candidates for pets. It takes tumor DNA sequencing data and produces an actionable report you can share with your veterinarian.

!!! warning "Important"
    Chloe OS is a research tool. It does not provide veterinary medical advice. Always consult a qualified veterinarian before making treatment decisions.

## How It Works

Your pet's tumor has unique mutations — "typos" in its DNA that make cancer cells different from healthy cells. Chloe OS reads those mutations and identifies which ones might make good targets for the immune system to attack.

```
Tumor Sequencing Data → Find Mutations → Predict Immune Targets → Rank Candidates → Report
```

## Getting Started

New to Chloe OS? Start here:

- **[Getting Started Guide](guide/getting-started.md)** — Install and run your first analysis
- **[Getting Sequencing Done](guide/getting-sequencing.md)** — How to get your pet's tumor sequenced
- **[Understanding Your Results](guide/understanding-results.md)** — What the report means

## The Science

Want to understand what's happening under the hood?

- **[DNA & Mutations](science/dna-and-mutations.md)** — Why cancer happens at the DNA level
- **[Neoantigens](science/neoantigens.md)** — What makes a good vaccine target
- **[mRNA Vaccines](science/mrna-vaccines.md)** — How personalized vaccines work

## For Developers

- **[Architecture](developer/architecture.md)** — How the pipeline is built
- **[API Reference](developer/api-reference.md)** — Use Chloe OS as a Python library
- **[Contributing](developer/contributing.md)** — Help improve Chloe OS
