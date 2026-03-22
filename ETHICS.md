# Ethical Framework

Chloe OS is an open-source research tool that helps identify potential neoantigen candidates from tumor sequencing data. It is designed to make the science of personalized cancer vaccines more accessible — but accessibility comes with responsibility.

## This Is a Tool, Not a Treatment

Chloe OS does **not** provide veterinary medical advice, diagnose disease, or deliver therapy. It processes genomic data through established bioinformatics algorithms and presents results in an accessible format. The output is a starting point for further expert review — not a prescription.

## Core Principles

### 1. Always Involve a Veterinarian

No output from Chloe OS should be acted upon without consultation with a qualified veterinarian, ideally a board-certified veterinary oncologist. The pipeline identifies *candidates* — whether those candidates are appropriate for a specific animal requires clinical judgment that software cannot provide.

### 2. Do Not Replace Standard of Care

Chloe OS is not a substitute for proven cancer treatments such as surgery, chemotherapy, or radiation therapy. Personalized neoantigen vaccines are experimental. They should be considered as a *complement* to standard treatment, not a replacement, and only under veterinary supervision.

### 3. Be Transparent About Limitations

- Neoantigen prediction algorithms have significant false positive and false negative rates
- MHC binding prediction for canine DLA alleles is less validated than for human HLA
- A predicted neoantigen may not actually be immunogenic in vivo
- mRNA construct designs produced by this tool have not been clinically validated
- Results from a single animal cannot be generalized

### 4. No False Hope

Cancer is devastating, and the desire to help a sick pet is powerful. We commit to:
- Clearly stating what the tool can and cannot do at every stage
- Not exaggerating the potential of personalized vaccines
- Presenting success rates and limitations honestly
- Acknowledging that most neoantigen candidates will not lead to clinical benefit

### 5. Open Science

All algorithms, scoring methods, and thresholds used by Chloe OS are documented and open source. Anyone can inspect, critique, and improve them. We believe transparency is essential for tools that touch healthcare, even veterinary healthcare.

## Data Privacy

### What We Collect

By default, Chloe OS collects **nothing**. All processing happens locally on your machine. No data is sent to external servers except:
- If you choose to use the AI interpretation feature, your neoantigen candidates (not raw genomic data) are sent to the configured LLM API
- If you choose to submit to the community dataset, only anonymized summary data is shared (see below)

### Community Data Sharing

The community dataset is strictly opt-in. Before any submission:
- You review exactly what will be shared
- No personally identifiable information is included
- No raw sequencing data is shared
- No pet names or owner details are included
- You can request removal of your submission at any time

### What We Never Do

- Sell or share individual genomic data
- Store API keys or credentials
- Track usage or collect analytics
- Require an account to use the tool

## Guidance for Researchers

If you are using Chloe OS in a research context:
- This tool has not been peer-reviewed as a complete system
- Individual components (VEP, NetMHCpan, MHCflurry) have their own validation literature
- The pipeline's integration of these tools and its scoring algorithm should be validated for your specific use case
- If publishing results, clearly state the tool version, parameters, and reference genome used

## Guidance for Pet Owners

1. **Talk to your vet first.** Before sequencing, before running the pipeline, before anything else.
2. **Understand the costs.** Tumor DNA sequencing typically costs $1,000–$5,000. mRNA vaccine synthesis requires a research lab partner. This is not free or cheap.
3. **Manage expectations.** Personalized cancer vaccines are experimental. There are documented cases of success, but also many cases of no response.
4. **Don't skip proven treatments.** If your vet recommends surgery or chemo, those have established track records. An experimental vaccine is an addition, not a replacement.
5. **Find the right partners.** The output of this tool needs expert review. Look for veterinary schools, research hospitals, or labs experienced in mRNA vaccine design.

## Reporting Misuse

If you become aware of Chloe OS being used in ways that could harm animals or mislead pet owners, please report it by opening an issue on GitHub or contacting the maintainers.

## Changes to This Document

This ethical framework is a living document. As the technology and our understanding evolve, we will update it. Significant changes will be noted in the repository's changelog.
