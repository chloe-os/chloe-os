# MHC Binding

## How the Immune System "Sees" Cancer

For the immune system to attack a cancer cell, it first needs to **see** the cancer. This happens through a molecular display system called **MHC** (Major Histocompatibility Complex).

In dogs, MHC molecules are called **DLA** (Dog Leukocyte Antigen).

## How MHC/DLA Works

1. Inside every cell, proteins are constantly being broken down into small fragments called **peptides** (typically 8-11 amino acids long)
2. **MHC/DLA molecules** grab these peptides and carry them to the cell surface
3. **T cells** (immune cells that kill abnormal cells) patrol the body, checking these displayed peptides
4. If a T cell recognizes a peptide as foreign or abnormal, it triggers an immune response to destroy the cell

## Why MHC Binding Matters for Vaccines

Not every peptide can bind to MHC molecules. The binding site on MHC has a specific shape, and only peptides with the right chemical properties will fit. This is why Chloe OS predicts MHC binding — if a mutant peptide can't bind to your dog's DLA molecules, the immune system will never see it.

## What Chloe OS Measures

- **IC50** — How strongly a peptide binds to a DLA molecule, measured in nanomolar (nM). Lower numbers mean stronger binding.
  - **< 50 nM** = Strong binder (very likely to be displayed)
  - **< 500 nM** = Weak binder (may be displayed)
  - **> 500 nM** = Non-binder (unlikely to be displayed)

- **Agretopicity** — Does the mutant peptide bind *better* than the normal version? If the mutant peptide binds more strongly, it's more likely to be recognized as new and dangerous by the immune system.

## DLA Diversity

Just like humans have different blood types, dogs have different DLA types. Different breeds tend to have different DLA alleles (versions). Chloe OS uses breed-specific DLA profiles when available to make more accurate predictions for your dog.
