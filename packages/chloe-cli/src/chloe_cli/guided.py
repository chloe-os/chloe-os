"""Interactive guided mode for the Chloe OS pipeline.

Walks the user through each pipeline stage step-by-step with
plain-language explanations, input prompts, and progress spinners.
"""

from __future__ import annotations

from pathlib import Path

import typer
from chloe_core.models import BREED_DLA_ALLELES, PipelineConfig, Species
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Confirm, Prompt
from rich.table import Table

from chloe_cli.explain import (
    explain_annotation,
    explain_mhc_binding,
    explain_neoantigens,
    explain_next_steps,
    explain_variants,
    explain_vcf,
)

console = Console()

# ---------------------------------------------------------------------------
# Styling constants
# ---------------------------------------------------------------------------

STEP_STYLE = "bold cyan"
HEADING_BORDER = "bright_cyan"


def _step_panel(step: int, total: int, title: str, body: str) -> Panel:
    """Create a consistently-styled panel for a pipeline step."""
    return Panel(
        body,
        title=f"[{STEP_STYLE}]Step {step}/{total} — {title}[/]",
        border_style=HEADING_BORDER,
        padding=(1, 2),
    )


def _learn_more(explanation: str) -> None:
    """Optionally show a deeper explanation if the user asks."""
    if Confirm.ask("  [dim]Would you like to learn more?[/]", default=False):
        console.print()
        console.print(Panel(explanation, border_style="dim", padding=(1, 2)))
        console.print()


# ---------------------------------------------------------------------------
# Guided pipeline
# ---------------------------------------------------------------------------

TOTAL_STEPS = 7


def guided() -> None:
    """Interactive step-by-step walkthrough of the neoantigen pipeline."""
    # Lazy imports for heavy libraries.
    from chloe_core.annotation import annotate_variants  # type: ignore[attr-defined]
    from chloe_core.prediction.base import predict_variants
    from chloe_core.ranking import rank_neoantigens  # type: ignore[attr-defined]
    from chloe_core.report import generate_report  # type: ignore[attr-defined]
    from chloe_core.variants import load_variants  # type: ignore[attr-defined]

    # ==================================================================
    # Step 1 — Welcome + disclaimer
    # ==================================================================
    console.print()
    console.print(
        Panel.fit(
            "[bold cyan]Welcome to Chloe OS[/]\n\n"
            "This guided walkthrough will help you identify personalised\n"
            "cancer vaccine targets from your pet's tumor sequencing data.\n\n"
            "[dim]Each step will explain what is happening in plain language\n"
            "and give you the option to learn more.[/]",
            border_style=HEADING_BORDER,
        )
    )
    console.print()
    console.print(
        Panel(
            "[yellow]Disclaimer:[/] Chloe OS is a research tool. Results must be "
            "reviewed by a qualified veterinary oncologist before any clinical "
            "decisions are made. This software does not provide medical advice.",
            border_style="yellow",
        )
    )
    console.print()
    if not Confirm.ask("Ready to begin?", default=True):
        raise typer.Exit()

    # ==================================================================
    # Step 2 — Load VCF file
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            2,
            TOTAL_STEPS,
            "Load VCF File",
            "We need the file containing your pet's tumor sequencing results.\n"
            "This is typically a [bold].vcf[/] or [bold].vcf.gz[/] file produced\n"
            "by a sequencing lab.",
        )
    )
    _learn_more(explain_vcf())

    tumor_vcf = _ask_file_path("Path to the [bold]tumor[/] VCF file")
    normal_vcf: str | None = None
    if Confirm.ask("  Do you also have a matched-normal VCF file?", default=False):
        normal_vcf = _ask_file_path("Path to the [bold]normal[/] VCF file")

    breed = Prompt.ask(
        "  Breed",
        default="staffy",
        choices=list(BREED_DLA_ALLELES.keys()),
    )

    alleles = BREED_DLA_ALLELES.get(breed, BREED_DLA_ALLELES["default"])
    config = PipelineConfig(
        species=Species.CANINE,
        breed=breed,
        alleles=alleles,
    )

    # ==================================================================
    # Step 3 — Variant parsing
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            3,
            TOTAL_STEPS,
            "Parse Variants",
            "Reading the VCF file and filtering for high-quality somatic variants...",
        )
    )

    with console.status("[bold cyan]Parsing variants...", spinner="dots"):
        variant_set = load_variants(
            tumor_vcf=tumor_vcf,
            normal_vcf=normal_vcf,
            species=config.species,
            min_quality=config.min_quality,
            min_read_depth=config.min_read_depth,
            min_allele_frequency=config.min_allele_frequency,
        )

    console.print(
        f"  [green]\u2713[/] Found [bold]{variant_set.somatic_count:,}[/] somatic variants "
        f"(from {variant_set.total_variants_raw:,} total in the file)"
    )
    _learn_more(explain_variants(variant_set.somatic_count))

    # ==================================================================
    # Step 4 — Annotation
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            4,
            TOTAL_STEPS,
            "Annotate Variants",
            "Looking up each variant to see which gene it affects and\n"
            "whether it changes a protein...",
        )
    )

    with console.status("[bold cyan]Annotating variants with VEP...", spinner="dots"):
        annotated_set = annotate_variants(
            variant_set=variant_set,
            assembly=config.assembly,
        )

    console.print(
        f"  [green]\u2713[/] Annotated [bold]{annotated_set.total_annotated:,}[/] variants"
    )
    console.print(
        f"  [green]\u2713[/] [bold]{annotated_set.total_with_protein_change:,}[/] "
        "change a protein (these are the interesting ones)"
    )
    _learn_more(explain_annotation())

    # ==================================================================
    # Step 5 — MHC binding prediction
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            5,
            TOTAL_STEPS,
            "Predict DLA Binding",
            f"Checking which mutant protein fragments fit into your pet's\n"
            f"DLA molecules (breed: [bold]{breed}[/], "
            f"alleles: {len(alleles)})...",
        )
    )
    _learn_more(explain_mhc_binding())

    with console.status("[bold cyan]Running MHC binding predictions...", spinner="dots"):
        from chloe_cli.commands import _get_predictor  # noqa: PLC0415

        predictor_backend = _get_predictor(config.predictor)
        prediction_results = predict_variants(
            annotated_set=annotated_set,
            predictor=predictor_backend,
            alleles=alleles,
            peptide_lengths=config.peptide_lengths,
        )

    console.print(
        f"  [green]\u2713[/] Evaluated [bold]{len(prediction_results.variant_predictions):,}[/] "
        "variants"
    )
    console.print(
        f"  [green]\u2713[/] Found [bold]{prediction_results.total_binders:,}[/] "
        "with at least one binding peptide"
    )

    # ==================================================================
    # Step 6 — Neoantigen ranking
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            6,
            TOTAL_STEPS,
            "Rank Neoantigens",
            "Scoring and ranking candidates by binding strength, mutation\n"
            "type, and other factors to find the best vaccine targets...",
        )
    )

    with console.status("[bold cyan]Ranking neoantigen candidates...", spinner="dots"):
        ranked = rank_neoantigens(
            prediction_results=prediction_results,
            top_n=config.top_n,
        )

    n_candidates = len(ranked.candidates)
    console.print(f"  [green]\u2713[/] Identified [bold]{n_candidates:,}[/] neoantigen candidates")
    _learn_more(explain_neoantigens(n_candidates))

    # Show top candidates in a table
    if ranked.candidates:
        _print_candidate_table(ranked.candidates[: config.top_n])

    # ==================================================================
    # Step 7 — Report generation
    # ==================================================================
    console.print()
    console.print(
        _step_panel(
            7,
            TOTAL_STEPS,
            "Generate Report",
            "Creating a detailed HTML report with all results...",
        )
    )

    output_path = Prompt.ask("  Output file path", default="report.html")

    with console.status("[bold cyan]Generating report...", spinner="dots"):
        generate_report(
            ranked_neoantigens=ranked,
            config=config,
            output_path=output_path,
        )

    console.print(f"  [green]\u2713[/] Report saved to [bold]{output_path}[/]")
    console.print()

    # ------------------------------------------------------------------
    # Closing
    # ------------------------------------------------------------------
    console.print(
        Panel(
            explain_next_steps(),
            title="[bold cyan]What happens next?[/]",
            border_style=HEADING_BORDER,
            padding=(1, 2),
        )
    )
    console.print()
    console.print("[dim]Thank you for using Chloe OS.[/]")
    console.print()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _ask_file_path(label: str) -> str:
    """Prompt for a file path, validating that the file exists."""
    while True:
        raw = Prompt.ask(f"  {label}")
        path = Path(raw).expanduser().resolve()
        if path.is_file():
            return str(path)
        console.print(f"  [red]File not found:[/] {path}")


def _print_candidate_table(candidates: list) -> None:  # noqa: ANN001
    """Display a Rich table with the top neoantigen candidates."""
    from chloe_core.models import NeoantigenCandidate

    typed: list[NeoantigenCandidate] = candidates
    console.print()

    table = Table(
        title="Top Neoantigen Candidates",
        show_lines=True,
        title_style="bold cyan",
    )
    table.add_column("#", style="bold", justify="right", width=4)
    table.add_column("Gene", style="cyan")
    table.add_column("Mutation", style="magenta")
    table.add_column("Best Peptide", style="green")
    table.add_column("DLA Allele", style="blue")
    table.add_column("IC50 (nM)", justify="right")
    table.add_column("Score", justify="right", style="bold yellow")

    for c in typed:
        ic50_str = f"{c.best_ic50:.1f}" if c.best_ic50 is not None else "-"
        binding_style = ""
        if c.best_ic50 is not None:
            if c.best_ic50 < 50:
                binding_style = "bold green"
            elif c.best_ic50 < 500:
                binding_style = "green"
            else:
                binding_style = "dim"

        table.add_row(
            str(c.rank),
            c.gene or "-",
            c.protein_change or "-",
            c.best_peptide or "-",
            c.best_allele or "-",
            f"[{binding_style}]{ic50_str}[/]" if binding_style else ic50_str,
            f"{c.composite_score:.3f}",
        )

    console.print(table)
    console.print()
