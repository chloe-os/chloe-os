"""Direct CLI commands for running the Chloe OS pipeline.

Provides the ``chloe run`` command that wires up the full pipeline:
load_variants -> annotate_variants -> predict_variants ->
rank_neoantigens -> generate_report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.table import Table

from chloe_core.models import BREED_DLA_ALLELES, PipelineConfig, Species

console = Console()


def run(
    vcf: Annotated[
        Path,
        typer.Option("--vcf", help="Path to the tumor VCF file.", exists=True, readable=True),
    ],
    normal: Annotated[
        Optional[Path],
        typer.Option("--normal", help="Path to the matched-normal VCF file.", exists=True),
    ] = None,
    species: Annotated[
        str,
        typer.Option("--species", help="Species (currently only 'canine')."),
    ] = "canine",
    breed: Annotated[
        Optional[str],
        typer.Option("--breed", help="Breed name for DLA allele selection (e.g. staffy, labrador)."),
    ] = None,
    predictor: Annotated[
        str,
        typer.Option("--predictor", help="MHC binding predictor backend."),
    ] = "mhcflurry",
    top_n: Annotated[
        int,
        typer.Option("--top-n", help="Number of top neoantigen candidates to include."),
    ] = 20,
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Output report file path."),
    ] = Path("report.html"),
) -> None:
    """Run the full neoantigen identification pipeline."""
    # Lazy imports so the CLI stays responsive even when heavy deps are loading.
    from chloe_core.annotation import annotate_variants  # type: ignore[attr-defined]
    from chloe_core.prediction.base import predict_variants
    from chloe_core.ranking import rank_neoantigens  # type: ignore[attr-defined]
    from chloe_core.report import generate_report  # type: ignore[attr-defined]
    from chloe_core.variants import load_variants  # type: ignore[attr-defined]

    # ------------------------------------------------------------------
    # Build pipeline configuration
    # ------------------------------------------------------------------
    try:
        species_enum = Species(species)
    except ValueError:
        console.print(f"[bold red]Error:[/] Unknown species '{species}'.")
        raise typer.Exit(code=1)

    alleles = BREED_DLA_ALLELES.get(breed or "default", BREED_DLA_ALLELES["default"])

    config = PipelineConfig(
        species=species_enum,
        breed=breed,
        predictor=predictor,
        alleles=alleles,
        top_n=top_n,
    )

    # ------------------------------------------------------------------
    # Banner
    # ------------------------------------------------------------------
    console.print()
    console.print(
        Panel.fit(
            "[bold cyan]Chloe OS[/] — Personalized cancer vaccine pipeline",
            subtitle=f"species={species}  breed={breed or 'default'}  predictor={predictor}",
        )
    )
    console.print()

    # ------------------------------------------------------------------
    # Pipeline execution with progress
    # ------------------------------------------------------------------
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
        transient=False,
    ) as progress:
        overall = progress.add_task("[bold]Pipeline", total=5)

        # Stage 1 — Load variants
        progress.update(overall, description="[bold]Loading variants...")
        variant_set = load_variants(
            tumor_vcf=str(vcf),
            normal_vcf=str(normal) if normal else None,
            species=species_enum,
            min_quality=config.min_quality,
            min_read_depth=config.min_read_depth,
            min_allele_frequency=config.min_allele_frequency,
        )
        progress.advance(overall)
        console.print(
            f"  [green]\u2713[/] Loaded [bold]{variant_set.somatic_count:,}[/] somatic variants "
            f"(from {variant_set.total_variants_raw:,} raw)"
        )

        # Stage 2 — Annotate variants
        progress.update(overall, description="[bold]Annotating variants...")
        annotated_set = annotate_variants(
            variant_set=variant_set,
            assembly=config.assembly,
        )
        progress.advance(overall)
        console.print(
            f"  [green]\u2713[/] Annotated [bold]{annotated_set.total_annotated:,}[/] variants "
            f"({annotated_set.total_with_protein_change:,} with protein changes)"
        )

        # Stage 3 — MHC binding prediction
        progress.update(overall, description="[bold]Predicting MHC binding...")
        prediction_results = predict_variants(
            annotated_set=annotated_set,
            predictor=_get_predictor(predictor),
            alleles=alleles,
            peptide_lengths=config.peptide_lengths,
        )
        progress.advance(overall)
        console.print(
            f"  [green]\u2713[/] Predicted binding for "
            f"[bold]{len(prediction_results.variant_predictions):,}[/] variants "
            f"({prediction_results.total_binders:,} binders)"
        )

        # Stage 4 — Rank neoantigens
        progress.update(overall, description="[bold]Ranking neoantigens...")
        ranked = rank_neoantigens(
            prediction_results=prediction_results,
            top_n=top_n,
        )
        progress.advance(overall)
        console.print(
            f"  [green]\u2713[/] Ranked [bold]{len(ranked.candidates):,}[/] neoantigen candidates"
        )

        # Stage 5 — Generate report
        progress.update(overall, description="[bold]Generating report...")
        generate_report(
            ranked_neoantigens=ranked,
            config=config,
            output_path=str(output),
        )
        progress.advance(overall)
        console.print(f"  [green]\u2713[/] Report saved to [bold]{output}[/]")

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    console.print()
    _print_summary_table(ranked, top_n)
    console.print()
    console.print(
        Panel(
            "[dim]This analysis is a research tool and is not a substitute for "
            "professional veterinary advice.[/]",
            title="[yellow]Disclaimer[/]",
            border_style="yellow",
        )
    )


def _get_predictor(name: str):  # noqa: ANN202
    """Instantiate the requested MHC predictor backend.

    Returns a concrete ``MHCPredictor`` instance.  Currently supports
    ``mhcflurry`` and ``netmhcpan``.  Falls back to a stub predictor
    until real backends are wired up.
    """
    # Future: import concrete predictors here.
    # from chloe_core.prediction.mhcflurry import MHCflurryPredictor
    # from chloe_core.prediction.netmhcpan import NetMHCpanPredictor
    raise typer.Exit(
        code=1,
    )


def _print_summary_table(ranked, top_n: int) -> None:  # noqa: ANN001
    """Print a Rich table of the top neoantigen candidates."""
    from chloe_core.models import RankedNeoantigens

    ranked_neo: RankedNeoantigens = ranked
    candidates = ranked_neo.candidates[:top_n]

    if not candidates:
        console.print("[yellow]No neoantigen candidates met the scoring threshold.[/]")
        return

    table = Table(
        title=f"Top {min(top_n, len(candidates))} Neoantigen Candidates",
        show_lines=True,
    )
    table.add_column("Rank", style="bold", justify="right", width=5)
    table.add_column("Gene", style="cyan")
    table.add_column("Mutation", style="magenta")
    table.add_column("Peptide", style="green")
    table.add_column("Allele", style="blue")
    table.add_column("IC50 (nM)", justify="right")
    table.add_column("Score", justify="right", style="bold yellow")

    for c in candidates:
        ic50_str = f"{c.best_ic50:.1f}" if c.best_ic50 is not None else "-"
        table.add_row(
            str(c.rank),
            c.gene or "-",
            c.protein_change or "-",
            c.best_peptide or "-",
            c.best_allele or "-",
            ic50_str,
            f"{c.composite_score:.3f}",
        )

    console.print(table)
