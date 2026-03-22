"""Chloe OS CLI — Typer application entrypoint.

Provides the ``chloe`` command with two primary sub-commands:

- ``chloe guided``  — interactive step-by-step walkthrough
- ``chloe run``     — direct pipeline execution with CLI flags
"""

from __future__ import annotations

import typer
from rich.console import Console
from rich.panel import Panel

from chloe_cli import __version__

app = typer.Typer(
    name="chloe",
    help="Chloe OS \u2014 Personalized cancer vaccine pipeline for pets",
    no_args_is_help=True,
    rich_markup_mode="rich",
)

console = Console()


# ------------------------------------------------------------------
# Sub-commands
# ------------------------------------------------------------------


@app.command()
def guided() -> None:
    """Interactive step-by-step walkthrough of the neoantigen pipeline.

    Walks you through each stage with plain-language explanations,
    input prompts, and progress spinners.
    """
    from chloe_cli.guided import guided as _guided  # noqa: PLC0415

    _guided()


@app.command()
def run(
    vcf: str = typer.Option(..., "--vcf", help="Path to the tumor VCF file."),
    normal: str | None = typer.Option(None, "--normal", help="Path to the matched-normal VCF."),
    species: str = typer.Option("canine", "--species", help="Species (currently only 'canine')."),
    breed: str | None = typer.Option(None, "--breed", help="Breed for DLA allele selection."),
    predictor: str = typer.Option("mhcflurry", "--predictor", help="MHC predictor backend."),
    top_n: int = typer.Option(20, "--top-n", help="Number of top candidates to report."),
    output: str = typer.Option("report.html", "--output", "-o", help="Output report file path."),
) -> None:
    """Run the full neoantigen identification pipeline.

    Executes load \u2192 annotate \u2192 predict \u2192 rank \u2192 report in a single
    pass, printing progress to the terminal.
    """
    from pathlib import Path  # noqa: PLC0415

    from chloe_cli.commands import run as _run  # noqa: PLC0415

    vcf_path = Path(vcf)
    if not vcf_path.is_file():
        console.print(f"[bold red]Error:[/] VCF file not found: {vcf}")
        raise typer.Exit(code=1)

    normal_path: Path | None = None
    if normal is not None:
        normal_path = Path(normal)
        if not normal_path.is_file():
            console.print(f"[bold red]Error:[/] Normal VCF file not found: {normal}")
            raise typer.Exit(code=1)

    _run(
        vcf=vcf_path,
        normal=normal_path,
        species=species,
        breed=breed,
        predictor=predictor,
        top_n=top_n,
        output=Path(output),
    )


@app.command()
def version() -> None:
    """Show the Chloe OS version."""
    console.print(
        Panel.fit(
            f"[bold cyan]Chloe OS[/]  v{__version__}",
            border_style="cyan",
        )
    )


@app.command()
def setup(
    species: str = typer.Option("canine", "--species", help="Species to download data for."),
    assembly: str | None = typer.Option(None, "--assembly", help="Genome assembly name."),
) -> None:
    """Download reference genome data and annotation databases.

    This prepares the local environment so that ``chloe run`` can
    execute without internet access.

    [dim]Note: This is a placeholder \u2014 genome download is not yet implemented.[/]
    """
    console.print()
    console.print(
        Panel(
            f"[bold cyan]Setup[/]\n\n"
            f"Species:  [bold]{species}[/]\n"
            f"Assembly: [bold]{assembly or 'default'}[/]\n\n"
            "[yellow]Genome data download is not yet implemented.[/]\n"
            "This command will eventually download reference genomes,\n"
            "VEP cache files, and DLA allele databases.",
            border_style="cyan",
            padding=(1, 2),
        )
    )
    console.print()


# ------------------------------------------------------------------
# Entrypoint
# ------------------------------------------------------------------

if __name__ == "__main__":
    app()
