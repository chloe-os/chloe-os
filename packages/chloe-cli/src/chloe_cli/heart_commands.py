"""CLI commands for chloe-heart cardiac monitoring analysis.

Provides the ``chloe heart`` subcommand group with guided and analyze modes.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.prompt import Confirm, Prompt
from rich.table import Table

heart_app = typer.Typer(
    name="heart",
    help="Cardiac monitoring analysis for dogs — ECG, arrhythmia detection, health scoring.",
    no_args_is_help=True,
)

console = Console()


@heart_app.command()
def analyze(
    ecg: str | None = typer.Option(None, "--ecg", help="Path to ECG file (CSV, EDF, or WFDB)."),
    hr_data: str | None = typer.Option(None, "--hr-data", help="Path to HR/HRV data file (CSV)."),
    breed: str | None = typer.Option(None, "--breed", help="Dog breed for baseline selection."),
    size: str | None = typer.Option(
        None, "--size", help="Dog size if breed unknown (giant/large/medium/small/toy)."
    ),
    output: str = typer.Option("cardiac_report.html", "--output", "-o", help="Output report path."),
) -> None:
    """Analyze ECG or heart rate data for cardiac abnormalities."""
    if ecg is None and hr_data is None:
        console.print("[bold red]Error:[/] Provide --ecg or --hr-data.")
        raise typer.Exit(code=1)

    source = ecg or hr_data
    source_path = Path(source)  # type: ignore[arg-type]
    if not source_path.is_file():
        console.print(f"[bold red]Error:[/] File not found: {source}")
        raise typer.Exit(code=1)

    from chloe_heart.models import DogSize, HeartConfig, get_baseline  # noqa: PLC0415

    # Resolve breed baseline
    dog_size = None
    if size:
        try:
            dog_size = DogSize(size.lower())
        except ValueError as err:
            console.print(
                f"[bold red]Error:[/] Unknown size '{size}'. Use: giant/large/medium/small/toy."
            )
            raise typer.Exit(code=1) from err

    baseline = get_baseline(breed=breed, size=dog_size)
    config = HeartConfig(
        breed=breed,
        size=dog_size,
        ecg_source=ecg,
        hr_source=hr_data,
        output_path=output,
    )

    console.print()
    console.print(
        Panel(
            "[bold red]Chloe Heart[/] — Cardiac Monitoring Analysis\n\n"
            f"Source:  [bold]{source}[/]\n"
            f"Breed:  [bold]{breed or 'Not specified'}[/] "
            f"(Size: {baseline.size.value.title()})\n"
            f"Normal HR: [bold]{baseline.hr_min}–{baseline.hr_max} BPM[/]\n\n"
            "[dim]This is a screening tool. Not veterinary medical advice.[/]",
            border_style="red",
            padding=(1, 2),
        )
    )
    console.print()

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        # Stage 1: Ingest
        task = progress.add_task("Loading ECG data...", total=None)
        if ecg:
            from chloe_heart.ingest import load_ecg  # noqa: PLC0415

            ecg_signal = load_ecg(ecg)
            progress.update(
                task,
                description=f"[green]Loaded {ecg_signal.num_samples} samples "
                f"({ecg_signal.duration_seconds:.1f}s at {ecg_signal.sample_rate}Hz)[/]",
            )
        else:
            console.print("[yellow]HR/HRV data analysis not yet implemented in MVP.[/]")
            raise typer.Exit(code=0)

        # Stage 2: Preprocess
        task = progress.add_task("Preprocessing ECG signal...", total=None)
        from chloe_heart.preprocess import preprocess_ecg  # noqa: PLC0415

        processed = preprocess_ecg(ecg_signal)
        progress.update(
            task,
            description=f"[green]Detected {processed.num_beats} heartbeats "
            f"(mean HR: {processed.mean_hr_bpm:.0f} BPM)[/]",
        )

        # Stage 3: Analyze
        task = progress.add_task("Analyzing cardiac rhythm...", total=None)
        from chloe_heart.analysis import analyze_cardiac  # noqa: PLC0415

        analysis = analyze_cardiac(processed, baseline)
        progress.update(
            task,
            description=f"[green]Analysis complete — "
            f"{analysis.abnormal_beat_count} abnormal beats "
            f"({analysis.arrhythmia_burden_pct:.1f}% burden)[/]",
        )

        # Stage 4: Score
        task = progress.add_task("Computing cardiac health score...", total=None)
        from chloe_heart.scoring import score_cardiac_health  # noqa: PLC0415

        score = score_cardiac_health(analysis, baseline)
        risk_color = {"green": "green", "yellow": "yellow", "red": "red"}[score.risk_level.value]
        progress.update(
            task,
            description=f"[{risk_color}]Health Score: {score.overall_score:.0f}/100 "
            f"(Risk: {score.risk_level.value.upper()})[/{risk_color}]",
        )

        # Stage 6: Report
        task = progress.add_task("Generating report...", total=None)
        from chloe_heart.report import generate_cardiac_report  # noqa: PLC0415

        report_path = generate_cardiac_report(processed, analysis, score, config, output)
        progress.update(
            task,
            description=f"[green]Report saved → {report_path}[/]",
        )

    console.print()

    # Summary table
    table = Table(title="Cardiac Analysis Summary", border_style="red")
    table.add_column("Metric", style="bold")
    table.add_column("Value")
    table.add_column("Status")

    hr_status = (
        "[green]Normal[/]"
        if baseline.hr_min <= analysis.mean_hr_bpm <= baseline.hr_max
        else "[yellow]Outside range[/]"
    )
    table.add_row("Mean Heart Rate", f"{analysis.mean_hr_bpm:.0f} BPM", hr_status)
    table.add_row("Total Beats", str(analysis.total_beats), "")
    table.add_row(
        "Abnormal Beats",
        f"{analysis.abnormal_beat_count} ({analysis.arrhythmia_burden_pct:.1f}%)",
        "[green]Normal[/]" if analysis.arrhythmia_burden_pct < 1 else "[yellow]Elevated[/]",
    )
    table.add_row("SDNN", f"{analysis.hrv_metrics.sdnn:.1f} ms", "")
    table.add_row("RMSSD", f"{analysis.hrv_metrics.rmssd:.1f} ms", "")

    risk_style = {"green": "bold green", "yellow": "bold yellow", "red": "bold red"}[
        score.risk_level.value
    ]
    table.add_row(
        "Health Score",
        f"{score.overall_score:.0f}/100",
        f"[{risk_style}]{score.risk_level.value.upper()}[/{risk_style}]",
    )

    console.print(table)
    console.print()
    console.print(
        Panel(
            "[bold]Important:[/] This is a screening tool. Share this report with your "
            "veterinarian for professional interpretation.",
            border_style="yellow",
        )
    )
    console.print()


@heart_app.command()
def guided() -> None:
    """Interactive step-by-step cardiac analysis walkthrough."""
    console.print()
    console.print(
        Panel(
            "[bold red]Chloe Heart[/] — Guided Cardiac Analysis\n\n"
            "This will walk you through analyzing your dog's heart data\n"
            "step by step, explaining everything along the way.\n\n"
            "[dim italic]This is a screening tool. Not veterinary medical advice.\n"
            "Always consult a veterinary cardiologist for cardiac concerns.[/]",
            border_style="red",
            padding=(1, 2),
        )
    )
    console.print()

    # Step 1: Get ECG file
    ecg_path = Prompt.ask(
        "[bold]Path to your ECG file[/] (CSV format)\n"
        "  [dim]This is the recording from your pet's heart monitor[/]"
    )
    if not Path(ecg_path).is_file():
        console.print(f"[bold red]Error:[/] File not found: {ecg_path}")
        raise typer.Exit(code=1)

    # Step 2: Get breed
    breed = Prompt.ask(
        "\n[bold]What breed is your dog?[/]\n"
        "  [dim]This helps set the right heart rate baselines.\n"
        "  Examples: labrador, boxer, cavalier, chihuahua\n"
        "  Press Enter to skip[/]",
        default="",
    )
    breed = breed.strip() or None

    # Run analysis using the analyze command
    console.print()
    analyze(ecg=ecg_path, breed=breed, output="cardiac_report.html")

    if Confirm.ask("\n[bold]Would you like to open the report in your browser?[/]"):
        import subprocess  # noqa: PLC0415

        subprocess.run(["open", "cardiac_report.html"], check=False)  # noqa: S603, S607
