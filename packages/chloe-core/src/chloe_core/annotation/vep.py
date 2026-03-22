"""Variant annotation using Ensembl VEP (Variant Effect Predictor).

Stage 2 of the Chloe OS personalized cancer vaccine pipeline. Takes a
``VariantSet`` produced by Stage 1 and annotates each variant with its
gene/protein impact, consequence class, and flanking peptide sequences
suitable for downstream MHC binding prediction.

VEP can be executed either via a local installation or through the official
Ensembl VEP Docker container (``ensemblorg/ensembl-vep``). Docker mode is
the default and recommended for reproducibility.

Supported canine assemblies:
    * CanFam_GSD  (GSD 1.0 — UU_Cfam_GSD_1.0)
    * CanFam3.1   (legacy Broad assembly)
"""

from __future__ import annotations

import csv
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from chloe_core.models import (
    AnnotatedVariant,
    AnnotatedVariantSet,
    Variant,
    VariantConsequence,
    VariantSet,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Default number of flanking amino acids on each side of a mutation when
#: extracting peptide context from the VEP ``WildtypeProtein`` plugin output.
PEPTIDE_FLANK_LENGTH: int = 15

#: VEP Docker image used when ``use_docker=True``.
VEP_DOCKER_IMAGE: str = "ensemblorg/ensembl-vep:release_112"

#: Mapping from VEP consequence strings to our internal enum.  VEP may
#: report multiple consequences separated by commas — we match in priority
#: order and fall back to ``OTHER``.
_CONSEQUENCE_MAP: dict[str, VariantConsequence] = {
    "missense_variant": VariantConsequence.MISSENSE,
    "frameshift_variant": VariantConsequence.FRAMESHIFT,
    "stop_gained": VariantConsequence.NONSENSE,
    "inframe_insertion": VariantConsequence.INFRAME_INSERTION,
    "inframe_deletion": VariantConsequence.INFRAME_DELETION,
    "splice_donor_variant": VariantConsequence.SPLICE_SITE,
    "splice_acceptor_variant": VariantConsequence.SPLICE_SITE,
    "synonymous_variant": VariantConsequence.SYNONYMOUS,
}

#: VEP consequence terms ordered by severity (most severe first).  Used to
#: select the "worst" consequence when VEP reports multiple terms.
_CONSEQUENCE_PRIORITY: list[str] = [
    "frameshift_variant",
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "synonymous_variant",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


class VEPNotInstalledError(RuntimeError):
    """Raised when neither Docker nor a local VEP install is found."""


def _check_vep_available(use_docker: bool) -> None:
    """Verify that VEP is available, raising a descriptive error if not.

    Parameters
    ----------
    use_docker:
        If ``True``, check that ``docker`` is on the ``PATH`` and the daemon
        is responsive.  If ``False``, check for the ``vep`` command.

    Raises
    ------
    VEPNotInstalledError
        With installation instructions tailored to the mode.
    """
    if use_docker:
        if shutil.which("docker") is None:
            raise VEPNotInstalledError(
                "Docker is required to run VEP but was not found on PATH.\n"
                "Install Docker from https://docs.docker.com/get-docker/ and "
                "then pull the VEP image:\n"
                f"  docker pull {VEP_DOCKER_IMAGE}"
            )
        # Verify the Docker daemon is running
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise VEPNotInstalledError(
                "Docker is installed but the daemon is not running.\n"
                "Start Docker Desktop or the Docker service and retry."
            )
    else:
        if shutil.which("vep") is None:
            raise VEPNotInstalledError(
                "Ensembl VEP is not installed or not on PATH.\n"
                "Install it with:\n"
                "  # Via conda (recommended):\n"
                "  conda install -c bioconda ensembl-vep\n\n"
                "  # Or use Docker mode (recommended):\n"
                f"  docker pull {VEP_DOCKER_IMAGE}\n"
                "  annotate_variants(variant_set, use_docker=True)"
            )


# ---------------------------------------------------------------------------
# VEP input / output
# ---------------------------------------------------------------------------


def _write_vep_input(variants: list[Variant], output_path: str) -> None:
    """Write variants in VEP default input format.

    Each line is tab-delimited with the fields::

        chromosome  start  end  allele  strand  identifier

    For SNVs ``start == end``.  For insertions the start is the position
    *before* the insertion and end is ``start``.  For deletions end is
    ``start + len(ref) - 1``.

    Parameters
    ----------
    variants:
        Somatic variants to convert.
    output_path:
        Filesystem path where the VEP input file will be written.
    """
    with open(output_path, "w") as fh:
        for v in variants:
            ref_len = len(v.ref)
            alt_len = len(v.alt)

            if ref_len == 1 and alt_len == 1:
                # SNV
                start = v.pos
                end = v.pos
                allele = f"{v.ref}/{v.alt}"
            elif ref_len > alt_len:
                # Deletion — VEP expects the *deleted* bases only.  If the
                # VCF anchor base is shared (REF[0] == ALT[0]), skip it.
                if v.ref[0] == v.alt[0]:
                    start = v.pos + 1
                    end = v.pos + ref_len - 1
                    allele = f"{v.ref[1:]}/{v.alt[1:] if alt_len > 1 else '-'}"
                else:
                    start = v.pos
                    end = v.pos + ref_len - 1
                    allele = f"{v.ref}/{v.alt}"
            else:
                # Insertion — VEP uses the positions *flanking* the insertion.
                if v.ref[0] == v.alt[0]:
                    start = v.pos
                    end = v.pos + ref_len - 1
                    allele = f"{v.ref[1:] if ref_len > 1 else '-'}/{v.alt[1:]}"
                else:
                    start = v.pos
                    end = v.pos
                    allele = f"{v.ref}/{v.alt}"

            line = "\t".join(
                [
                    v.chrom,
                    str(start),
                    str(end),
                    allele,
                    "+",
                    v.variant_id,
                ]
            )
            fh.write(line + "\n")

    logger.info("Wrote %d variants to VEP input file: %s", len(variants), output_path)


def _build_vep_command(
    input_path: str,
    output_path: str,
    assembly: str,
    use_docker: bool,
    vep_cache_dir: str | None,
) -> list[str]:
    """Construct the VEP command line.

    Parameters
    ----------
    input_path:
        Path to the VEP input file.
    output_path:
        Path where VEP should write its output.
    assembly:
        Genome assembly name (``CanFam_GSD`` or ``CanFam3.1``).
    use_docker:
        Whether to wrap the command in ``docker run``.
    vep_cache_dir:
        Path to the VEP annotation cache.  If ``None`` a default is chosen
        (``~/.vep`` for local, bind-mounted into Docker).

    Returns
    -------
    list[str]
        Command tokens suitable for :func:`subprocess.run`.
    """
    cache_dir = vep_cache_dir or str(Path.home() / ".vep")

    # Core VEP flags — tab output with the fields we need.
    vep_flags: list[str] = [
        "--input_file",
        input_path,
        "--output_file",
        output_path,
        "--species",
        "canis_lupus_familiaris",
        "--assembly",
        assembly,
        "--cache",
        "--dir_cache",
        cache_dir,
        "--offline",
        "--tab",
        "--no_stats",
        "--force_overwrite",
        "--symbol",
        "--protein",
        "--biotype",
        "--canonical",
        "--hgvs",
        "--fields",
        (
            "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,"
            "Consequence,cDNA_position,CDS_position,Protein_position,"
            "Amino_acids,Codons,SYMBOL,BIOTYPE,IMPACT,HGVSp,CANONICAL"
        ),
    ]

    if use_docker:
        # All paths inside the container live under /data.
        container_input = f"/data/{Path(input_path).name}"
        container_output = f"/data/{Path(output_path).name}"
        container_cache = "/opt/vep/.vep"

        # Replace host paths with container paths in vep_flags.
        docker_vep_flags = []
        for flag in vep_flags:
            if flag == input_path:
                docker_vep_flags.append(container_input)
            elif flag == output_path:
                docker_vep_flags.append(container_output)
            elif flag == cache_dir:
                docker_vep_flags.append(container_cache)
            else:
                docker_vep_flags.append(flag)

        data_dir = str(Path(input_path).parent)
        cmd = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{data_dir}:/data",
            "-v",
            f"{cache_dir}:/opt/vep/.vep",
            VEP_DOCKER_IMAGE,
            "vep",
            *docker_vep_flags,
        ]
    else:
        cmd = ["vep", *vep_flags]

    return cmd


def _run_vep(
    input_path: str,
    output_path: str,
    assembly: str = "CanFam_GSD",
    use_docker: bool = True,
    vep_cache_dir: str | None = None,
) -> str:
    """Execute Ensembl VEP and return the path to the output file.

    Parameters
    ----------
    input_path:
        Path to VEP-formatted input file.
    output_path:
        Desired path for VEP tab-delimited output.
    assembly:
        Reference genome assembly name.
    use_docker:
        Run VEP inside the Docker container if ``True``.
    vep_cache_dir:
        Optional override for the VEP annotation cache directory.

    Returns
    -------
    str
        Path to the VEP output file (same as *output_path*).

    Raises
    ------
    VEPNotInstalledError
        If VEP / Docker is not available.
    subprocess.CalledProcessError
        If VEP exits with a non-zero status.
    """
    _check_vep_available(use_docker)

    cmd = _build_vep_command(input_path, output_path, assembly, use_docker, vep_cache_dir)

    logger.info("Running VEP: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=3600,  # 1 hour timeout
    )

    if result.returncode != 0:
        logger.error("VEP stderr:\n%s", result.stderr)
        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            output=result.stdout,
            stderr=result.stderr,
        )

    logger.info("VEP completed successfully. Output: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Consequence mapping
# ---------------------------------------------------------------------------


def _map_consequence(vep_consequence: str) -> VariantConsequence:
    """Map a VEP consequence string to a :class:`VariantConsequence` enum value.

    VEP may report multiple consequence terms separated by commas (e.g.
    ``"missense_variant,splice_region_variant"``).  This function picks the
    most severe term according to :data:`_CONSEQUENCE_PRIORITY` and maps it
    via :data:`_CONSEQUENCE_MAP`, falling back to
    :attr:`VariantConsequence.OTHER`.

    Parameters
    ----------
    vep_consequence:
        Raw consequence string from VEP output (may contain multiple
        comma-separated terms).

    Returns
    -------
    VariantConsequence
        The mapped enum member.
    """
    terms = [t.strip() for t in vep_consequence.split(",")]

    # Walk through priority list and return the first match.
    for priority_term in _CONSEQUENCE_PRIORITY:
        if priority_term in terms:
            return _CONSEQUENCE_MAP[priority_term]

    # Check the map for any term not in the priority list.
    for term in terms:
        if term in _CONSEQUENCE_MAP:
            return _CONSEQUENCE_MAP[term]

    return VariantConsequence.OTHER


# ---------------------------------------------------------------------------
# Peptide extraction
# ---------------------------------------------------------------------------


def _extract_peptide_context(
    protein_position: str,
    amino_acids: str,
    consequence: VariantConsequence,
    flank: int = PEPTIDE_FLANK_LENGTH,
    full_protein: str | None = None,
) -> tuple[str | None, str | None]:
    """Derive wildtype and mutant peptide windows around a mutation.

    When a full protein sequence is available (e.g. from the VEP
    ``WildtypeProtein`` plugin), the flanking context is sliced directly.
    Otherwise a minimal representation is constructed from the ``Amino_acids``
    field (e.g. ``"V/E"``).

    Parameters
    ----------
    protein_position:
        The 1-based protein position string from VEP (e.g. ``"600"``).
    amino_acids:
        The amino-acid change string from VEP (e.g. ``"V/E"``).
    consequence:
        Mapped consequence enum — only protein-altering variants are handled.
    flank:
        Number of flanking amino acids on each side of the mutation.
    full_protein:
        Optional complete wildtype protein sequence.

    Returns
    -------
    tuple[str | None, str | None]
        ``(wildtype_peptide, mutant_peptide)`` or ``(None, None)`` if the
        variant does not produce a meaningful peptide change.
    """
    if consequence in (VariantConsequence.SYNONYMOUS, VariantConsequence.OTHER):
        return None, None

    if not protein_position or protein_position == "-":
        return None, None

    if not amino_acids or amino_acids == "-":
        return None, None

    # Parse the position — may be a range for indels (e.g. "600-602").
    pos_match = re.match(r"(\d+)", protein_position)
    if not pos_match:
        return None, None
    pos = int(pos_match.group(1)) - 1  # 0-based

    # Parse amino-acid change.
    parts = amino_acids.split("/")
    if len(parts) == 2:
        wt_aa, mut_aa = parts
    elif len(parts) == 1:
        # Frameshifts sometimes only show the WT residue.
        wt_aa = parts[0]
        mut_aa = "X"  # unknown mutant
    else:
        return None, None

    # Handle the "-" placeholder for insertions/deletions.
    wt_aa = wt_aa if wt_aa != "-" else ""
    mut_aa = mut_aa if mut_aa != "-" else ""

    if full_protein is not None and len(full_protein) > pos:
        # Slice from the full protein sequence.
        start = max(0, pos - flank)
        end_wt = min(len(full_protein), pos + len(wt_aa) + flank)

        wt_peptide = full_protein[start:end_wt]

        # Build the mutant peptide by substituting in the mutant residue(s).
        mut_protein = full_protein[:pos] + mut_aa + full_protein[pos + len(wt_aa):]
        end_mut = min(len(mut_protein), pos + len(mut_aa) + flank)
        mut_peptide = mut_protein[start:end_mut]

        return wt_peptide, mut_peptide

    # Fallback: construct a minimal context using only the amino-acid change.
    # Pad with 'X' to represent unknown flanking residues.
    wt_peptide = "X" * flank + wt_aa + "X" * flank
    mut_peptide = "X" * flank + mut_aa + "X" * flank

    return wt_peptide, mut_peptide


# ---------------------------------------------------------------------------
# VEP output parsing
# ---------------------------------------------------------------------------


def _parse_vep_output(
    output_path: str,
    variants: list[Variant],
) -> list[AnnotatedVariant]:
    """Parse VEP tab-delimited output into :class:`AnnotatedVariant` objects.

    VEP may emit multiple annotation rows per input variant (one per
    transcript).  This function keeps only the *canonical* transcript row
    for each variant.  If no canonical annotation is present the first row
    is used.

    Parameters
    ----------
    output_path:
        Path to the VEP tab-delimited output file.
    variants:
        Original variant list — used to look up :class:`Variant` objects by
        their ``variant_id``.

    Returns
    -------
    list[AnnotatedVariant]
        One annotated variant per input variant that VEP produced output for.
    """
    variant_lookup: dict[str, Variant] = {v.variant_id: v for v in variants}

    # Collect all rows keyed by uploaded variant id.
    rows_by_id: dict[str, list[dict[str, str]]] = {}

    with open(output_path) as fh:
        # Skip VEP header comment lines.
        header_line: str | None = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                # This is the column header line.
                header_line = line.lstrip("#").strip()
                break

        if header_line is None:
            logger.warning("VEP output file has no header: %s", output_path)
            return []

        reader = csv.DictReader(
            fh,
            fieldnames=header_line.split("\t"),
            delimiter="\t",
        )

        for row in reader:
            var_id = row.get("Uploaded_variation", "").strip()
            if var_id:
                rows_by_id.setdefault(var_id, []).append(row)

    # Build annotated variants — prefer canonical transcript.
    annotated: list[AnnotatedVariant] = []

    for var_id, rows in rows_by_id.items():
        variant = variant_lookup.get(var_id)
        if variant is None:
            logger.debug("VEP output references unknown variant: %s", var_id)
            continue

        # Pick the canonical row if available.
        chosen = rows[0]
        for row in rows:
            if row.get("CANONICAL", "").upper() == "YES":
                chosen = row
                break

        consequence = _map_consequence(chosen.get("Consequence", ""))
        amino_acids = chosen.get("Amino_acids", "")
        protein_position = chosen.get("Protein_position", "")
        codons = chosen.get("Codons", "")

        # Build protein change string (e.g. "V600E").
        protein_change: str | None = None
        if amino_acids and "/" in amino_acids and protein_position and protein_position != "-":
            aa_parts = amino_acids.split("/")
            if len(aa_parts) == 2:
                pos_match = re.match(r"(\d+)", protein_position)
                if pos_match:
                    protein_change = f"{aa_parts[0]}{pos_match.group(1)}{aa_parts[1]}"

        # Extract peptide context.
        wt_peptide, mut_peptide = _extract_peptide_context(
            protein_position=protein_position,
            amino_acids=amino_acids,
            consequence=consequence,
        )

        # Parse HGVSp for a more precise protein change notation if available.
        hgvsp = chosen.get("HGVSp", "")
        if hgvsp and ":" in hgvsp:
            # e.g. "ENSCAFP00000012345:p.Val600Glu"
            protein_change = protein_change or hgvsp.split(":")[-1]

        av = AnnotatedVariant(
            variant=variant,
            gene_symbol=chosen.get("SYMBOL") or None,
            gene_id=chosen.get("Gene") or None,
            transcript_id=chosen.get("Feature") or None,
            consequence=consequence,
            protein_change=protein_change,
            codon_change=codons if codons and codons != "-" else None,
            wildtype_peptide=wt_peptide,
            mutant_peptide=mut_peptide,
            impact=chosen.get("IMPACT") or None,
        )
        annotated.append(av)

    # Include variants that VEP did not annotate (intergenic, etc.).
    annotated_ids = {av.variant.variant_id for av in annotated}
    for v in variants:
        if v.variant_id not in annotated_ids:
            annotated.append(
                AnnotatedVariant(
                    variant=v,
                    consequence=VariantConsequence.OTHER,
                    impact="MODIFIER",
                )
            )

    logger.info(
        "Parsed %d annotated variants (%d with protein change) from VEP output",
        len(annotated),
        sum(1 for a in annotated if a.has_protein_change),
    )
    return annotated


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def annotate_variants(
    variant_set: VariantSet,
    *,
    use_docker: bool = True,
    vep_cache_dir: str | None = None,
) -> AnnotatedVariantSet:
    """Annotate somatic variants using Ensembl VEP (Stage 2).

    This is the main entry point for variant annotation.  It:

    1. Writes the variants in VEP input format to a temporary file.
    2. Invokes VEP (via Docker or local install) with canine-specific settings.
    3. Parses the tab-delimited VEP output.
    4. Returns an :class:`AnnotatedVariantSet` with gene, protein, and peptide
       information for each variant.

    Parameters
    ----------
    variant_set:
        Stage 1 output containing filtered somatic variants.
    use_docker:
        If ``True`` (default), run VEP inside the official Docker container.
        If ``False``, expect ``vep`` to be on ``PATH``.
    vep_cache_dir:
        Path to the VEP annotation cache directory.  Defaults to
        ``~/.vep``.  The cache must contain data for the canine species
        and the target assembly.

    Returns
    -------
    AnnotatedVariantSet
        Annotated variants ready for Stage 3 (MHC binding prediction).

    Raises
    ------
    VEPNotInstalledError
        If VEP or Docker is not available.
    subprocess.CalledProcessError
        If VEP fails during execution.
    ValueError
        If the variant set is empty.

    Example
    -------
    >>> from chloe_core.annotation import annotate_variants
    >>> annotated = annotate_variants(variant_set, use_docker=True)
    >>> for av in annotated.protein_changing_variants:
    ...     print(av.gene_symbol, av.protein_change)
    """
    if not variant_set.variants:
        raise ValueError("VariantSet contains no variants to annotate.")

    assembly = variant_set.assembly
    logger.info(
        "Annotating %d variants with VEP (assembly=%s, docker=%s)",
        len(variant_set.variants),
        assembly,
        use_docker,
    )

    with tempfile.TemporaryDirectory(prefix="chloe_vep_") as tmpdir:
        input_path = str(Path(tmpdir) / "vep_input.txt")
        output_path = str(Path(tmpdir) / "vep_output.txt")

        # Step 1: Write VEP input.
        _write_vep_input(variant_set.variants, input_path)

        # Step 2: Run VEP.
        _run_vep(
            input_path=input_path,
            output_path=output_path,
            assembly=assembly,
            use_docker=use_docker,
            vep_cache_dir=vep_cache_dir,
        )

        # Step 3: Parse output.
        annotated = _parse_vep_output(output_path, variant_set.variants)

    # Detect VEP version from Docker image tag or command output.
    vep_version = VEP_DOCKER_IMAGE.split(":")[-1] if use_docker else None

    result = AnnotatedVariantSet(
        annotated_variants=annotated,
        total_annotated=len(annotated),
        total_with_protein_change=sum(1 for a in annotated if a.has_protein_change),
        vep_version=vep_version,
        assembly=assembly,
    )

    logger.info(
        "Annotation complete: %d total, %d with protein changes",
        result.total_annotated,
        result.total_with_protein_change,
    )
    return result
