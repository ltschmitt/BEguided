"""Genome I/O operations."""

import logging
import re
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

logger = logging.getLogger("beguided.genome.io")


# Records whose header description contains any of these tokens are skipped.
# ``:scaffold`` covers both Ensembl's ``dna:scaffold`` (human primary
# assembly fixtures) and the EBI ``ensemblorganisms`` ``unmasked:scaffold``
# (mouse), dropping unplaced scaffolds so only chromosomes land in the dict.
SKIP_HEADER_TOKENS = (":scaffold",)


def load_genome(genome_file, *, skip_header_tokens: Iterable[str] = SKIP_HEADER_TOKENS):
    """Read genome file and return dictionary of chromosomes.

    FASTA records whose description contains any of the tokens in
    ``skip_header_tokens`` (default: ``":scaffold"``) are skipped. This
    drops assembly scaffolds and patches (e.g. ``KI270728.1``,
    ``GL000225.1``) so only the canonical chromosomes land in the dict.

    Args:
        genome_file: Path to FASTA file.
        skip_header_tokens: Substrings; if any appears in a record's full
            header line (after the ``>`` and id), that record is skipped.
            Pass an empty iterable to disable filtering.

    Returns:
        dict: chromosome_id -> sequence string.
    """
    skip_set = tuple(skip_header_tokens)
    genome_dict = {}
    skipped = 0
    with open(genome_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if skip_set and any(tok in record.description for tok in skip_set):
                skipped += 1
                continue
            chrom_id = str(record.id).split()[0]
            genome_dict[chrom_id] = str(record.seq)
    if skipped:
        logger.info(
            "Loaded %d chromosomes from %s (skipped %d scaffold/patch record(s))",
            len(genome_dict),
            genome_file,
            skipped,
        )
    return genome_dict


_ASSEMBLY_RE = re.compile(
    r"^(?:Homo_sapiens|Mus_musculus)\.([A-Za-z0-9._]+?)(?:\.dna\.|\.fa$|\.fasta$|$)"
)


def derive_genome_label(fasta_path: str) -> str:
    """Best-effort extraction of the assembly label from an Ensembl-style
    FASTA filename.

    Examples:
        ``Homo_sapiens.GRCh38.dna.primary_assembly.fa`` -> ``"GRCh38"``
        ``Homo_sapiens.GRCh38.112.dna.primary_assembly.fa`` -> ``"GRCh38.112"``
        ``Mus_musculus.GRCm39.2026_04.fa`` -> ``"GRCm39.2026_04"``
        ``/data/genomes/anything.fa`` -> ``"GRCh38"`` (fallback)

    The fallback is ``"GRCh38"`` because the on-disk fixtures shipped with
    BEguided use that assembly; it provides a sane default for any FASTA
    that doesn't match the Ensembl naming convention.
    """
    name = Path(fasta_path).name
    m = _ASSEMBLY_RE.match(name)
    if m:
        return m.group(1)
    return "GRCh38"
