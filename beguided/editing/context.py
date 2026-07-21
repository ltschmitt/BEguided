"""Trinucleotide edit-context lookup for ABE base editors.

The bundled table maps a sense-strand ``XAY`` trinucleotide (middle base = A,
the ABE target) to a normalised editing proportion derived from Marquart
et al. 2021 (ABE8e data, scaled so the maximum = 1.0). For ``-``-strand ABE
edits the genome 3-mer is reverse-complemented before the lookup, since the
mRNA-sense context is what the paper measured.

CBE rows return ``None`` because the table is ABE-specific. Approximating
other ABE editors from the ABE8e data is supported but should be flagged as
such to downstream consumers (see :attr:`EditContextTable.source_label`).

Citation: Marquart et al. 2021, https://doi.org/10.1038/s41467-021-25375-z.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Dict, Optional, Union

import pandas as pd

from beguided.genome.sequence import reverse_complement

logger = logging.getLogger("beguided.editing.context")

# ABE editing direction on each strand (matches load_base_editor_profiles):
#   sense     A -> G      (+ strand guide, middle of NAN = A)
#   antisense T -> C      (- strand guide, middle of NAN (mRNA) = A)
_ABE_SENSE = ("A", "G")
_ABE_ANTISENSE = ("T", "C")

_SENSE_KEY_PATTERN = re.compile(r"^[ACGT]A[ACGT]$")
_EXPECTED_SENSE_KEYS = frozenset({
    "AAA", "AAC", "AAG", "AAT",
    "CAA", "CAC", "CAG", "CAT",
    "GAA", "GAC", "GAG", "GAT",
    "TAA", "TAC", "TAG", "TAT",
})

_TABLE_FILENAME = "edit_context_norm_proportion.tsv"


def _default_table_path() -> Path:
    """Return the path to the bundled edit-context table TSV."""
    return Path(
        resources.files("beguided.editing.data")
    ) / _TABLE_FILENAME


@dataclass
class EditContextTable:
    """In-memory trinucleotide -> normalised editing proportion table.

    Attributes:
        table: Sense-strand ``XAY`` keys (middle = A) -> proportions in [0, 1].
        antisense_table: Precomputed reverse-complement keys for fast lookup
            on ``-``-strand ABE guides.
        source_label: Human-readable provenance string. Surfaced in the
            logger so downstream operators see where the numbers came from.
        approximate_for_non_ABE8e_ABE: When True, the scores should be treated
            as an ABE8e-derived approximation for every ABE profile (the
            current default; the paper itself only measured ABE8e).
    """

    table: Dict[str, float]
    antisense_table: Dict[str, float] = field(init=False)
    source_label: str = "Marquart et al. 2021 (ABE8e, normalised)"
    approximate_for_non_ABE8e_ABE: bool = True

    def __post_init__(self) -> None:
        self.antisense_table = {
            reverse_complement(k): v for k, v in self.table.items()
        }


def load_edit_context_table(path: Union[str, Path, None] = None) -> EditContextTable:
    """Load the trinucleotide context TSV from *path* (or the bundled copy).

    The file must have ``trinuc`` and ``proportion`` columns, contain exactly
    the 16 expected sense keys (``AAA``, ``AAC``, ...), and every key must
    match the pattern ``^[ACGT]A[ACGT]$`` (middle base = A).

    Args:
        path: Optional path to a custom TSV. Defaults to the bundled copy
            (``beguided/editing/data/edit_context_norm_proportion.tsv``).

    Returns:
        Loaded :class:`EditContextTable`.
    """
    if path is None:
        path = _default_table_path()
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Edit-context table not found: {path}")

    df = pd.read_csv(path, sep="\t")
    if "trinuc" not in df.columns or "proportion" not in df.columns:
        raise ValueError(
            f"Edit-context table must have 'trinuc' and 'proportion' columns, "
            f"got: {', '.join(df.columns)}"
        )

    keys = df["trinuc"].astype(str).str.upper().tolist()
    bad = [k for k in keys if not _SENSE_KEY_PATTERN.match(k)]
    if bad:
        raise ValueError(
            f"Edit-context table contains invalid sense keys: {', '.join(bad)}. "
            f"Expected 16 sense keys of the form ^[ACGT]A[ACGT]$ (middle base = A)."
        )

    seen = set(keys)
    missing = _EXPECTED_SENSE_KEYS - seen
    if missing:
        raise ValueError(
            f"Edit-context table missing sense keys: {', '.join(sorted(missing))}"
        )

    table = {
        k: float(v)
        for k, v in zip(df["trinuc"].astype(str).str.upper(), df["proportion"])
    }
    logger.info(
        "Loaded %d trinucleotides from edit-context table: %s",
        len(table), path,
    )
    return EditContextTable(table=table)


def is_abe(base_from: str, base_to: str) -> bool:
    """Return True if the (from, to) editing direction is ABE (sense or antisense).

    Matches the per-strand ``BaseEditorConfig.target_base``/``desired_base``
    values that ``load_base_editor_profiles`` writes: ``A>G`` for ``+`` strand
    guides and ``T>C`` for ``-`` strand guides from the same ABE profile.
    """
    fb = str(base_from).upper()
    tb = str(base_to).upper()
    return (fb, tb) in (_ABE_SENSE, _ABE_ANTISENSE)


def lookup(ect: EditContextTable, trinuc: str) -> Optional[float]:
    """Return the normalised proportion for a sense-strand *trinuc*."""
    if not trinuc:
        return None
    return ect.table.get(trinuc.upper())


def lookup_from_genome(
    ect: EditContextTable,
    genome_seq: str,
    position: int,
    strand: str,
    target_base: str,
) -> Optional[float]:
    """Return the normalised proportion for the edit at *position* on *strand*.

    The 3-mer is built from the two genomic neighbours of *position*
    (``genome_seq[position-2]`` and ``genome_seq[position]``) with
    *target_base* as the middle base. Using *target_base* (rather than the
    actual reference base) ensures the lookup still works when the user's
    target file already encodes the desired edit (e.g. a ``G`` at an ``A>G``
    target) — the Marquart context is defined by the ABE target adenine,
    not by whatever the reference happens to carry.

    For the ``-`` strand the 3-mer is reverse-complemented before lookup so
    the context is in mRNA sense, matching how Marquart et al. 2021 reported
    their data.

    Returns ``None`` when:
    - one of the neighbour bases is outside the chromosome (chromosome end),
    - any of the three bases contains a non-ACGT character (e.g. ``N``),
    - the resulting trinucleotide is missing from the table (defensive;
      should not happen given :func:`load_edit_context_table` validation).

    The caller is responsible for verifying the edit is ABE (sense or
    antisense) before calling this helper; see :func:`is_abe`.
    """
    if position - 2 < 0 or position + 1 > len(genome_seq):
        return None
    left = genome_seq[position - 2].upper()
    right = genome_seq[position].upper()
    middle = str(target_base).upper()
    if left not in "ACGT" or right not in "ACGT" or middle not in "ACGT":
        return None
    raw = left + middle + right
    key = reverse_complement(raw) if strand == "-" else raw
    return ect.table.get(key)


__all__ = [
    "EditContextTable",
    "load_edit_context_table",
    "is_abe",
    "lookup",
    "lookup_from_genome",
]
