"""Parser for the single-target input mode.

The single input mode accepts one manually-entered target: a chromosome
name, a 1-based genomic position, and an editing direction such as
``A>G``. It produces the same 4-column targets ``DataFrame`` as the
spreadsheet and sequence paths so the existing pipeline
(``run_guidefinder`` → ``run_pipeline``) can consume it without changes.

Unlike sequence mode, single mode supplies no custom genome: the target
points at a position on the operator-configured reference genome, so the
job behaves exactly like a spreadsheet-mode job downstream (bystander
annotation and off-target screening stay available).

Example:
    ``chromosome="chr1", position="100", edit="A>G"`` yields the row
    ``{"chromosome": "chr1", "position": 100, "target_base": "A",
    "desired_base": "G"}``.
"""

from __future__ import annotations

import re

import pandas as pd

from beguided.web.sequence_input import VALID_EDITS

# Core targets columns produced for the pipeline (identical to the
# spreadsheet and sequence paths).
TARGETS_COLUMNS = ["chromosome", "position", "target_base", "desired_base"]

# A position must be a plain positive integer (1-based genomic coordinate).
# 12 digits is far beyond any chromosome length but keeps the value cheap
# to parse and rejects decimals, signs and empty strings by construction.
_POSITION_RE = re.compile(r"^\d{1,12}$")

# Human-readable list of the accepted editing directions (message text).
_EDIT_HINT = "Allowed: A>G, C>T, G>A, T>C."


class SingleTargetError(ValueError):
    """Raised when the single-target form fields cannot be parsed."""


def parse_single_target(chromosome: str, position: str, edit: str) -> dict:
    """Validate the three form fields and return one targets row.

    Args:
        chromosome: Chromosome name as used by the loaded reference genome
            (e.g. ``chr1``, ``1``, ``X``, ``Y``).
        position: 1-based genomic coordinate of the base to edit.
        edit: Editing direction ``X>Y`` (e.g. ``A>G``). ``X`` becomes
            ``target_base`` and ``Y`` ``desired_base``. Case-insensitive;
            whitespace around the fields (and around ``>``) is tolerated.

    Returns:
        ``{"chromosome": str, "position": int, "target_base": str,
        "desired_base": str}``.

    Raises:
        SingleTargetError: On any validation failure, with a
            human-readable message.
    """
    name = (chromosome or "").strip()
    if not name:
        raise SingleTargetError("Chromosome is required.")
    if re.search(r"\s", name):
        raise SingleTargetError(
            f"Chromosome must not contain whitespace: {name!r}. "
            "Use a single token like 'chr1', '1', 'X' or 'Y'."
        )

    raw_position = (position or "").strip()
    if not raw_position:
        raise SingleTargetError("Position is required.")
    if not _POSITION_RE.match(raw_position) or int(raw_position) < 1:
        raise SingleTargetError(
            f"Position must be a positive integer (1-based genomic "
            f"coordinate), got {position!r}."
        )

    raw_edit = (edit or "").strip().upper()
    parts = [p.strip() for p in raw_edit.split(">")]
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise SingleTargetError(
            f"Editing type must be of the form 'X>Y', got {edit!r}. "
            + _EDIT_HINT
        )
    target_base, desired_base = parts[0], parts[1]
    if (target_base, desired_base) not in VALID_EDITS:
        raise SingleTargetError(
            f"Editing direction {target_base}>{desired_base} is not "
            f"supported. {_EDIT_HINT}"
        )

    return {
        "chromosome": name,
        "position": int(raw_position),
        "target_base": target_base,
        "desired_base": desired_base,
    }


def row_to_targets_dataframe(row: dict) -> pd.DataFrame:
    """Build the standard 1-row targets DataFrame from a parsed row.

    Mirrors the spreadsheet input columns exactly so the existing
    ``run_guidefinder`` pipeline works without modification.
    """
    return pd.DataFrame(
        [
            {
                "chromosome": row["chromosome"],
                "position": row["position"],
                "target_base": row["target_base"],
                "desired_base": row["desired_base"],
            }
        ],
        columns=TARGETS_COLUMNS,
    )


__all__ = [
    "TARGETS_COLUMNS",
    "SingleTargetError",
    "parse_single_target",
    "row_to_targets_dataframe",
]
