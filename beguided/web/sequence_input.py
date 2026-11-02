"""Parser for the sequence input mode.

The sequence input mode accepts raw sequence text (one or many FASTA records)
where each sequence contains a single ``[X>Y]`` bracket marking the desired
edit position. The parser produces:

  * A targets ``DataFrame`` with the same shape as the spreadsheet path so the
    existing pipeline (``run_guidefinder`` → ``run_pipeline``) can consume it
    without changes.
  * A custom-genome ``dict[str, str]`` / FASTA string that overrides the
    operator-configured genome FASTA, because the user's custom sequences have
    no relation to any installed reference genome.

Position conventions:
  * ``edit_position`` is 1-based and points at the base represented by the
    ``[X>Y]`` bracket. The bracket is removed from the sequence and ``X``
    (the target_base) is inserted at that position; the text originally after
    ``]`` shifts one position later in the resulting sequence. This matches
    the convention the rest of the pipeline expects (1-based genomic
    positions).

Example:
    Input:  ``ATAG[T>C]GAC``
    Edit position: 5 (1-based, where ``[`` sat)
    Stripped sequence: ``ATAGTGAC`` (8 nt)
    Position 5: ``T`` (the bracket's target_base)
    Position 6: ``G`` (was originally right after ``]``)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, List

import pandas as pd


# A bracket of the form [X>Y] where X and Y are A/C/G/T (uppercase or lowercase).
# Whitespace inside the brackets is tolerated.
_BRACKET_RE = re.compile(r"\[\s*([ACGTacgt])\s*>\s*([ACGTacgt])\s*\]")

# Accepted editing directions. The pipeline's base editor profile loader
# already creates sense + antisense complement entries automatically (see
# editing/presets.py), so all four user-friendly spellings are valid input.
# Public so other input modes (e.g. the single-target web input) can reuse
# the same set.
VALID_EDITS = {
    ("A", "G"),
    ("T", "C"),
    ("C", "T"),
    ("G", "A"),
}

# Hard cap on the length of a single sequence (after bracket removal).
# 1000 bp is more than enough for the editing windows supported by all base editor profiles
# (max ~30 nt around the edit). Longer sequences waste server resources and
# rarely yield more guides.
MAX_SEQUENCE_LENGTH = 1000

# Minimum length after bracket removal. Anything shorter cannot host a 20-nt
# guide + PAM and the guidefinder will report "no guide possible" anyway;
# rejecting upfront gives a clearer error message.
MIN_SEQUENCE_LENGTH = 20

DEFAULT_CHROMOSOME = "custom"


class SequenceParseError(ValueError):
    """Raised when the sequence text cannot be parsed.

    The ``record_index`` field points at the offending sequence (0-based) so
    callers can surface ``"Sequence #3: ..."`` style messages.
    """

    def __init__(self, message: str, record_index: int = -1) -> None:
        super().__init__(message)
        self.record_index = record_index


@dataclass
class SequenceRecord:
    """A single parsed sequence with its edit annotation."""

    chromosome: str
    sequence: str
    edit_position: int  # 1-based, position of the bracket's `[` + 1
    target_base: str
    desired_base: str
    genome_base: str  # base at edit_position in the stripped sequence
                      # (equal to target_base by construction)
    raw_text: str = field(repr=False)


def _parse_bracket(text: str) -> tuple[int, str, str] | None:
    """Find the first ``[X>Y]`` bracket in ``text``.

    Returns ``(start_index, target_base, desired_base)`` or ``None`` if not
    found. ``start_index`` is the absolute index in ``text`` where the ``[``
    sits, so callers can compute the 1-based edit position by counting
    non-bracket characters before it.
    """
    m = _BRACKET_RE.search(text)
    if m is None:
        return None
    target = m.group(1).upper()
    desired = m.group(2).upper()
    return m.start(), target, desired


def _strip_brackets(text: str) -> str:
    """Remove every ``[X>Y]`` substring from ``text``.

    Used to obtain the bare DNA sequence (the "genome") after the edit
    annotation has been extracted.
    """
    return _BRACKET_RE.sub("", text)


def _sanitize_chromosome_name(raw: str, record_index: int) -> str:
    """Return a chromosome id safe for downstream TSV/CSV output.

    Mirrors the convention in ``beguided.genome.io.load_genome``: take the
    first whitespace-separated token from the FASTA header. Reject names
    containing tab or newline characters because they would corrupt TSV
    output.
    """
    name = raw.strip()
    if not name:
        raise SequenceParseError(
            "FASTA header is empty; cannot derive chromosome name.",
            record_index=record_index,
        )
    if "\t" in name or "\n" in name or "\r" in name:
        raise SequenceParseError(
            f"FASTA header contains tab/newline characters: {name!r}.",
            record_index=record_index,
        )
    return name.split()[0]


def _parse_records(text: str) -> List[tuple[str | None, str]]:
    """Parse raw text into ``[(header, sequence), ...]``.

    Header semantics:
      * ``""``   → plain single sequence with no FASTA header at all
                    (use the default chromosome name).
      * ``None`` → user wrote a ``>`` line but left the name empty
                    (treated as an error in the caller).
      * non-empty string → the FASTA header text (first whitespace token
                    becomes the chromosome name).

    Empty lines and whitespace-only lines are ignored.
    """
    if text is None:
        raise SequenceParseError("No sequence text supplied.")

    raw = text.replace("\r\n", "\n").replace("\r", "\n")
    if not raw.strip():
        raise SequenceParseError("No sequence text supplied.")

    lines = [ln for ln in raw.split("\n") if ln.strip() != ""]
    has_fasta_header = any(ln.startswith(">") for ln in lines)

    if has_fasta_header:
        records: list[tuple[str | None, str]] = []
        current_header: str | None = None
        current_chunks: list[str] = []
        for ln in lines:
            if ln.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_chunks)))
                # ``None`` marks an explicit empty header (">" with no name);
                # ``""`` would be a plain sequence.
                name = ln[1:].strip()
                current_header = name if name else None
                current_chunks = []
            else:
                current_chunks.append(ln.strip())
        if current_header is not None or current_chunks:
            records.append((current_header, "".join(current_chunks)))
        if not records:
            raise SequenceParseError("FASTA input did not contain any sequences.")
        return records

    # Plain single sequence (no FASTA header). Empty header → "" → caller
    # uses the default chromosome name.
    sequence = "".join(ln.strip() for ln in lines)
    return [("", sequence)]


def _validate_dna(seq: str, record_index: int) -> None:
    """Reject any character outside ``A/C/G/T`` (case-insensitive).

    IUPAC ambiguity codes (``N``, ``R``, ...) are rejected here so the
    guidefinder gets a clean sequence. This matches the ``chromosome``
    validation stance taken elsewhere in the web layer.
    """
    cleaned = re.sub(r"[\sACGTacgt]", "", seq)
    if cleaned:
        bad = sorted({c for c in cleaned})[:5]
        raise SequenceParseError(
            f"Sequence contains non-ACGT characters: {''.join(bad)!r}.",
            record_index=record_index,
        )


def parse_sequence_text(
    text: str, *, default_chromosome: str = DEFAULT_CHROMOSOME
) -> List[SequenceRecord]:
    """Parse raw sequence text into a list of ``SequenceRecord``.

    Args:
        text: Raw user input. Plain sequence or FASTA.
        default_chromosome: Chromosome name used when a plain sequence has no
            FASTA header. Defaults to ``"custom"``.

    Returns:
        List of ``SequenceRecord`` (one per input sequence), each carrying the
        chromosome name, the bracket-stripped DNA sequence, and the edit
        annotation.

    Raises:
        SequenceParseError: On any parse failure, with a human-readable
            message and the offending record's index.
    """
    raw_records = _parse_records(text)
    out: list[SequenceRecord] = []
    seen_chromosomes: set[str] = set()

    for idx, (header, sequence) in enumerate(raw_records):
        # Header validation first — failing here is a structural problem
        # with the input that's easier to fix than a missing bracket.
        if header is None:
            raise SequenceParseError(
                "FASTA header is empty (a '>' line with no name). "
                "Provide a header or omit the '>' entirely to use the "
                "default chromosome name.",
                record_index=idx,
            )
        if header:
            chromosome = _sanitize_chromosome_name(header, idx)
        else:
            chromosome = default_chromosome

        if chromosome in seen_chromosomes:
            raise SequenceParseError(
                f"Duplicate sequence identifier {chromosome!r}: each input "
                "sequence must have a unique header (FASTA '>name') so it "
                "can be tracked individually through the pipeline.",
                record_index=idx,
            )
        seen_chromosomes.add(chromosome)

        # Reject any extra brackets: each sequence must contain exactly one.
        bracket_matches = list(_BRACKET_RE.finditer(sequence))
        if len(bracket_matches) == 0:
            raise SequenceParseError(
                "Sequence is missing an edit bracket. "
                "Use the syntax [A>G] (or [T>C], [C>T], [G>A]) "
                "to mark the position of the desired edit.",
                record_index=idx,
            )
        if len(bracket_matches) > 1:
            raise SequenceParseError(
                "Sequence contains multiple edit brackets; "
                "exactly one [X>Y] per sequence is allowed.",
                record_index=idx,
            )

        bracket_match = bracket_matches[0]
        target_base = bracket_match.group(1).upper()
        desired_base = bracket_match.group(2).upper()

        if (target_base, desired_base) not in VALID_EDITS:
            raise SequenceParseError(
                f"Edit direction [{target_base}>{desired_base}] is not supported. "
                "Allowed: [A>G], [T>C], [C>T], [G>A].",
                record_index=idx,
            )

        # 1-based position: the bracket's `[` index + 1. After bracket
        # removal + target_base insertion, this is the position of the new
        # base in the stripped sequence.
        edit_position = bracket_match.start() + 1

        # Build the stripped sequence by replacing the bracket with the
        # target_base. The text originally after `]` shifts one position
        # later. This means the bracket IS the base at the edit position:
        # whatever the user typed on the left of `>` becomes the genome
        # base at edit_position in the stripped sequence.
        stripped = (
            sequence[:bracket_match.start()]
            + target_base
            + sequence[bracket_match.end():]
        ).upper()

        if len(stripped) < MIN_SEQUENCE_LENGTH:
            raise SequenceParseError(
                f"Sequence is too short ({len(stripped)} nt after "
                f"inserting the target base); need at least "
                f"{MIN_SEQUENCE_LENGTH} nt.",
                record_index=idx,
            )
        if len(stripped) > MAX_SEQUENCE_LENGTH:
            raise SequenceParseError(
                f"Sequence is too long ({len(stripped)} nt); "
                f"cap is {MAX_SEQUENCE_LENGTH} nt.",
                record_index=idx,
            )

        _validate_dna(stripped, idx)

        genome_base = stripped[edit_position - 1]
        # By construction ``genome_base == target_base`` — the bracket was
        # inserted at the edit position. No mismatch check needed.

        out.append(
            SequenceRecord(
                chromosome=chromosome,
                sequence=stripped,
                edit_position=edit_position,
                target_base=target_base,
                desired_base=desired_base,
                genome_base=genome_base,
                raw_text=sequence,
            )
        )

    return out


def records_to_targets_dataframe(records: Iterable[SequenceRecord]) -> pd.DataFrame:
    """Convert ``SequenceRecord``s into the standard targets DataFrame.

    Mirrors the spreadsheet input columns exactly so the existing
    ``run_guidefinder`` pipeline works without modification:

        chromosome, position, target_base, desired_base

    Extra columns (``sequence``, ``sequence_preview``) are included for the
    web UI preview but excluded from ``passthrough_extras`` so they don't
    leak into the final TSV.
    """
    rows = []
    for r in records:
        rows.append({
            "chromosome": r.chromosome,
            "position": r.edit_position,
            "target_base": r.target_base,
            "desired_base": r.desired_base,
            "sequence": r.sequence,
            "sequence_preview": _make_preview(r.sequence, r.edit_position),
        })
    df = pd.DataFrame(
        rows,
        columns=[
            "chromosome",
            "position",
            "target_base",
            "desired_base",
            "sequence",
            "sequence_preview",
        ],
    )
    return df


def records_to_fasta(records: Iterable[SequenceRecord]) -> str:
    """Serialize records to a single multi-FASTA string.

    The custom genome FASTA mirrors the operator-configured FASTA convention
    (``beguided.genome.io.load_genome``): the chromosome id is the first
    whitespace-separated token of the header.
    """
    parts: list[str] = []
    for r in records:
        parts.append(f">{r.chromosome}\n{r.sequence}\n")
    return "".join(parts)


def _make_preview(sequence: str, edit_position: int, window: int = 15) -> str:
    """Build a short preview string centred on the edit position.

    Returns a snippet like ``...ATAGACCCT[C]CCTATCC...`` with the bracketed
    base wrapped in ``[...]``. ``window`` controls how many bases are shown
    on each side.
    """
    idx = edit_position - 1
    start = max(0, idx - window)
    end = min(len(sequence), idx + window + 1)
    prefix = "..." if start > 0 else ""
    suffix = "..." if end < len(sequence) else ""
    left = sequence[start:idx]
    base = sequence[idx]
    right = sequence[idx + 1:end]
    return f"{prefix}{left}[{base}]{right}{suffix}"


__all__ = [
    "DEFAULT_CHROMOSOME",
    "MAX_SEQUENCE_LENGTH",
    "MIN_SEQUENCE_LENGTH",
    "VALID_EDITS",
    "SequenceParseError",
    "SequenceRecord",
    "parse_sequence_text",
    "records_to_fasta",
    "records_to_targets_dataframe",
]
