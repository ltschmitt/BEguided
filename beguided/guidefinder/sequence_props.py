"""Sequence property helpers used by the guidefinder.

These functions are deliberately tiny and dependency-free so they can be
unit-tested in isolation and reused outside the main guidefinder pipeline
(for example in the web UI).
"""

import re


_GC_PATTERN = re.compile("[CG]")
_RUN_PATTERN = re.compile(r"(.)\1+")


def gc_content(seq: str) -> float:
    """Return the GC content of *seq* as a percentage in the range 0..100.

    Only A/C/G/T characters are considered: ambiguous bases (N, R, Y, ...)
    are ignored entirely rather than treated as half-GC. The result is
    rounded to one decimal place so the column is stable for equality
    comparisons in tests and so the web UI can display it without further
    rounding.
    """
    if not seq:
        return 0.0
    cg = len(_GC_PATTERN.findall(seq))
    return round(100.0 * cg / len(seq), 1)


def longest_homopolymer_run(seq: str) -> int:
    """Return the length of the longest run of identical bases in *seq*.

    A return value of 1 means every base differs from its neighbour
    (no homopolymer at all). For an empty string the return value is 0.
    """
    if not seq:
        return 0
    longest = 1
    for match in _RUN_PATTERN.finditer(seq):
        longest = max(longest, len(match.group(0)))
    return longest


def has_homopolymer(seq: str, min_run: int = 4) -> bool:
    """Return True if *seq* contains any run of identical bases of
    length >= *min_run* (default 4, the BEguided cutoff for an
    "unwanted" homopolymer)."""
    return longest_homopolymer_run(seq) >= min_run
