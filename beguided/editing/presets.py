"""Base editor profile loading and PAM pattern management."""

import re
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from beguided.models import BaseEditorConfig


_EDITING_BASE_MAP = {
    "ABE": ("A", "G"),
    "CBE": ("C", "T"),
}

# All bundled base editor profiles use 20-nt guides. This constant is the
# single source of truth for guide length; the CRISPRater on-target score in
# beguided/guidefinder/core.py is only computed for guides of exactly this
# length. If a future profile needs a different length, expose it via a
# config key or re-introduce a TSV column.
_GUIDE_LENGTH = 20


def pam_to_regex(pam_pattern):
    """Convert PAM pattern with wildcards to regex with lookahead.

    N = A, C, T, G
    R = A, G (purines)
    Y = C, T (pyrimidines)
    W = A, T
    S = C, G
    M = A, C
    K = G, T
    B = C, G, T
    D = A, G, T
    H = A, C, T
    V = A, C, G
    """
    translation = {
        "N": "[ACGT]",
        "R": "[AG]",
        "Y": "[CT]",
        "W": "[AT]",
        "S": "[CG]",
        "M": "[AC]",
        "K": "[GT]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
    }
    pattern = "".join(
        translation.get(c.upper(), c.upper()) for c in pam_pattern
    )
    # Use lookahead to allow overlapping matches
    return f"(?=({pattern}))"


def load_base_editor_profiles(filepath):
    """Load base editor profiles from TSV file.

    Generates both sense and antisense variants for each profile.

    Args:
        filepath: Path to base editor profiles TSV.

    Returns:
        dict: {(from_base, to_base): {profile_name: BaseEditorConfig, ...}, ...}
    """
    df = pd.read_csv(filepath, sep="\t")

    name_counts = df["name"].value_counts()
    duplicates = name_counts[name_counts > 1]
    if not duplicates.empty:
        duplicate_names = ", ".join(duplicates.index.tolist())
        raise ValueError(
            f"Duplicate base editor names found in {filepath}: {duplicate_names}."
        )

    profiles = {}

    for _, row in df.iterrows():
        name = row["name"]

        editing_type = row["editing"]
        if editing_type not in _EDITING_BASE_MAP:
            continue
        from_base, to_base = _EDITING_BASE_MAP[editing_type]

        pam_pattern = row["pam"]
        pam_regex = re.compile(pam_to_regex(pam_pattern))
        guide_length = _GUIDE_LENGTH
        pam_length = len(pam_pattern)

        general_edit_start = int(row["general_edit_start"])
        general_edit_end = int(row["general_edit_end"])

        if general_edit_start < 1 or general_edit_start > guide_length:
            raise ValueError(
                f"Profile '{name}' has general_edit_start={general_edit_start} outside valid range [1, {guide_length}]"
            )
        if general_edit_end < 1 or general_edit_end > guide_length:
            raise ValueError(
                f"Profile '{name}' has general_edit_end={general_edit_end} outside valid range [1, {guide_length}]"
            )
        if general_edit_start > general_edit_end:
            raise ValueError(
                f"Profile '{name}' has general_edit_start={general_edit_start} > general_edit_end={general_edit_end}"
            )

        sense_key = (from_base, to_base)
        if sense_key not in profiles:
            profiles[sense_key] = {}

        profiles[sense_key][name] = BaseEditorConfig(
            name=name,
            pam_pattern=pam_pattern,
            pam_regex=pam_regex,
            pam_length=pam_length,
            guide_length=guide_length,
            high_edit_start=int(row["high_edit_start"]),
            high_edit_end=int(row["high_edit_end"]),
            general_edit_start=general_edit_start,
            general_edit_end=general_edit_end,
            strand="+",
            reverse_complement=False,
            target_base=from_base,
            desired_base=to_base,
        )

        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        antisense_from = complement[from_base]
        antisense_to = complement[to_base]
        antisense_pam = str(Seq(pam_pattern).reverse_complement())
        antisense_pam_regex = re.compile(pam_to_regex(antisense_pam))

        antisense_key = (antisense_from, antisense_to)
        if antisense_key not in profiles:
            profiles[antisense_key] = {}

        profiles[antisense_key][name] = BaseEditorConfig(
            name=name,
            pam_pattern=antisense_pam,
            pam_regex=antisense_pam_regex,
            pam_length=len(antisense_pam),
            guide_length=guide_length,
            high_edit_start=int(row["high_edit_start"]),
            high_edit_end=int(row["high_edit_end"]),
            general_edit_start=general_edit_start,
            general_edit_end=general_edit_end,
            strand="-",
            reverse_complement=True,
            target_base=antisense_from,
            desired_base=antisense_to,
        )

    return profiles


def get_available_base_editor_names(profiles):
    """Return a set of all unique base editor names across all editing directions."""
    available = set()
    for editor_profiles in profiles.values():
        available.update(editor_profiles.keys())
    return available


def get_base_editor_names_by_type(profiles):
    """Return ABE and CBE base editor names from the loaded profiles dict.

    Classification uses the sense-strand (``+") editing direction:
    A→G  -> ABE,  C→T -> CBE.

    Args:
        profiles: dict from load_base_editor_profiles.

    Returns:
        dict: {"abe": [sorted names], "cbe": [sorted names]}
    """
    abe = set()
    cbe = set()
    for editor_key, editor_profiles in profiles.items():
        for name, config in editor_profiles.items():
            if config.strand == "+":
                if config.target_base == "A" and config.desired_base == "G":
                    abe.add(name)
                elif config.target_base == "C" and config.desired_base == "T":
                    cbe.add(name)
    return {"abe": sorted(abe), "cbe": sorted(cbe)}


def get_supported_editor_keys(profiles, base_editor_names):
    """Return the set of (target_base, desired_base) keys supported by the
    selected base editor names, including antisense complements automatically.

    Args:
        profiles: dict from load_base_editor_profiles.
        base_editor_names: list of base editor name strings.

    Returns:
        set of (target_base, desired_base) tuples.
    """
    supported = set()
    for editor_key, editor_profiles in profiles.items():
        for name in base_editor_names:
            if name in editor_profiles:
                supported.add(editor_key)
                break
    return supported


def filter_compatible_targets(targets_df, profiles, base_editor_names):
    """Filter targets DataFrame to only rows compatible with selected base editors.

    Complement editing directions (e.g. T->C for an ABE A->G profile) are
    automatically accepted because load_base_editor_profiles generates
    antisense variants.

    Args:
        targets_df: pandas DataFrame with target_base and desired_base cols.
        profiles: dict from load_base_editor_profiles.
        base_editor_names: list of base editor name strings.

    Returns:
        tuple: (compatible_df, incompatible_df)
    """
    supported = get_supported_editor_keys(profiles, base_editor_names)
    mask = targets_df.apply(
        lambda row: (str(row["target_base"]), str(row["desired_base"])) in supported,
        axis=1,
    )
    return targets_df[mask].copy(), targets_df[~mask].copy()
