"""Preset loading and PAM pattern management."""

import re
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from beguided.models import PresetConfig


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


def load_presets(filepath):
    """Load base editor presets from TSV file.

    Generates both sense and antisense variants for each preset.

    Args:
        filepath: Path to presets TSV.

    Returns:
        dict: {(from_base, to_base): {preset_name: PresetConfig, ...}, ...}
    """
    df = pd.read_csv(filepath, sep="\t")

    # Check for duplicate preset names
    name_counts = df["name"].value_counts()
    duplicates = name_counts[name_counts > 1]
    if not duplicates.empty:
        duplicate_names = ", ".join(duplicates.index.tolist())
        raise ValueError(
            f"Duplicate preset names found in {filepath}: {duplicate_names}."
        )

    presets = {}

    for _, row in df.iterrows():
        name = row["name"]

        # Parse editing direction (e.g., "A-to-G")
        editing_parts = row["editing"].split("-to-")
        if len(editing_parts) != 2:
            continue
        from_base, to_base = editing_parts[0], editing_parts[1]

        pam_pattern = row["pam"]
        pam_regex = re.compile(pam_to_regex(pam_pattern))
        guide_length = int(row["guide_length"])
        pam_length = len(pam_pattern)

        general_edit_start = int(row["general_edit_start"])
        general_edit_end = int(row["general_edit_end"])

        if general_edit_start < 1 or general_edit_start > guide_length:
            raise ValueError(
                f"Preset '{name}' has general_edit_start={general_edit_start} outside valid range [1, {guide_length}]"
            )
        if general_edit_end < 1 or general_edit_end > guide_length:
            raise ValueError(
                f"Preset '{name}' has general_edit_end={general_edit_end} outside valid range [1, {guide_length}]"
            )
        if general_edit_start > general_edit_end:
            raise ValueError(
                f"Preset '{name}' has general_edit_start={general_edit_start} > general_edit_end={general_edit_end}"
            )

        # Sense strand
        sense_key = (from_base, to_base)
        if sense_key not in presets:
            presets[sense_key] = {}

        presets[sense_key][name] = PresetConfig(
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

        # Antisense strand (complement)
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        antisense_from = complement[from_base]
        antisense_to = complement[to_base]
        antisense_pam = str(Seq(pam_pattern).reverse_complement())
        antisense_pam_regex = re.compile(pam_to_regex(antisense_pam))

        antisense_key = (antisense_from, antisense_to)
        if antisense_key not in presets:
            presets[antisense_key] = {}

        presets[antisense_key][name] = PresetConfig(
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

    return presets


def get_available_preset_names(presets):
    """Return a set of all unique preset names across all editing directions."""
    available = set()
    for editor_presets in presets.values():
        available.update(editor_presets.keys())
    return available


def get_preset_names_by_type(presets):
    """Return ABE and CBE preset names from the loaded presets dict.

    Classification uses the sense-strand (``+") editing direction:
    A→G  -> ABE,  C→T -> CBE.

    Args:
        presets: dict from load_presets.

    Returns:
        dict: {"abe": [sorted names], "cbe": [sorted names]}
    """
    abe = set()
    cbe = set()
    for editor_key, editor_presets in presets.items():
        for name, config in editor_presets.items():
            if config.strand == "+":
                if config.target_base == "A" and config.desired_base == "G":
                    abe.add(name)
                elif config.target_base == "C" and config.desired_base == "T":
                    cbe.add(name)
    return {"abe": sorted(abe), "cbe": sorted(cbe)}


def get_supported_editor_keys(presets, preset_names):
    """Return the set of (target_base, desired_base) keys supported by the
    selected preset names, including antisense complements automatically.

    Args:
        presets: dict from load_presets.
        preset_names: list of preset name strings.

    Returns:
        set of (target_base, desired_base) tuples.
    """
    supported = set()
    for editor_key, editor_presets in presets.items():
        for name in preset_names:
            if name in editor_presets:
                supported.add(editor_key)
                break
    return supported


def filter_compatible_targets(targets_df, presets, preset_names):
    """Filter targets DataFrame to only rows compatible with selected presets.

    Complement editing directions (e.g. T->C for an ABE A->G preset) are
    automatically accepted because load_presets generates antisense variants.

    Args:
        targets_df: pandas DataFrame with target_base and desired_base cols.
        presets: dict from load_presets.
        preset_names: list of preset name strings.

    Returns:
        tuple: (compatible_df, incompatible_df)
    """
    supported = get_supported_editor_keys(presets, preset_names)
    mask = targets_df.apply(
        lambda row: (str(row["target_base"]), str(row["desired_base"])) in supported,
        axis=1,
    )
    return targets_df[mask].copy(), targets_df[~mask].copy()
