"""Guide filtering logic."""

import sys
import logging

import pandas as pd
from pathlib import Path

from beguided.models import (
    GUIDES_REMOVED_BY_FILTER,
    PLACEHOLDER_MARKERS,
)


logger = logging.getLogger("beguided.filtering")


VALID_GC_FILTERS = ("off", "high", "medium")
_GC_RANGES = {
    "high": (40.0, 60.0),
    "medium": (30.0, 70.0),
}
# Deprecated pre-rename names, kept working for existing configs and API
# callers. ``good`` and ``acceptable`` had the same GC ranges as ``high``
# and ``medium`` respectively.
_GC_ALIASES = {
    "good": "high",
    "acceptable": "medium",
}


def normalize_gc_filter(mode):
    """Map a ``gc_content_filter`` value to its canonical name.

    Accepts the canonical values (``"off"``, ``"high"``, ``"medium"``) and
    the deprecated aliases (``"good"`` → ``"high"``, ``"acceptable"`` →
    ``"medium"``, logged as deprecated).

    Raises:
        ValueError: for any other value.
    """
    if mode is None or mode == "":
        return "off"
    if mode in _GC_ALIASES:
        canonical = _GC_ALIASES[mode]
        logger.warning(
            "gc_content_filter=%r is deprecated; use %r instead (same GC range).",
            mode, canonical,
        )
        return canonical
    if mode in VALID_GC_FILTERS:
        return mode
    raise ValueError(
        f"Invalid gc_content_filter: {mode!r}. "
        f"Expected one of: {', '.join(VALID_GC_FILTERS)}"
        f" (deprecated aliases: {', '.join(sorted(_GC_ALIASES))})."
    )


def load_filter_criteria(
    config_file=None,
    mode=None,
    guide_score_threshold=None,
    gc_content_filter=None,
    filter_homopolymers=None,
    trinuc_edit_score_min=None,
):
    """Load filter criteria from config file or command line.

    Returns:
        dict with 'mode', 'guide_score_threshold', 'gc_content_filter',
        'filter_homopolymers' and 'trinuc_edit_score_min' keys. The new keys
        default to "off", False and None respectively so legacy callers see
        no behaviour change.
    """
    criteria = {
        "mode": "all",
        "guide_score_threshold": guide_score_threshold,
        "gc_content_filter": "off",
        "filter_homopolymers": False,
        "trinuc_edit_score_min": None,
    }

    if config_file:
        try:
            import yaml
            with open(config_file, "r") as f:
                config = yaml.safe_load(f)
            if config:
                if "keep_bystanders" in config:
                    criteria["mode"] = config["keep_bystanders"]
                elif "filtering" in config:
                    criteria["mode"] = config["filtering"].get("mode", criteria["mode"])

                # guide_score_threshold can live at root or under filtering
                if "guide_score_threshold" in config:
                    criteria["guide_score_threshold"] = config["guide_score_threshold"]
                elif "filtering" in config and "guide_score_threshold" in config["filtering"]:
                    criteria["guide_score_threshold"] = config["filtering"]["guide_score_threshold"]

                if "gc_content_filter" in config:
                    criteria["gc_content_filter"] = config["gc_content_filter"]
                elif "filtering" in config and "gc_content_filter" in config["filtering"]:
                    criteria["gc_content_filter"] = config["filtering"]["gc_content_filter"]

                if "filter_homopolymers" in config:
                    criteria["filter_homopolymers"] = bool(config["filter_homopolymers"])
                elif "filtering" in config and "filter_homopolymers" in config["filtering"]:
                    criteria["filter_homopolymers"] = bool(config["filtering"]["filter_homopolymers"])

                if "trinuc_edit_score_min" in config:
                    raw = config["trinuc_edit_score_min"]
                    criteria["trinuc_edit_score_min"] = (
                        None if raw in (None, "", "null", "None") else float(raw)
                    )
                elif "filtering" in config and "trinuc_edit_score_min" in config["filtering"]:
                    raw = config["filtering"]["trinuc_edit_score_min"]
                    criteria["trinuc_edit_score_min"] = (
                        None if raw in (None, "", "null", "None") else float(raw)
                    )
        except FileNotFoundError:
            pass

    if mode:
        criteria["mode"] = mode

    if guide_score_threshold is not None:
        criteria["guide_score_threshold"] = guide_score_threshold

    if gc_content_filter is not None:
        criteria["gc_content_filter"] = normalize_gc_filter(gc_content_filter)

    if filter_homopolymers is not None:
        criteria["filter_homopolymers"] = bool(filter_homopolymers)

    if trinuc_edit_score_min is not None:
        try:
            trinuc_min = float(trinuc_edit_score_min)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"trinuc_edit_score_min must be a number, "
                f"got {trinuc_edit_score_min!r}"
            ) from exc
        if not (0.0 <= trinuc_min <= 1.0):
            raise ValueError(
                f"trinuc_edit_score_min must be between 0.0 and 1.0, "
                f"got {trinuc_min}"
            )
        criteria["trinuc_edit_score_min"] = trinuc_min

    return criteria


def _apply_threshold(df, threshold):
    """Remove guides with guide_score < threshold.

    Guides without a score (NaN) are kept so that non-20-nt base editors or
    pre-existing data are not silently discarded.
    """
    if "guide_score" not in df.columns:
        logger.info(
            "'guide_score' column not found; skipping on-target score filter"
        )
        return df

    if threshold is None:
        return df

    mask = df["guide_score"].isna() | (df["guide_score"] >= threshold)
    removed = (~mask).sum()
    if removed > 0:
        logger.info(
            "On-target score filter (score >= %.2f) removed %d guide(s)",
            threshold, int(removed)
        )
    return df[mask]


def _apply_gc_filter(df, mode):
    """Restrict guides by GC content (percent 0-100).

    ``mode`` is one of ``"off"`` (no filter, the default), ``"high"`` (keep
    40-60%) or ``"medium"`` (keep 30-70%). The deprecated aliases
    ``"good"``/``"acceptable"`` are accepted and mapped to
    ``"high"``/``"medium"``. Rows with a missing ``gc_content`` value are
    kept, matching the guide-score policy of favouring recall over silent
    dropouts.
    """
    mode = normalize_gc_filter(mode)
    if mode == "off":
        return df

    if "gc_content" not in df.columns:
        logger.info(
            "'gc_content' column not found; skipping GC content filter"
        )
        return df

    lo, hi = _GC_RANGES[mode]
    mask = df["gc_content"].isna() | ((df["gc_content"] >= lo) & (df["gc_content"] <= hi))
    removed = (~mask).sum()
    if removed > 0:
        logger.info(
            "GC content filter (%.0f-%.0f%%, mode=%r) removed %d guide(s)",
            lo, hi, mode, int(removed)
        )
    return df[mask]


def _apply_homopolymer_filter(df, enabled):
    """Drop guides whose ``homopolymer`` flag is True when *enabled*."""
    if not enabled:
        return df

    if "homopolymer" not in df.columns:
        logger.info(
            "'homopolymer' column not found; skipping homopolymer filter"
        )
        return df

    mask = df["homopolymer"].isna() | (df["homopolymer"] == False)  # noqa: E712
    removed = (~mask).sum()
    if removed > 0:
        logger.info(
            "Homopolymer filter (>=4-mer) removed %d guide(s)",
            int(removed)
        )
    return df[mask]


def _apply_trinuc_score_filter(df, threshold):
    """Drop ABE guides whose ``trinuc_edit_score`` is below *threshold*.

    The trinucleotide context score is the Marquart et al. 2021 normalised
    ABE8e proportion (0..1). CBE rows have ``NaN`` (the table is ABE-only)
    and are always kept, matching the policy used for ``gc_content`` and
    ``homopolymer``: missing values never silently drop a guide. ``None``
    disables the filter entirely.
    """
    if threshold is None:
        return df

    if "trinuc_edit_score" not in df.columns:
        logger.info(
            "'trinuc_edit_score' column not found; skipping trinucleotide score filter"
        )
        return df

    mask = df["trinuc_edit_score"].isna() | (df["trinuc_edit_score"] >= threshold)
    removed = (~mask).sum()
    if removed > 0:
        logger.info(
            "Trinucleotide edit-score filter (>= %.3f, ABE only) removed %d guide(s)",
            threshold, int(removed),
        )
    return df[mask]


def filter_guides(input_df, criteria, fill_missing_targets=False):
    """Filter guides based on criteria.

    Args:
        input_df: DataFrame with guide results.
        criteria: dict with 'mode' and optional 'guide_score_threshold' and
            optional 'skip_bystander' (when True, no bystander-based filter
            is applied regardless of ``mode``).
        fill_missing_targets: When True, every target present in ``input_df``
            that has no surviving row after filtering gets one placeholder row
            at its input position (identity + bases columns copied from the
            input, other result fields empty); rows keep the target order of
            ``input_df``. The ``guide_rna_sequence`` marker is
            preserved when the input row was already a placeholder (e.g.
            ``**No guide possible**``) and set to
            ``**Guides removed by filter**`` when the input row was a real
            guide removed by the filters. Used by the web UI's dynamic
            filters so a target whose guides are all removed still shows up.
            Off by default (CLI pipeline step 3 behaviour is unchanged).

    Returns:
        Filtered DataFrame.
    """
    df = input_df.copy()
    mode = criteria["mode"]
    skip_bystander = bool(criteria.get("skip_bystander", False))

    if not skip_bystander:
        if mode == "none":
            if "bystander_effect" in df.columns:
                df = df[df["bystander_effect"] == "none"]
            elif "bystander_positions" in df.columns:
                df = df[df["bystander_positions"].isna() | (df["bystander_positions"] == "")]
            # No bystander metadata at all: nothing to filter on, rows pass.

        elif mode == "synonymous":
            if "bystander_effect" not in df.columns:
                # No bystander metadata at all (e.g. sequence-mode jobs where
                # step 2 is skipped): nothing to filter on. Fall back to the
                # ``bystander_positions`` column when present, mirroring the
                # "none" branch, and otherwise let every row pass.
                if "bystander_positions" in df.columns:
                    df = df[
                        df["bystander_positions"].isna()
                        | (df["bystander_positions"] == "")
                    ]
            else:
                effect = df["bystander_effect"]
                # Keep synonymous and bystander-free guides, plus rows with a
                # missing/empty annotation (e.g. sequence mode). This favours
                # recall over silent dropouts, like the other filters.
                df = df[
                    effect.isin(["synonymous", "none"])
                    | effect.isna()
                    | (effect == "")
                ]

        elif mode != "all":
            raise ValueError(f"Unknown filter mode: {mode}")

    # Apply on-target score filter if configured
    threshold = criteria.get("guide_score_threshold")
    df = _apply_threshold(df, threshold)

    # Apply GC content tier filter if configured
    gc_mode = criteria.get("gc_content_filter", "off")
    df = _apply_gc_filter(df, gc_mode)

    # Apply homopolymer filter if configured
    df = _apply_homopolymer_filter(df, criteria.get("filter_homopolymers", False))

    # Apply trinucleotide edit-score filter if configured (ABE only; CBE rows kept)
    df = _apply_trinuc_score_filter(df, criteria.get("trinuc_edit_score_min"))

    # Remove the pipeline's internal "not a real guide" marker rows; these are
    # re-added by the placeholder fill below (with their marker preserved) so
    # the web UI still shows one row per target.
    if "guide_rna_sequence" in df.columns:
        df = df[df["guide_rna_sequence"] != "**No guide possible**"]

    # Historic short-circuit: an empty frame with no matching data columns is
    # returned as-is unless the placeholder fill is requested (the fill can
    # still produce rows from a non-empty source frame).
    if len(df) == 0 and not fill_missing_targets:
        return df

    if fill_missing_targets:
        df = _fill_missing_target_rows(df, input_df)

    return df


# Columns that identify a target row. These are copied verbatim from the
# source row when building a placeholder; every other result column is
# emptied. The list covers spreadsheet mode and sequence mode alike
# (chromosome/position, target_base/desired_base); ``sequence`` and
# ``desired_edit`` are legacy display aliases kept as harmless fallbacks.
_PLACEHOLDER_IDENTITY_COLUMNS = (
    "chromosome",
    "sequence",
    "position",
    "target_base",
    "desired_base",
    "desired_edit",
)


def _fill_missing_target_rows(df, source_df):
    """Interleave one placeholder row per target that lost all its guides.

    ``source_df`` is the pre-filter frame and ``df`` the already-filtered
    frame. The returned frame keeps every surviving row of ``df`` and follows
    the target order of ``source_df`` (which is the pipeline's input order):
    each target keeps its surviving rows, and a target with no surviving row
    gets one placeholder row at its input position:

    * identity columns (chromosome/sequence, position, target/desired bases)
      are copied from the first source row for that target;
    * the ``guide_rna_sequence`` marker is preserved when the source row was
      already a placeholder (e.g. ``**No guide possible**`` or
      ``**Base change incompatible**``) and set to ``**Guides removed by
      filter**`` when the source row was a real guide removed by the filters;
    * every other result column is left empty/``None``.

    The target key matches the frontend's ``targetKey`` helper: the
    chromosome (or sequence in sequence mode) plus the 1-based position
    (empty when the column is absent).
    """
    if "guide_rna_sequence" not in source_df.columns:
        return df

    def _target_key(row):
        chrom = row.get("chromosome") or row.get("sequence")
        pos = row.get("position")
        return ("" if chrom is None else str(chrom), "" if pos is None else str(pos))

    # First-appearance order of every target in the pre-filter frame (which
    # follows the operator's input order), plus the first row per target as
    # the placeholder template.
    order: dict = {}
    templates: dict = {}
    for _, row in source_df.iterrows():
        key = _target_key(row)
        if key not in order:
            order[key] = len(order)
            templates[key] = row

    # Surviving rows inherit their target's input position; keys missing from
    # the source frame (should not happen) sort last, order preserved.
    row_keys = [_target_key(row) for _, row in df.iterrows()]
    present = set(row_keys)
    fallback = len(order)
    df = df.assign(_target_order=[order.get(k, fallback) for k in row_keys])

    placeholder_rows = []
    for key, template in templates.items():
        if key in present:
            continue
        row = {}
        for col in df.columns:
            if col in _PLACEHOLDER_IDENTITY_COLUMNS and col in template.index:
                row[col] = template[col]
            else:
                row[col] = None
        source_seq = template.get("guide_rna_sequence")
        if source_seq in PLACEHOLDER_MARKERS:
            row["guide_rna_sequence"] = source_seq
        else:
            row["guide_rna_sequence"] = GUIDES_REMOVED_BY_FILTER
        row["_target_order"] = order[key]
        placeholder_rows.append(row)

    if not placeholder_rows:
        return df.drop(columns="_target_order")

    new_rows = pd.DataFrame(placeholder_rows, columns=df.columns)
    merged = pd.concat([df, new_rows], ignore_index=True)
    merged = merged.sort_values("_target_order", kind="stable")
    return merged.drop(columns="_target_order").reset_index(drop=True)
