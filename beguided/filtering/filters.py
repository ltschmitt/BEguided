"""Guide filtering logic."""

import sys
import logging

import pandas as pd
from pathlib import Path


logger = logging.getLogger("beguided.filtering")


def load_filter_criteria(config_file=None, mode=None, guide_score_threshold=None):
    """Load filter criteria from config file or command line.

    Returns:
        dict with 'mode' and 'guide_score_threshold' keys.
    """
    criteria = {
        "mode": "all",
        "guide_score_threshold": guide_score_threshold,
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
        except FileNotFoundError:
            pass

    if mode:
        criteria["mode"] = mode

    if guide_score_threshold is not None:
        criteria["guide_score_threshold"] = guide_score_threshold

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


def filter_guides(input_df, criteria):
    """Filter guides based on criteria.

    Args:
        input_df: DataFrame with guide results.
        criteria: dict with 'mode' and optional 'guide_score_threshold' and
            optional 'skip_bystander' (when True, no bystander-based filter
            is applied regardless of ``mode``).

    Returns:
        Filtered DataFrame.
    """
    df = input_df.copy()
    mode = criteria["mode"]
    skip_bystander = bool(criteria.get("skip_bystander", False))

    if not skip_bystander:
        if mode == "no_bystanders":
            if "bystander_effect" in df.columns:
                df = df[df["bystander_effect"] == "no_bystanders"]
            elif len(df) == 0:
                return df
            else:
                df = df[df["bystander_positions"].isna() | (df["bystander_positions"] == "")]

        elif mode == "silent":
            if "bystander_effect" not in df.columns:
                if len(df) == 0:
                    return df
                raise ValueError(
                    "'bystander_effect' column not found. Run bystander analysis first."
                )
            df = df[df["bystander_effect"].isin(["silent", "no_bystanders"])]

        elif mode != "all":
            raise ValueError(f"Unknown filter mode: {mode}")

    # Apply on-target score filter if configured
    threshold = criteria.get("guide_score_threshold")
    df = _apply_threshold(df, threshold)

    df = df[df["guide_rna_sequence"] != "**No guide possible**"]

    if len(df) == 0:
        return df

    return df
