"""Merge off-target summaries with guide data."""

import pandas as pd


OTB_COLUMNS = [
    "guide_seq",
    "worst_offtarget_score",
    "worst_offtarget_mismatch_pos",
    "offtarget_count",
    "worst_offtarget_chr",
    "worst_offtarget_pos",
    "worst_offtarget_strand",
    "worst_offtarget_sequence",
]


def aggregate_beguidedoff_results(ot_df):
    """Aggregate BEguidedOFF long-format output into per-guide summary.

    BEguidedOFF returns a DataFrame with columns including:
        guide_rna_sequence, ot_rank, ot_chromosome, ot_position,
        ot_strand, ot_sequence, ot_pam, n_mismatches,
        mismatch_positions, score, n_offtargets_below_threshold

    This function returns one row per guide with the single worst
    (most concerning = lowest CCTOP score) off-target.

    Args:
        ot_df: BEguidedOFF result DataFrame.

    Returns:
        DataFrame with columns: guide_seq, worst_offtarget_score,
        offtarget_count, worst_offtarget_chr, worst_offtarget_pos,
        worst_offtarget_sequence.
    """
    if ot_df.empty:
        return pd.DataFrame(columns=OTB_COLUMNS)

    required = {"guide_rna_sequence", "score"}
    if not required.issubset(ot_df.columns):
        raise ValueError(
            f"BEguidedOFF DataFrame missing required columns. "
            f"Expected at least {required}, got {set(ot_df.columns)}"
        )

    def _pick_worst(group):
        """Pick the single worst (lowest score) off-target for a guide."""
        worst = group.loc[group["score"].idxmin()]
        mm_raw = worst.get("mismatch_positions", "")
        if pd.isna(mm_raw):
            mm_raw = ""
        return pd.Series({
            "worst_offtarget_score": float(worst["score"]),
            "worst_offtarget_mismatch_pos": str(mm_raw),
            "offtarget_count": int(worst["n_offtargets_below_threshold"]),
            "worst_offtarget_chr": str(worst.get("ot_chromosome", "")),
            "worst_offtarget_pos": int(worst.get("ot_position", 0)),
            "worst_offtarget_strand": str(worst.get("ot_strand", "")),
            "worst_offtarget_sequence": str(worst.get("ot_sequence", "")),
        })

    summary = (
        ot_df.groupby("guide_rna_sequence", sort=False)
        .apply(_pick_worst, include_groups=False)
        .reset_index()
        .rename(columns={"guide_rna_sequence": "guide_seq"})
    )

    summary["worst_offtarget_pos"] = summary["worst_offtarget_pos"].astype("Int64")
    return summary[["guide_seq", "worst_offtarget_score",
                    "worst_offtarget_mismatch_pos",
                    "offtarget_count",
                    "worst_offtarget_chr", "worst_offtarget_pos",
                    "worst_offtarget_strand",
                    "worst_offtarget_sequence"]]


def merge_with_input(input_df, summary_df):
    """Merge off-target summary with input DataFrame.

    Args:
        input_df: DataFrame with guide results.
        summary_df: DataFrame from aggregate_beguidedoff_results.

    Returns:
        Merged DataFrame.
    """
    if "guide_rna_sequence" not in input_df.columns:
        raise ValueError("Input DataFrame missing 'guide_rna_sequence' column")

    merged = input_df.merge(
        summary_df,
        left_on="guide_rna_sequence",
        right_on="guide_seq",
        how="left",
    )
    merged = merged.drop(columns=["guide_seq"])
    # Guides with no off-target rows keep ``worst_offtarget_score`` = NaN:
    # the step-5 removal rule treats NaN explicitly as passing (nothing
    # concerning was found), the TSV renders it as an empty cell, and the web
    # UI shows it as "passes filter". (Filling with +inf here would leak
    # ``Infinity`` into the results JSON, which browsers cannot parse.)
    merged["worst_offtarget_mismatch_pos"] = merged["worst_offtarget_mismatch_pos"].fillna("").astype(str)
    merged["worst_offtarget_strand"] = merged["worst_offtarget_strand"].fillna("").astype(str)
    merged["offtarget_count"] = merged["offtarget_count"].fillna(0).astype(int)
    return merged
