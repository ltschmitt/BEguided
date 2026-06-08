"""Bystander effect analysis for a single guide."""

from Bio.Seq import Seq

from beguided.genome.coordinates import guide_coordinates
from beguided.genome.context import reconstruct_genomic_region
from beguided.editing.mutations import (
    find_bystander_edits,
    mark_mutations,
    simulate_editing,
)
from beguided.bystander.annotation import (
    get_cds_overlap,
    extract_coding_sequence_from_region,
)
from beguided.bystander.protein import translate_with_codons


def _normalize_bases(base_from, base_to):
    """Return (forward_from, forward_to) based on desired_base."""
    return str(base_from), str(base_to)


def analyze_guide(row, interval_trees):
    """Analyze a single guide for bystander effects.

    Args:
        row: pandas Series or dict-like with required columns.
        interval_trees: dict from load_gtf_to_interval_tree.

    Returns:
        tuple: (effect, simulated_outcome, aa_wanted, aa_bystander, effect_detail)
    """
    chrom = str(row["chromosome"])
    position = int(row["position"])
    strand_used = row["strand_used"]
    editing_position = int(row["editing_position"])
    bystander_positions = row.get("bystander_positions", "")
    desired_base = row["desired_base"]

    guide_start, guide_end = guide_coordinates(
        position, editing_position, strand_used, guide_length=20
    )

    cds_info, is_multiple = get_cds_overlap(
        chrom, guide_start, guide_end, interval_trees
    )

    if not cds_info:
        return "non-coding", "", "", "", "no_cds_overlap"

    effect = "multi-transcript" if is_multiple else ""

    region, guide_offset = reconstruct_genomic_region(row)

    coding_seq, codon_offset_start, codon_offset_end, cds_strand = (
        extract_coding_sequence_from_region(
            chrom, guide_start, guide_end, cds_info, region, guide_offset
        )
    )

    if coding_seq is None:
        return "non-coding", "", "", "", "cds_reconstruction_failed"

    # Parse bystanders
    bystander_list = []
    if (
        bystander_positions
        and isinstance(bystander_positions, str)
        and bystander_positions.strip()
    ):
        for pos_str in bystander_positions.split(","):
            ps = pos_str.strip()
            if ps:
                bystander_list.append(int(ps))

    has_bystanders = len(bystander_list) > 0

    # Check for bystanders outside CDS
    cds_start, cds_end = cds_info[0][0], cds_info[0][1]
    for pos in bystander_list:
        if strand_used == "+":
            genomic_pos = guide_start + (pos - 1)
        else:
            genomic_pos = guide_end - pos
        if genomic_pos < cds_start or genomic_pos >= cds_end:
            return "bystander_outside_cds", "", "", "", "bystander_outside_cds"

    # Check editing position in CDS
    if strand_used == "+":
        genomic_editing_pos = guide_start + (editing_position - 1)
    else:
        genomic_editing_pos = guide_end - editing_position
    if genomic_editing_pos < cds_start or genomic_editing_pos >= cds_end:
        return "non-coding", "", "", "", "editing_position_non_coding"

    codon_offset = codon_offset_end if strand_used == "-" else codon_offset_start

    seq_wanted = simulate_editing(
        coding_seq, codon_offset, [editing_position], desired_base, strand_used
    )
    all_edit_positions = [editing_position] + bystander_list
    seq_bystander = simulate_editing(
        coding_seq, codon_offset, all_edit_positions, desired_base, strand_used
    )

    if cds_strand == "-":
        seq_wanted = str(Seq(seq_wanted).reverse_complement())
        seq_bystander = str(Seq(seq_bystander).reverse_complement())

    _, aa_wanted = translate_with_codons(seq_wanted)
    _, aa_bystander = translate_with_codons(seq_bystander)

    if not effect:
        if not has_bystanders:
            effect = "no_bystanders"
        elif aa_wanted == aa_bystander:
            effect = "silent"
        else:
            effect = "missense"

    seq_bystander = mark_mutations(
        seq_bystander, editing_position, bystander_list,
        strand_used, cds_strand, codon_offset
    )

    aa_editing_pos = ((editing_position + codon_offset - 1) // 3) + 1
    aa_wanted = mark_mutations(aa_wanted, aa_editing_pos, [], strand_used, cds_strand, 0)
    aa_bystander = mark_mutations(aa_bystander, aa_editing_pos, [], strand_used, cds_strand, 0)

    effect_detail = effect
    return effect, seq_bystander, aa_wanted, aa_bystander, effect_detail


def run_all(guides_df, interval_trees, progress_callback=None):
    """Run bystander analysis on all guides in a DataFrame.

    Args:
        guides_df: DataFrame with guide results.
        interval_trees: dict from load_gtf_to_interval_tree.
        progress_callback: Optional callable(idx, total) for progress.

    Returns:
        DataFrame with added bystander columns.
    """
    import pandas as pd

    results = []
    total = len(guides_df)

    for idx, (_, row) in enumerate(guides_df.iterrows()):
        if progress_callback and idx % 100 == 0:
            progress_callback(idx, total)

        effect, simulated_outcome, aa_wanted, aa_bystander, effect_detail = analyze_guide(
            row, interval_trees
        )
        results.append({
            "bystander_effect": effect,
            "simulated_outcome": simulated_outcome,
            "aa_wanted": aa_wanted,
            "aa_bystander": aa_bystander,
            "bystander_effect_detail": effect_detail,
        })

    results_df = pd.DataFrame(results)
    df = pd.concat([guides_df.reset_index(drop=True), results_df], axis=1)

    cols = [c for c in df.columns if c not in [
        "bystander_effect", "simulated_outcome", "aa_wanted", "aa_bystander",
        "bystander_effect_detail"
    ]]
    cols.extend(["bystander_effect", "bystander_effect_detail",
                 "simulated_outcome", "aa_wanted", "aa_bystander"])
    df = df[cols]
    return df
