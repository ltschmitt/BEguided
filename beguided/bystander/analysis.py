"""Bystander effect analysis for a single guide."""

from Bio.Seq import Seq

from beguided.genome.coordinates import guide_coordinates
from beguided.genome.context import reconstruct_genomic_region
from beguided.genome.sequence import COMPLEMENT
from beguided.editing.mutations import (
    find_bystander_edits,
    mark_mutations,
    simulate_editing,
)
from beguided.bystander.annotation import (
    get_cds_overlap,
    extract_coding_sequence_from_region,
    absolute_codon_number,
)
from beguided.bystander.protein import translate_with_codons, genetic_code
from beguided.models import CORE_COLUMN_ORDER


def _normalize_bases(base_from, base_to):
    """Return (forward_from, forward_to) based on desired_base."""
    return str(base_from), str(base_to)


def analyze_guide(row, interval_trees, genome=None):
    """Analyze a single guide for bystander effects.

    Args:
        row: pandas Series or dict-like with required columns.
        interval_trees: dict from load_gtf_to_interval_tree.
        genome: Optional ``{chrom: sequence_string}`` dict used as a
            fallback when the codon-aligned window crosses an exon
            boundary. When ``None`` the analysis uses the local region
            only; in that case guides whose codon spans an exon
            junction will be reported as ``cds_reconstruction_failed``.

    Returns:
        tuple: (bystander_effect, simulated_outcome, aa_wanted, aa_bystander,
                gene_edit_desired)
            bystander_effect is one of the labels documented in the module
            (e.g. ``no_bystanders``, ``silent``, ``missense``, ``multi-transcript``,
            ``bystander_outside_cds``, or one of the granular ``*_non_coding``
            / ``*_failed`` labels when CDS handling fails).
            gene_edit_desired is ``"<gene_name> <Target_AA><codon_number><Desired_AA>"``
            (e.g. ``"TP53 *306R"``) when the target lies in a CDS with a resolvable
            gene name and translatable codon, otherwise ``""``. Target_AA is the
            amino acid when ``target_base`` is placed at the editing position;
            Desired_AA is the amino acid when ``desired_base`` is placed there.
            The codon number is the absolute position within the full protein
            (not the local CDS window). Only the target edit is considered
            (bystanders are ignored).
    """
    chrom = str(row["chromosome"])
    position = int(row["position"])
    strand_used = row["strand_used"]
    editing_position = int(row["editing_position"])
    bystander_positions = row.get("bystander_positions", "")
    target_base = row["target_base"]
    desired_base = row["desired_base"]

    guide_start, guide_end = guide_coordinates(
        position, editing_position, strand_used, guide_length=20
    )

    cds_info, is_multiple = get_cds_overlap(
        chrom, guide_start, guide_end, interval_trees
    )

    if not cds_info:
        return "no_cds_overlap", "", "", "", ""

    if is_multiple:
        effect = "multi-transcript"
    else:
        effect = ""

    gene_name = cds_info[0][4] if len(cds_info[0]) > 4 else ""
    transcript_id = cds_info[0][3] if len(cds_info[0]) > 3 else ""

    region, guide_offset = reconstruct_genomic_region(row)

    coding_seq, codon_offset_start, codon_offset_end, cds_strand = (
        extract_coding_sequence_from_region(
            chrom, guide_start, guide_end, cds_info, region, guide_offset,
            genome=genome, interval_trees=interval_trees,
        )
    )

    if coding_seq is None:
        return "cds_reconstruction_failed", "", "", "", ""

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

    cds_start, cds_end = cds_info[0][0], cds_info[0][1]
    codon_offset = codon_offset_end if strand_used == "-" else codon_offset_start

    # Check editing position in CDS FIRST (before bystander check), so
    # that gene_edit_desired can be computed whenever the target itself is
    # in the CDS — even if a bystander later turns out to be outside.
    if strand_used == "+":
        genomic_editing_pos = guide_start + (editing_position - 1)
    else:
        genomic_editing_pos = guide_end - editing_position
    target_in_cds = cds_start <= genomic_editing_pos < cds_end

    # Compute gene_edit_desired: "<gene> <Target_AA><abs_codon><Desired_AA>"
    # Target_AA: translate the codon with target_base at the editing position.
    # Desired_AA: translate the codon with desired_base at the editing position.
    # Only the target edit is considered (bystanders are ignored). The other
    # two bases of the codon come from the genome context.
    gene_edit_desired = ""
    # The codon-aligned extract starts at a codon boundary in the full
    # transcript frame, so the position of the editing base within the
    # local extract is just ``codon_offset + editing_position - 1`` for
    # sense guides (or its reverse-complement mirror for antisense).
    if strand_used != cds_strand:
        mrna_edit_pos = len(coding_seq) - codon_offset - editing_position
    else:
        mrna_edit_pos = codon_offset + editing_position - 1
    local_codon = (mrna_edit_pos // 3) + 1
    if target_in_cds and gene_name:
        abs_codon = absolute_codon_number(
            transcript_id, genomic_editing_pos, interval_trees
        )
        if abs_codon is not None:
            mrna_seq = coding_seq
            if cds_strand == "-":
                mrna_seq = str(Seq(coding_seq).reverse_complement())
            codon_start = (local_codon - 1) * 3
            mrna_codon = mrna_seq[codon_start:codon_start + 3]
            offset_in_codon = mrna_edit_pos % 3
            mrna_target_base = (
                COMPLEMENT[target_base]
                if cds_strand == "-"
                else target_base
            )
            mrna_desired_base = (
                COMPLEMENT[desired_base]
                if cds_strand == "-"
                else desired_base
            )
            if len(mrna_codon) == 3:
                target_codon = (
                    mrna_codon[:offset_in_codon]
                    + mrna_target_base
                    + mrna_codon[offset_in_codon + 1:]
                )
                desired_codon = (
                    mrna_codon[:offset_in_codon]
                    + mrna_desired_base
                    + mrna_codon[offset_in_codon + 1:]
                )
                target_aa = genetic_code.get(target_codon, "?")
                desired_aa = genetic_code.get(desired_codon, "?")
            else:
                target_aa = "?"
                desired_aa = "?"
            gene_edit_desired = (
                f"{gene_name} {target_aa}{abs_codon}{desired_aa}"
            )

    if not target_in_cds:
        return "editing_position_non_coding", "", "", "", ""

    # Check for bystanders outside CDS (target is already confirmed in CDS).
    # A bystander is "outside CDS" if its genomic position falls outside
    # the CDS exon that the guide overlaps (i.e. it would live in a
    # different exon of the same transcript). Such bystanders are
    # reported separately because translating them in the local extract
    # frame would mix codon boundaries across exons.
    for pos in bystander_list:
        if strand_used == "+":
            genomic_pos = guide_start + (pos - 1)
        else:
            genomic_pos = guide_end - pos
        if genomic_pos < cds_start or genomic_pos >= cds_end:
            return "bystander_outside_cds", "", "", "", gene_edit_desired

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

    aa_editing_pos = local_codon
    aa_wanted = mark_mutations(aa_wanted, aa_editing_pos, [], strand_used, cds_strand, 0)
    aa_bystander = mark_mutations(aa_bystander, aa_editing_pos, [], strand_used, cds_strand, 0)

    return effect, seq_bystander, aa_wanted, aa_bystander, gene_edit_desired


def run_all(guides_df, interval_trees, progress_callback=None, genome=None):
    """Run bystander analysis on all guides in a DataFrame.

    Args:
        guides_df: DataFrame with guide results.
        interval_trees: dict from load_gtf_to_interval_tree.
        progress_callback: Optional callable(idx, total) for progress.
        genome: Optional ``{chrom: sequence_string}`` dict. Forwarded to
            :func:`analyze_guide` so the codon-aligned window can splice
            across exon junctions when needed.

    Returns:
        DataFrame with added bystander columns.
    """
    import pandas as pd

    results = []
    total = len(guides_df)

    for idx, (_, row) in enumerate(guides_df.iterrows()):
        if progress_callback and idx % 100 == 0:
            progress_callback(idx, total)

        effect, simulated_outcome, aa_wanted, aa_bystander, gene_edit_desired = (
            analyze_guide(row, interval_trees, genome=genome)
        )
        results.append({
            "bystander_effect": effect,
            "simulated_outcome": simulated_outcome,
            "aa_wanted": aa_wanted,
            "aa_bystander": aa_bystander,
            "gene_edit_desired": gene_edit_desired,
        })

    results_df = pd.DataFrame(results)
    df = pd.concat([guides_df.reset_index(drop=True), results_df], axis=1)

    new_cols = [
        "bystander_effect",
        "gene_edit_desired",
        "simulated_outcome",
        "aa_wanted",
        "aa_bystander",
    ]
    # Column order: CORE_COLUMN_ORDER columns first (in CORE_COLUMN_ORDER), then
    # bystander-specific columns, then any extra columns from the input TSV.
    seen = set()
    cols = []
    for c in CORE_COLUMN_ORDER:
        if c in df.columns and c not in seen:
            cols.append(c)
            seen.add(c)
    for c in new_cols:
        if c in df.columns and c not in seen:
            cols.append(c)
            seen.add(c)
    for c in df.columns:
        if c not in seen:
            cols.append(c)
            seen.add(c)
    df = df[cols]
    return df
