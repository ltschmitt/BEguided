"""Bystander effect analysis for a single guide."""

from Bio.Seq import Seq

from beguided.genome.coordinates import guide_coordinates
from beguided.genome.context import reconstruct_genomic_region
from beguided.genome.sequence import COMPLEMENT
from beguided.editing.mutations import (
    find_bystander_edits,
    mark_mutations,
)
from beguided.bystander.annotation import (
    cumulative_mrna_pos,
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

    (
        coding_seq, codon_offset_start, codon_offset_end, cds_strand,
        extract_start, extract_end,
    ) = extract_coding_sequence_from_region(
        chrom, guide_start, guide_end, cds_info, region, guide_offset,
        genome=genome, interval_trees=interval_trees,
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

    # Check editing position in CDS FIRST (before bystander check), so
    # that gene_edit_desired can be computed whenever the target itself is
    # in the CDS — even if a bystander later turns out to be outside.
    if strand_used == "+":
        genomic_editing_pos = guide_start + (editing_position - 1)
    else:
        genomic_editing_pos = guide_end - editing_position
    target_in_cds = cds_start <= genomic_editing_pos < cds_end

    # ``extract_start`` / ``extract_end`` are the actual genomic boundaries
    # of the codon-aligned window (returned by extract_coding_sequence_from_region
    # via _codon_aligned_bounds). The OLD formula's ``codon_offset_start`` /
    # ``codon_offset_end`` are genomic distances only on the local-region path;
    # on the genome-fallback path they are splice-aware indices into
    # ``coding_seq``. Using the real ``extract_start`` / ``extract_end``
    # ensures the splicing-aware mep lookup below uses the correct anchor.

    # mRNA-side base for the target / desired edit at the editing position.
    mrna_target_base_mrna = (
        COMPLEMENT[target_base] if cds_strand == "-" else target_base
    )
    mrna_desired_base_mrna = (
        COMPLEMENT[desired_base] if cds_strand == "-" else desired_base
    )

    # mRNA-side local sequence (already reverse-complemented for - strand CDS).
    mrna_seq_unmutated = coding_seq
    if cds_strand == "-":
        mrna_seq_unmutated = str(Seq(coding_seq).reverse_complement())

    # Compute the 0-based index of the editing base in ``mrna_seq``.
    #
    # The previous implementation used the genomic distances
    # ``codon_offset_start`` / ``codon_offset_end`` as indices into the
    # *spliced* ``mrna_seq``. Those two values agree only when the
    # codon-aligned window does not span an intron; when it does (e.g. a
    # 44 bp window that splices down to 24 nt) the genomic distance is
    # much larger than the corresponding ``mrna_seq`` index, so the
    # computed edit position lands out of bounds and ``gene_edit_desired``
    # ends up ``"?"`` with no brackets in the AA / simulated-outcome
    # columns.
    #
    # For the common case where the OLD formula gives an in-bounds ``mep``
    # we use it directly (it works correctly for non-splicing windows and
    # for cross-strand ``-`` CDS strand guides whose codon-aligned window
    # happens to splice but where the OLD local extraction still falls
    # inside the right codon). For windows where the OLD ``mep`` is
    # out-of-bounds we recompute with ``cumulative_mrna_pos``, which is
    # splicing-aware and gives the correct index into ``mrna_seq``.
    #
    # ``codon_offset_end`` is used for ``-`` strand guides and
    # ``codon_offset_start`` for ``+`` strand guides, matching the original
    # gene_edit_desired computation.
    if strand_used == "-":
        codon_offset_for_old = codon_offset_end
    else:
        codon_offset_for_old = codon_offset_start

    if strand_used != cds_strand:
        old_mep = len(coding_seq) - codon_offset_for_old - editing_position
    else:
        old_mep = codon_offset_for_old + editing_position - 1

    # Use the OLD ``mep`` if it falls inside the in-bounds slice for the
    # codon (accounting for Python's negative-index wrapping, which the
    # OLD formula's ``mrna_seq[codon_start:codon_start + 3]`` relies on).
    old_local_codon = (old_mep // 3) + 1 if old_mep >= -len(coding_seq) else None
    if old_local_codon is not None:
        old_codon_start = (old_local_codon - 1) * 3
        old_mrna_codon = mrna_seq_unmutated[old_codon_start:old_codon_start + 3]
        use_old_mep = len(old_mrna_codon) == 3
    else:
        use_old_mep = False

    def _mep_for_genomic_pos(genomic_pos):
        mrna_pos = cumulative_mrna_pos(
            transcript_id, genomic_pos, interval_trees
        )
        if mrna_pos is None:
            return None
        anchor_mrna_pos = (
            cumulative_mrna_pos(transcript_id, extract_start, interval_trees)
            if cds_strand == "+"
            else cumulative_mrna_pos(transcript_id, extract_end - 1, interval_trees)
        )
        if anchor_mrna_pos is None:
            return None
        return mrna_pos - anchor_mrna_pos

    if use_old_mep:
        target_mep = old_mep
    else:
        target_mep = _mep_for_genomic_pos(genomic_editing_pos)

    bystander_meps = []
    for bp in bystander_list:
        bp_genomic = (
            guide_start + (bp - 1) if strand_used == "+" else guide_end - bp
        )
        # Same OLD-style formula with ``codon_offset_for_old`` so the
        # bystander positions stay aligned with the target bracket when
        # the OLD ``mep`` is in bounds.
        if strand_used != cds_strand:
            bp_old = len(coding_seq) - codon_offset_for_old - bp
        else:
            bp_old = codon_offset_for_old + bp - 1
        # Only accept the OLD bystander position if its slice gives a
        # valid 1-char substring (accounting for negative-index wrapping).
        if -len(coding_seq) <= bp_old < len(coding_seq):
            bystander_meps.append(bp_old)
        else:
            bystander_meps.append(_mep_for_genomic_pos(bp_genomic))

    # Compute gene_edit_desired: "<gene> <Target_AA><abs_codon><Desired_AA>"
    # Target_AA: translate the codon with target_base at the editing position.
    # Desired_AA: translate the codon with desired_base at the editing position.
    # Only the target edit is considered (bystanders are ignored). The other
    # two bases of the codon come from the local mRNA context.
    gene_edit_desired = ""
    local_codon = None
    if target_in_cds and gene_name:
        abs_codon = absolute_codon_number(
            transcript_id, genomic_editing_pos, interval_trees
        )
        if abs_codon is not None:
            if target_mep is not None:
                local_codon = (target_mep // 3) + 1
                offset_in_codon = target_mep % 3
                codon_start = (local_codon - 1) * 3
                mrna_codon = mrna_seq_unmutated[codon_start:codon_start + 3]
                if len(mrna_codon) == 3:
                    target_codon = (
                        mrna_codon[:offset_in_codon]
                        + mrna_target_base_mrna
                        + mrna_codon[offset_in_codon + 1:]
                    )
                    desired_codon = (
                        mrna_codon[:offset_in_codon]
                        + mrna_desired_base_mrna
                        + mrna_codon[offset_in_codon + 1:]
                    )
                    target_aa = genetic_code.get(target_codon, "?")
                    desired_aa = genetic_code.get(desired_codon, "?")
                else:
                    # ``mrna_codon`` shorter than 3 (e.g. the codon-aligned
                    # window ran off the end of the CDS). Fall back to the
                    # ``?`` placeholder so downstream tooling still sees the
                    # codon number rather than a silently-empty string.
                    target_aa = "?"
                    desired_aa = "?"
            else:
                # Splicing-aware ``mep`` could not be resolved (e.g. the
                # edit lands in an intron). Keep the abs_codon but mark the
                # AA letters as unknown; downstream tools still see a useful
                # ``<gene> ?<abs_codon>?`` string rather than an empty cell.
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

    # Build mutated mRNA sequences by inserting ``mrna_desired_base_mrna``
    # directly at the splicing-aware mep positions. This replaces the
    # previous ``simulate_editing`` call which was tied to the same
    # buggy genomic-distance ``codon_offset`` and silently dropped
    # mutations that landed out of bounds.
    def _mutate(positions):
        s = list(mrna_seq_unmutated)
        for p in positions:
            if p is not None and 0 <= p < len(s):
                s[p] = mrna_desired_base_mrna
        return "".join(s)

    seq_wanted = _mutate([target_mep])
    seq_bystander = _mutate([target_mep] + bystander_meps)

    _, aa_wanted = translate_with_codons(seq_wanted)
    _, aa_bystander = translate_with_codons(seq_bystander)

    if not effect:
        if not has_bystanders:
            effect = "no_bystanders"
        elif aa_wanted == aa_bystander:
            effect = "silent"
        else:
            effect = "missense"

    # Bracket the editing base (target [] and bystanders ()) in
    # simulated_outcome directly at the splice-aware mep positions. We
    # do this directly instead of going through ``mark_mutations`` because
    # ``mark_mutations`` reverses the sequence for the cross-strand case
    # (strand_used != cds_strand), which is unnecessary when the sequence
    # is already in mRNA order.
    def _bracket(seq, pos, kind="target"):
        if pos is None or not (0 <= pos < len(seq)):
            return seq
        open_b, close_b = ("[", "]") if kind == "target" else ("(", ")")
        return (
            seq[:pos]
            + open_b + seq[pos] + close_b
            + seq[pos + 1:]
        )

    if target_mep is not None:
        seq_bystander = _bracket(seq_bystander, target_mep, kind="target")
    for bp_mep in bystander_meps:
        seq_bystander = _bracket(seq_bystander, bp_mep, kind="bystander")

    # Bracket the editing codon in AA sequences. ``local_codon`` may be
    # None if the codon identification fell through; in that case we leave
    # the AA strings un-bracketed (preserves prior fallback behavior).
    if local_codon is not None:
        aa_wanted = mark_mutations(aa_wanted, local_codon, [], strand_used, cds_strand, 0)
        aa_bystander = mark_mutations(aa_bystander, local_codon, [], strand_used, cds_strand, 0)

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
