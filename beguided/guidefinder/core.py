"""Guide RNA generation core logic."""

import logging

from beguided.genome.coordinates import (
    calculate_editing_position,
    calculate_search_range,
)
from beguided.genome.sequence import apply_reverse_complement_if_needed, insert_mutation
from beguided.editing.context import EditContextTable, is_abe, lookup_from_genome
from beguided.editing.mutations import find_bystander_edits
from beguided.editing.search import find_pams_regex
from beguided.guidefinder.crisprater import get_score as get_crisprater_score
from beguided.guidefinder.sequence_props import gc_content, has_homopolymer
from beguided.models import GuideResult, CORE_COLUMN_ORDER

logger = logging.getLogger("beguided.guidefinder")


def process_target(
    genome_seq, chromosome, position, target_base, desired_base,
    genome_base, base_editor_config, extra_columns=None,
    edit_context_table=None,
):
    """Process a single target and return list of GuideResult objects.

    Args:
        extra_columns: Optional dict of additional input columns to preserve.
        edit_context_table: Optional :class:`EditContextTable`. When supplied
            and the profile is ABE (sense or antisense), each valid guide
            gets a ``trinuc_edit_score`` from the Marquart et al. 2021
            normalised ABE8e context table. CBE profiles and ``None`` leave
            the field at ``None``.

    Returns:
        list of GuideResult (one per valid PAM, or a single special result).
    """
    results = []
    extra = extra_columns if extra_columns else {}

    if (
        base_editor_config.target_base != target_base
        or base_editor_config.desired_base != desired_base
    ):
        return [
            GuideResult(
                chromosome=chromosome,
                position=position,
                target_base=target_base,
                desired_base=desired_base,
                guide_rna_sequence="**Base-edit doesn't match base editor profile**",
                genome_base=genome_base,
                pam_sequence="XXX",
                editing_position=0,
                strand_used=base_editor_config.strand,
                extra_columns=extra.copy(),
            )
        ]

    guide_length = base_editor_config.guide_length
    pam_length = base_editor_config.pam_length
    min_pos = base_editor_config.high_edit_start
    max_pos = base_editor_config.high_edit_end
    strand = base_editor_config.strand

    start_search, end_search = calculate_search_range(
        position, min_pos, max_pos, strand, guide_length, pam_length
    )

    if strand == "+":
        end_search = min(len(genome_seq) - pam_length, end_search)

    pams_found = find_pams_regex(
        genome_seq, start_search, end_search,
        base_editor_config.pam_regex, pam_length
    )

    if not pams_found:
        return [
            GuideResult(
                chromosome=chromosome,
                position=position,
                target_base=target_base,
                desired_base=desired_base,
                guide_rna_sequence="**No guide possible**",
                genome_base=genome_base,
                pam_sequence="XXX",
                editing_position=0,
                strand_used=strand,
                extra_columns=extra.copy(),
            )
        ]

    for pam_pos in pams_found:
        if strand == "-":
            guide_start = pam_pos + pam_length
        else:
            guide_start = pam_pos - guide_length

        guide_end = guide_start + guide_length

        if guide_start < 0 or guide_end > len(genome_seq):
            continue

        guide_seq = genome_seq[guide_start:guide_end]
        if len(guide_seq) < guide_length:
            continue

        editing_pos = calculate_editing_position(
            position, guide_start, guide_length, strand
        )

        if not (min_pos <= editing_pos <= max_pos):
            continue
        if editing_pos < 1 or editing_pos > guide_length:
            continue

        modified_guide = insert_mutation(
            guide_seq, editing_pos, target_base, strand, guide_length
        )

        final_guide, pam_seq = apply_reverse_complement_if_needed(
            modified_guide,
            genome_seq[pam_pos:pam_pos + pam_length],
            base_editor_config.reverse_complement,
        )

        # Calculate CRISPRater on-target score (requires exactly 20 nt)
        guide_score = None
        if base_editor_config.guide_length == 20:
            try:
                guide_score = get_crisprater_score(final_guide)
            except Exception:
                logger.warning(
                    "Failed to calculate CRISPRater score for guide %s",
                    final_guide,
                )
                guide_score = None
        else:
            logger.warning(
                "CRISPRater score skipped: guide_length=%d != 20 "
                "for base editor %s",
                base_editor_config.guide_length,
                base_editor_config.name,
            )

        # Compute GC content (percent 0-100) and homopolymer flag for the
        # 20-nt protospacer. Both are length-agnostic, unlike CRISPRater,
        # so they are always populated for real guides.
        guide_gc_content = gc_content(final_guide)
        guide_has_homopolymer = has_homopolymer(final_guide)

        # Trinucleotide edit-context score (Marquart et al. 2021, ABE8e).
        # Only meaningful for ABE profiles; CBE rows stay None. The lookup
        # uses the actual genomic context around the target, so it is the
        # same value for every ABE profile selected for this target — the
        # editor identity only matters for whether the field is populated.
        guide_trinuc_score = None
        if edit_context_table is not None and is_abe(target_base, desired_base):
            guide_trinuc_score = lookup_from_genome(
                edit_context_table, genome_seq, position, strand, target_base,
            )

        bystander_positions = find_bystander_edits(
            final_guide, target_base, editing_pos,
            base_editor_config.general_edit_start,
            base_editor_config.general_edit_end,
            strand,
        )

        # Extract genomic context for downstream steps
        if strand == "+":
            up_start = max(0, guide_start - 4)
            context_up = genome_seq[up_start:guide_start]
            down_start = guide_end + pam_length
            context_down = genome_seq[down_start:down_start + 3]
        else:
            down_start = max(0, guide_start - pam_length - 3)
            context_down = genome_seq[down_start:guide_start - pam_length]
            up_start = guide_end
            context_up = genome_seq[up_start:up_start + 4]

        results.append(
            GuideResult(
                chromosome=chromosome,
                position=position,
                target_base=target_base,
                desired_base=desired_base,
                guide_rna_sequence=final_guide,
                genome_base=genome_base,
                pam_sequence=pam_seq,
                editing_position=editing_pos,
                strand_used=strand,
                guide_score=guide_score,
                gc_content=guide_gc_content,
                homopolymer=guide_has_homopolymer,
                trinuc_edit_score=guide_trinuc_score,
                bystander_positions=bystander_positions,
                genomic_context_upstream=context_up,
                genomic_context_downstream=context_down,
                extra_columns=extra.copy(),
            )
        )

    if results:
        return results
    else:
        return [
            GuideResult(
                chromosome=chromosome,
                position=position,
                target_base=target_base,
                desired_base=desired_base,
                guide_rna_sequence="**No guide possible**",
                genome_base=genome_base,
                pam_sequence="XXX",
                editing_position=0,
                strand_used=strand,
                extra_columns=extra.copy(),
            )
        ]


# Columns expected as core inputs to guidefinder
_CORE_INPUT_COLS = {
    "chromosome", "position", "target_base", "desired_base"
}


def run_genome_wide(
    genome, targets_df, profiles, base_editor_names,
    progress_callback=None, passthrough_extras=True,
    edit_context_table=None,
):
    """Run guidefinder across all targets and base editor profiles.

    Args:
        genome: dict of chromosome -> sequence.
        targets_df: pandas DataFrame with columns:
            chromosome, position, target_base, desired_base.
            Any additional columns are passed through to output unless
            ``passthrough_extras`` is ``False``.
        profiles: Profiles dict from load_base_editor_profiles.
        base_editor_names: List of base editor names to use.
        progress_callback: Optional callable(chrom, pos) for progress reporting.
        passthrough_extras: If ``True`` (default), any input columns beyond the
            four core columns (``chromosome``, ``position``, ``target_base``,
            ``desired_base``) are carried through to each ``GuideResult`` and
            later expanded into DataFrame columns. If ``False``, extra input
            columns are dropped.
        edit_context_table: Optional :class:`EditContextTable` forwarded to
            :func:`process_target` so each valid ABE guide gets a
            ``trinuc_edit_score`` from the Marquart et al. 2021 normalised
            ABE8e context table. ``None`` (default) leaves the field at
            ``None`` for every guide.

    Returns:
        list of GuideResult objects.
    """
    results = []

    # Identify extra columns to preserve from input
    if passthrough_extras:
        extra_cols = [c for c in targets_df.columns if c not in _CORE_INPUT_COLS]
    else:
        extra_cols = []

    modified_count = 0
    total_targets = 0

    for _, row in targets_df.iterrows():
        chromosome = str(row["chromosome"])
        position = int(row["position"])
        target_base = row["target_base"]
        desired_base = row["desired_base"]
        extra_columns = ({c: row[c] for c in extra_cols} if extra_cols else {})

        if progress_callback:
            progress_callback(chromosome, position)

        if chromosome not in genome:
            continue

        total_targets += 1
        genome_seq = genome[chromosome]
        genome_base = genome_seq[position - 1]

        if genome_base != target_base and genome_base != desired_base:
            logger.warning(
                "The genome base at location %s:%d "
                "neither matches the target_base nor the desired_base",
                chromosome, position
            )
        if genome_base != target_base:
            modified_count += 1

        editor_key = (target_base, desired_base)

        for base_editor_name in base_editor_names:
            if editor_key in profiles and base_editor_name in profiles[editor_key]:
                base_editor_config = profiles[editor_key][base_editor_name]
                target_results = process_target(
                    genome_seq, chromosome, position,
                    target_base, desired_base,
                    genome_base, base_editor_config, extra_columns=extra_columns,
                    edit_context_table=edit_context_table,
                )
                for result in target_results:
                    result.base_editor = base_editor_name
                results.extend(target_results)
            else:
                results.append(
                    GuideResult(
                        chromosome=chromosome,
                        position=position,
                        target_base=target_base,
                        desired_base=desired_base,
                        guide_rna_sequence="**No guide possible**",
                        genome_base=genome_base,
                        pam_sequence="XXX",
                        editing_position=0,
                        strand_used="",
                        base_editor=base_editor_name,
                        extra_columns=extra_columns.copy(),
                    )
                )

    if modified_count > 0:
        logger.warning(
            "From %d targets the genome_base was modified to match the target_base",
            modified_count
        )

    return results


def results_to_dataframe(results, passthrough_extras=True):
    """Convert list of GuideResult to a pandas DataFrame with proper column ordering.

    Any ``extra_columns`` stored on each GuideResult are expanded into columns
    so that original input metadata is preserved all the way to the final output
    (unless ``passthrough_extras`` is ``False``).

    Args:
        results: list of GuideResult.
        passthrough_extras: If ``True`` (default), the ``extra_columns`` dict on
            each ``GuideResult`` is expanded into DataFrame columns. If
            ``False``, extra input columns are dropped here as well.

    Returns:
        pandas DataFrame.
    """
    import pandas as pd

    if not results:
        return pd.DataFrame(columns=CORE_COLUMN_ORDER)

    # Build rows while expanding extra_columns
    rows = []
    extra_keys = set()
    for r in results:
        row = vars(r).copy()
        if passthrough_extras:
            extra = row.pop("extra_columns", {})
            extra_keys.update(extra.keys())
            row.update(extra)
        else:
            row.pop("extra_columns", None)
        rows.append(row)

    df = pd.DataFrame(rows)
    extra_keys_sorted = sorted(extra_keys)

    # Column order: CORE_COLUMN_ORDER columns first, then any remaining
    # internal columns, then user-provided extra_columns at the very end.
    internal_cols = [
        "genomic_context_upstream",
        "genomic_context_downstream",
    ]
    seen = set()
    col_order = []
    for c in CORE_COLUMN_ORDER:
        if c in df.columns and c not in seen:
            col_order.append(c)
            seen.add(c)
    for c in internal_cols:
        if c in df.columns and c not in seen:
            col_order.append(c)
            seen.add(c)
    for c in extra_keys_sorted:
        if c in df.columns and c not in seen:
            col_order.append(c)
            seen.add(c)
    df = df[col_order]
    df = df.sort_values(
        by=["chromosome", "position", "base_editor", "editing_position"],
        ascending=[True, True, True, True]
    )
    return df
