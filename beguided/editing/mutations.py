"""Mutation simulation and bystander detection."""

from beguided.genome.sequence import COMPLEMENT


def find_bystander_edits(
    guide_seq, target_base, editing_position,
    general_edit_start, general_edit_end, strand
):
    """Find bystander edit positions within the general editing window.

    Args:
        guide_seq: Guide RNA sequence string.
        target_base: Base that gets edited by the base editor.
        editing_position: 1-based target editing position (to exclude).
        general_edit_start: Start of general edit window (1-based).
        general_edit_end: End of general edit window.
        strand: '+' or '-'.

    Returns:
        Comma-separated string of bystander positions, or empty string.
    """
    if strand == "-":
        guide_base = COMPLEMENT[target_base]
    else:
        guide_base = target_base

    bystanders = []
    for pos in range(general_edit_start, general_edit_end + 1):
        if pos < 1 or pos > len(guide_seq):
            continue
        if pos == editing_position:
            continue
        if guide_seq[pos - 1] == guide_base:
            bystanders.append(str(pos))
    return ",".join(bystanders)


def simulate_editing(
    coding_seq, codon_offset, edit_positions, desired_base, strand_used
):
    """Simulate editing outcome at specified positions.

    Args:
        coding_seq: Coding sequence string (may include flanking).
        codon_offset: Offset in coding_seq where guide starts.
        edit_positions: List of positions to edit (1-based from 5' end).
        desired_base: Base to change to.
        strand_used: '+' or '-'.

    Returns:
        Edited sequence string.
    """
    seq_list = list(coding_seq)
    if strand_used == "-":
        seq_list = seq_list[::-1]

    for pos in edit_positions:
        pos_in_coding = codon_offset + (pos - 1)
        if 0 <= pos_in_coding < len(seq_list):
            seq_list[pos_in_coding] = desired_base

    if strand_used == "-":
        seq_list = seq_list[::-1]

    return "".join(seq_list)


def mark_mutations(
    seq, editing_position, bystander_list, strand_used, cds_strand, codon_offset
):
    """Add brackets around mutated positions.

    [] for target edit, () for bystanders.

    Args:
        seq: Sequence string.
        editing_position: 1-based target position.
        bystander_list: List of 1-based bystander positions.
        strand_used: Strand of guide ('+' or '-').
        cds_strand: Strand of CDS ('+' or '-').
        codon_offset: Offset where guide starts in seq.

    Returns:
        Bracketed sequence string.
    """
    seq_list = list(seq)

    if strand_used != cds_strand:
        seq_list = seq_list[::-1]

    editing_pos_adj = codon_offset + editing_position - 1
    if 0 <= editing_pos_adj < len(seq_list):
        seq_list[editing_pos_adj] = "[" + seq_list[editing_pos_adj] + "]"

    for pos in bystander_list:
        pos_in_coding = codon_offset + (pos - 1)
        if 0 <= pos_in_coding < len(seq_list):
            seq_list[pos_in_coding] = "(" + seq_list[pos_in_coding] + ")"

    if strand_used != cds_strand:
        seq_list = seq_list[::-1]

    return "".join(seq_list)
