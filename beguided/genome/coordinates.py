"""Genomic coordinate calculations."""


def calculate_search_range(
    position, min_pos, max_pos, strand, guide_length, pam_length
):
    """Calculate PAM search range for a given editing window.

    For antisense strand: PAM is upstream of target.
    For sense strand: PAM is downstream of target.

    Args:
        position: 1-based genomic position of target.
        min_pos: High-edit window start (1-based from 5' end).
        max_pos: High-edit window end.
        strand: '+' or '-'.
        guide_length: Length of guide RNA.
        pam_length: Length of PAM sequence.

    Returns:
        tuple: (start_search, end_search) in 0-based coordinates.
    """
    if strand == "-":
        base = position - 1
        offset = guide_length + pam_length
        start = base - (offset - min_pos)
        end = base - (offset - max_pos) + pam_length
    else:
        start = position + (guide_length - max_pos)
        end = position + (guide_length - min_pos) + pam_length
    return max(0, start), max(0, end)


def calculate_editing_position(position, guide_start, guide_length, strand):
    """Calculate editing position (1-based from 5' end of guide).

    Args:
        position: 1-based genomic target position.
        guide_start: 0-based start of guide in genome.
        guide_length: Length of guide RNA.
        strand: '+' or '-'.

    Returns:
        int: 1-based editing position from 5' end.
    """
    if strand == "-":
        return guide_length - ((position - 1) - guide_start)
    else:
        return (position - 1) - guide_start + 1


def guide_coordinates(position, editing_position, strand, guide_length=20):
    """Calculate guide start/end genomic coordinates.

    Args:
        position: 1-based genomic target position.
        editing_position: 1-based from 5' end.
        strand: '+' or '-'.
        guide_length: Length of guide (default 20).

    Returns:
        tuple: (guide_start, guide_end) in 0-based coordinates.
    """
    if strand == "+":
        guide_start = position - editing_position
        guide_end = guide_start + guide_length
    else:
        guide_end = position + (editing_position - 1)
        guide_start = guide_end - guide_length
    return guide_start, guide_end
