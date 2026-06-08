"""Genomic region reconstruction from context columns."""

from beguided.genome.sequence import reverse_complement


def reconstruct_genomic_region(row):
    """Reconstruct the genomic region around a guide from context columns.

    Uses ``genomic_context_upstream`` and ``genomic_context_downstream``
    written by guidefinder, together with ``guide_rna_sequence``,
    ``pam_sequence`` and strand information, to rebuild the genomic
    DNA string that spans the guide + PAM + minimal flanking.

    Args:
        row: pandas Series or dict-like with required columns:
            strand_used, guide_rna_sequence, pam_sequence,
            editing_position, genome_base, desired_base,
            genomic_context_upstream, genomic_context_downstream.

    Returns:
        tuple: (region_seq, guide_offset) where *guide_offset* is the
        index in *region_seq* where the guide starts.
    """
    strand = row["strand_used"]
    guide_rna = row["guide_rna_sequence"]
    editing_pos = int(row["editing_position"])
    genome_base = row["genome_base"]
    desired_base = row["desired_base"]
    pam = row["pam_sequence"]
    context_up = str(row.get("genomic_context_upstream", "") or "")
    context_down = str(row.get("genomic_context_downstream", "") or "")

    if strand == "+":
        # guide_rna IS the genome guide (with desired base inserted)
        genome_guide = list(guide_rna)
        genome_guide[editing_pos - 1] = genome_base
        genome_guide = "".join(genome_guide)
        region = context_up + genome_guide + pam + context_down
        guide_offset = len(context_up)
    else:
        # guide_rna is reverse-complement of genome guide.
        # Reverse it back, then undo the mutation.
        modified_guide = reverse_complement(guide_rna)
        genome_guide = list(modified_guide)
        idx = len(genome_guide) - editing_pos
        genome_guide[idx] = genome_base
        genome_guide = "".join(genome_guide)
        genome_pam = reverse_complement(pam)
        # Genomic order for antisense:
        #   context_down (upstream of PAM) → PAM → guide → context_up (downstream of guide)
        region = context_down + genome_pam + genome_guide + context_up
        guide_offset = len(context_down) + len(genome_pam)

    return region, guide_offset
