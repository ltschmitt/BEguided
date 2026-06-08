"""Genome annotation operations for bystander analysis."""

import re

from intervaltree import IntervalTree


def load_gtf_to_interval_tree(gtf_file):
    """Load GTF file and build interval trees for MANE_Select CDS per chromosome.

    Returns:
        dict: {chromosome: IntervalTree([(start, end, (strand, transcript_id)), ...])}
    """
    trees = {}
    cds_count = 0
    skipped_count = 0

    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue

            attributes = parts[8]
            if "MANE_Select" not in attributes:
                skipped_count += 1
                continue

            chrom = parts[0]
            start = int(parts[3])  # 1-based, inclusive
            end = int(parts[4])    # 1-based, inclusive
            strand = parts[6]

            transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
            transcript_id = transcript_match.group(1) if transcript_match else "unknown"

            if chrom not in trees:
                trees[chrom] = IntervalTree()

            # intervaltree uses 0-based, half-open [start, end)
            trees[chrom][start - 1:end] = (strand, transcript_id)
            cds_count += 1

    return trees


def get_cds_overlap(chrom, guide_start, guide_end, interval_trees):
    """Find CDS overlaps for a guide region.

    Returns:
        tuple: (list of (cds_start, cds_end, strand, transcript_id), is_multiple)
    """
    if chrom not in interval_trees:
        return [], False

    overlaps = interval_trees[chrom][guide_start:guide_end]
    if not overlaps:
        return [], False

    cds_list = []
    for interval in overlaps:
        strand, transcript_id = interval.data
        cds_list.append((interval.begin, interval.end, strand, transcript_id))

    is_multiple = len(cds_list) > 1
    return cds_list, is_multiple


def extract_coding_sequence(chrom, guide_start, guide_end, cds_info, genome):
    """Extract the coding sequence covering the guide region with flanking.

    Returns:
        tuple: (coding_seq, codon_offset_start, codon_offset_end, cds_strand)
    """
    if not cds_info:
        return None, None, None, None

    cds_start, cds_end, cds_strand, transcript_id = cds_info[0]

    offset_from_cds_start = guide_start - cds_start
    needed_upstream = offset_from_cds_start % 3
    extract_start = max(guide_start - needed_upstream, cds_start)

    offset_from_cds_end = guide_end - extract_start
    needed_downstream = 3 - (offset_from_cds_end % 3)
    extract_end = min(guide_end + needed_downstream, cds_end)

    seq = genome[chrom][extract_start:extract_end].upper()
    codon_offset_start = guide_start - extract_start
    codon_offset_end = extract_end - guide_end

    return seq, codon_offset_start, codon_offset_end, cds_strand


def extract_coding_sequence_from_region(chrom, guide_start, guide_end,
                                          cds_info, region, guide_offset):
    """Extract coding sequence from a pre-reconstructed region string.

    The *region* string is built by ``reconstruct_genomic_region`` and covers
    the guide + PAM + a few flanking bases.  *guide_offset* tells where the
    guide (not the PAM) starts inside *region*.

    Because the region does **not** include the whole genome, we can only
    extract bases that fall within the CDS bounds AND are present in
    *region*.  If the needed codon-flanking would extend beyond the region
    boundaries, we return ``None`` (the context is too short to perform
    reliable translation).

    Args:
        chrom: Chromosome name.
        guide_start: 0-based genomic start of guide.
        guide_end: 0-based genomic end of guide.
        cds_info: List of (cds_start, cds_end, cds_strand, transcript_id).
        region: Reconstructed genomic region string.
        guide_offset: Index in *region* where the guide starts.

    Returns:
        Same tuple as ``extract_coding_sequence``; coding_seq is None if
        required flanking extends beyond the region.
    """
    if not cds_info:
        return None, None, None, None

    cds_start, cds_end, cds_strand, transcript_id = cds_info[0]

    offset_from_cds_start = guide_start - cds_start
    needed_upstream = offset_from_cds_start % 3
    extract_start = max(guide_start - needed_upstream, cds_start)

    offset_from_cds_end = guide_end - extract_start
    needed_downstream = 3 - (offset_from_cds_end % 3)
    extract_end = min(guide_end + needed_downstream, cds_end)

    region_guide_start = guide_offset
    region_guide_end = guide_offset + (guide_end - guide_start)

    region_extract_start = region_guide_start - (guide_start - extract_start)
    region_extract_end = region_guide_end + (extract_end - guide_end)

    if region_extract_start < 0 or region_extract_end > len(region):
        return None, None, None, None

    seq = region[region_extract_start:region_extract_end].upper()
    codon_offset_start = guide_start - extract_start
    codon_offset_end = extract_end - guide_end

    return seq, codon_offset_start, codon_offset_end, cds_strand
