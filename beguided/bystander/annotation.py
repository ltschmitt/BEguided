"""Genome annotation operations for bystander analysis."""

import re

from intervaltree import IntervalTree
from Bio.Seq import Seq


# Reserved key in the returned dict for the transcript → codon index.
# A chromosome literally named this would collide; guard against it.
_TX_CODONS_KEY = "__tx_codons__"


def load_gtf_to_interval_tree(gtf_file, select_tag="MANE_Select"):
    """Load GTF file and build interval trees for the canonical CDS per chromosome.

    Only CDS records whose attributes contain ``select_tag`` are kept. The
    default ``"MANE_Select"`` matches the Ensembl human annotation; the EBI
    mouse annotation has no MANE tags and uses ``"Ensembl_canonical"``
    instead, so callers pass that per genome.

    Each interval stores ``(strand, transcript_id, gene_name)``. ``gene_name``
    is extracted from the ``gene_name`` attribute and falls back to ``""``
    when the attribute is absent.

    The returned dict also contains a reserved key ``"__tx_codons__"`` mapping
    ``transcript_id`` to a precomputed cumulative-codon index used by
    :func:`absolute_codon_number`.

    Args:
        gtf_file: Path to the GTF annotation file.
        select_tag: Attribute tag identifying the canonical transcript to keep
            (e.g. ``"MANE_Select"`` or ``"Ensembl_canonical"``).

    Returns:
        dict: {chromosome: IntervalTree(...), "__tx_codons__": {transcript_id: {...}}}
    """
    trees = {}
    cds_count = 0
    skipped_count = 0

    # Collect CDS records per transcript for the codon index.
    tx_cds = {}

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
            if select_tag not in attributes:
                skipped_count += 1
                continue

            chrom = parts[0]
            start = int(parts[3])  # 1-based, inclusive
            end = int(parts[4])    # 1-based, inclusive
            strand = parts[6]

            transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
            transcript_id = transcript_match.group(1) if transcript_match else "unknown"

            gene_match = re.search(r'gene_name "([^"]+)"', attributes)
            gene_name = gene_match.group(1) if gene_match else ""

            exon_match = re.search(r'exon_number "(\d+)"', attributes)
            exon_number = int(exon_match.group(1)) if exon_match else None

            if chrom not in trees:
                trees[chrom] = IntervalTree()

            # intervaltree uses 0-based, half-open [start, end)
            begin_0 = start - 1
            end_0 = end
            trees[chrom][begin_0:end_0] = (strand, transcript_id, gene_name)
            cds_count += 1

            # Collect for codon index
            if transcript_id not in tx_cds:
                tx_cds[transcript_id] = {"strand": strand, "exons": []}
            tx_cds[transcript_id]["exons"].append((begin_0, end_0, exon_number))

    # Build cumulative codon index per transcript.
    trees[_TX_CODONS_KEY] = _build_tx_codon_index(tx_cds)

    return trees


def _build_tx_codon_index(tx_cds):
    """Build the transcript → cumulative-codon lookup.

    Args:
        tx_cds: {transcript_id: {"strand": str, "exons": [(begin_0, end_0, exon_number_or_None), ...]}}

    Returns:
        {transcript_id: {"strand": str, "exons": [...], "cum_nt": np.array}}
    """
    import numpy as np

    index = {}
    for tx_id, info in tx_cds.items():
        strand = info["strand"]
        exons = info["exons"]

        # Sort exons in transcript order (5'→3' on mRNA).
        # Use exon_number when present for every exon; otherwise fall back
        # to genomic order.
        if all(e[2] is not None for e in exons):
            sorted_exons = sorted(exons, key=lambda e: e[2])
        else:
            if strand == "-":
                sorted_exons = sorted(exons, key=lambda e: -e[0])
            else:
                sorted_exons = sorted(exons, key=lambda e: e[0])

        lengths = [e[1] - e[0] for e in sorted_exons]
        cum_nt = np.cumsum(lengths, dtype=np.int64)

        index[tx_id] = {
            "strand": strand,
            "exons": sorted_exons,
            "cum_nt": cum_nt,
        }

    return index


def absolute_codon_number(transcript_id, genomic_pos, interval_trees):
    """Return the 1-based absolute codon number of *genomic_pos* in the protein.

    Uses the precomputed cumulative-codon index stored under
    ``"__tx_codons__"`` in *interval_trees* (built by
    :func:`load_gtf_to_interval_tree`).

    Args:
        transcript_id: Transcript identifier from the GTF.
        genomic_pos: 0-based genomic coordinate of the edited base.
        interval_trees: The dict returned by :func:`load_gtf_to_interval_tree`.

    Returns:
        int: 1-based codon number, or ``None`` if the transcript is not
        indexed or *genomic_pos* falls outside every exon.
    """
    mrna_pos = cumulative_mrna_pos(transcript_id, genomic_pos, interval_trees)
    if mrna_pos is None:
        return None
    return mrna_pos // 3 + 1


def cumulative_mrna_pos(transcript_id, genomic_pos, interval_trees):
    """Return the 0-based cumulative mRNA position of *genomic_pos*.

    The cumulative mRNA position is the 0-based offset (5'→3') of the
    given base within the spliced mRNA — i.e. zero at the first base
    of the first exon, regardless of CDS phase. The codon containing
    *genomic_pos* is ``cumulative_mrna_pos // 3 + 1`` (see
    :func:`absolute_codon_number`).

    This is the right frame of reference for codon-aligned extraction:
    the CDS is split across multiple GTF exons that may each have a
    non-zero ``phase`` (the number of bases that complete the previous
    codon at the 5' end of the exon). Anchoring the codon frame to a
    single exon's start/end gives wrong codon boundaries for any exon
    with ``phase != 0``.

    Args:
        transcript_id: Transcript identifier from the GTF.
        genomic_pos: 0-based genomic coordinate of the edited base.
        interval_trees: The dict returned by :func:`load_gtf_to_interval_tree`.

    Returns:
        int: 0-based cumulative mRNA position, or ``None`` if the
        transcript is not indexed or *genomic_pos* falls outside every
        exon.
    """
    tx_index = interval_trees.get(_TX_CODONS_KEY, {})
    if transcript_id not in tx_index:
        return None

    entry = tx_index[transcript_id]
    strand = entry["strand"]
    exons = entry["exons"]
    cum_nt = entry["cum_nt"]

    for i, (begin, end, _exon_no) in enumerate(exons):
        if begin <= genomic_pos < end:
            if strand == "-":
                pos_in_exon = (end - 1) - genomic_pos
            else:
                pos_in_exon = genomic_pos - begin
            prev_cum = int(cum_nt[i - 1]) if i > 0 else 0
            return prev_cum + pos_in_exon

    return None


def mrna_pos_to_genomic(transcript_id, mrna_pos, interval_trees):
    """Inverse of :func:`cumulative_mrna_pos`.

    Given a 0-based cumulative mRNA position, return the 0-based genomic
    coordinate of that base. Useful for converting codon-aligned
    boundaries back to genomic coordinates when extracting sequence.

    Args:
        transcript_id: Transcript identifier from the GTF.
        mrna_pos: 0-based cumulative mRNA position.
        interval_trees: The dict returned by :func:`load_gtf_to_interval_tree`.

    Returns:
        tuple: ``(chrom, genomic_pos_0based)`` if the position falls inside
        the transcript, or ``None`` if the transcript is not indexed or
        ``mrna_pos`` is outside the spliced mRNA.
    """
    tx_index = interval_trees.get(_TX_CODONS_KEY, {})
    if transcript_id not in tx_index:
        return None

    entry = tx_index[transcript_id]
    strand = entry["strand"]
    exons = entry["exons"]
    cum_nt = entry["cum_nt"]

    if mrna_pos < 0:
        return None
    for i, (begin, end, _exon_no) in enumerate(exons):
        length = int(cum_nt[i]) - (int(cum_nt[i - 1]) if i > 0 else 0)
        if mrna_pos < int(cum_nt[i]):
            offset_in_exon = mrna_pos - (int(cum_nt[i - 1]) if i > 0 else 0)
            if strand == "-":
                genomic_pos = (end - 1) - offset_in_exon
            else:
                genomic_pos = begin + offset_in_exon
            return genomic_pos

    return None


def get_cds_overlap(chrom, guide_start, guide_end, interval_trees):
    """Find CDS overlaps for a guide region.

    Returns:
        tuple: (list of (cds_start, cds_end, strand, transcript_id, gene_name), is_multiple)
    """
    if chrom not in interval_trees or chrom == _TX_CODONS_KEY:
        return [], False

    overlaps = interval_trees[chrom][guide_start:guide_end]
    if not overlaps:
        return [], False

    cds_list = []
    for interval in overlaps:
        strand, transcript_id, gene_name = interval.data
        cds_list.append((interval.begin, interval.end, strand, transcript_id, gene_name))

    is_multiple = len(cds_list) > 1
    return cds_list, is_multiple


def extract_coding_sequence(chrom, guide_start, guide_end, cds_info, genome,
                              interval_trees=None):
    """Extract the coding sequence covering the guide region with flanking.

    Codon boundaries are aligned to the **transcript** reading frame
    (using :func:`cumulative_mrna_pos`), not the start/end of the
    single exon that the guide happens to overlap. This matters when
    the overlapping exon has a non-zero ``phase`` (the first
    ``phase`` bases of that exon complete a codon started in the
    previous exon) — anchoring to the exon would shift every codon
    boundary in that exon.

    The extract is widened to the next codon boundary upstream and
    downstream of the guide so that ``coding_seq[0]`` is always the
    first base of a codon and ``len(coding_seq) % 3 == 0``. If the
    codon-aligned window crosses an exon junction, the splice is
    joined from the genome.

    Args:
        chrom: Chromosome name.
        guide_start: 0-based genomic start of guide.
        guide_end: 0-based genomic end of guide (exclusive).
        cds_info: List of (cds_start, cds_end, cds_strand, transcript_id, gene_name).
        genome: Dict-like ``{chrom: sequence_string}`` supporting
            ``genome[chrom][start:end]``.
        interval_trees: Dict from :func:`load_gtf_to_interval_tree`. Used
            to look up the transcript's exon index. If ``None``, falls
            back to single-exon frame alignment (legacy behaviour).

    Returns:
        tuple: ``(coding_seq, codon_offset_start, codon_offset_end, cds_strand)``.
        ``coding_seq`` is ``None`` if the window cannot be reconstructed.
    """
    if not cds_info:
        return None, None, None, None

    cds_start, cds_end, cds_strand, transcript_id, _gene_name = cds_info[0]

    bounds = _codon_aligned_bounds(
        transcript_id, guide_start, guide_end, cds_start, cds_end,
        interval_trees,
    )
    if bounds is None:
        return None, None, None, None
    extract_start, extract_end = bounds

    seq = _splice_genome(
        genome, interval_trees, chrom, extract_start, extract_end,
        transcript_id=transcript_id,
    )
    if seq is None:
        return None, None, None, None
    seq = seq.upper()

    codon_offset_start = guide_start - extract_start
    codon_offset_end = extract_end - guide_end

    return seq, codon_offset_start, codon_offset_end, cds_strand


def extract_coding_sequence_from_region(chrom, guide_start, guide_end,
                                          cds_info, region, guide_offset,
                                          genome=None, interval_trees=None):
    """Extract coding sequence from a pre-reconstructed region string.

    The *region* string is built by ``reconstruct_genomic_region`` and
    covers the guide + PAM + a few flanking bases.  *guide_offset*
    tells where the guide (not the PAM) starts inside *region*.

    Codon alignment uses the **transcript** reading frame (see
    :func:`extract_coding_sequence`). The function first tries to
    satisfy the codon-aligned window from ``region`` alone; if the
    window extends beyond the region (e.g. the codon containing the
    edit crosses an exon junction), it falls back to splicing from
    ``genome`` when ``genome`` is provided. Returns ``None`` if the
    window cannot be reconstructed.

    Args:
        chrom: Chromosome name.
        guide_start: 0-based genomic start of guide.
        guide_end: 0-based genomic end of guide (exclusive).
        cds_info: List of (cds_start, cds_end, cds_strand, transcript_id, gene_name).
        region: Reconstructed genomic region string.
        guide_offset: Index in *region* where the guide starts.
        genome: Optional ``{chrom: sequence_string}`` used as a
            fallback when the codon-aligned window exceeds ``region``.
        interval_trees: Dict from :func:`load_gtf_to_interval_tree`.

    Returns:
        Same tuple as :func:`extract_coding_sequence`.
    """
    if not cds_info:
        return None, None, None, None, None, None

    cds_start, cds_end, cds_strand, transcript_id, _gene_name = cds_info[0]

    bounds = _codon_aligned_bounds(
        transcript_id, guide_start, guide_end, cds_start, cds_end,
        interval_trees,
    )
    if bounds is None:
        return None, None, None, None, None, None
    extract_start, extract_end = bounds

    region_guide_start = guide_offset
    region_guide_end = guide_offset + (guide_end - guide_start)

    region_extract_start = region_guide_start - (guide_start - extract_start)
    region_extract_end = region_guide_end + (extract_end - guide_end)

    if region_extract_start >= 0 and region_extract_end <= len(region):
        seq = region[region_extract_start:region_extract_end].upper()
        codon_offset_start = guide_start - extract_start
        codon_offset_end = extract_end - guide_end
    elif genome is not None and chrom in genome:
        seq = _splice_genome(
            genome, interval_trees, chrom, extract_start, extract_end,
            transcript_id=transcript_id,
        )
        if seq is None:
            return None, None, None, None, None, None
        seq = seq.upper()
        codon_offset_start = _splice_aware_codon_offset(
            interval_trees, transcript_id, extract_start, extract_end,
            guide_start, guide_end, cds_strand,
        )
        codon_offset_end = codon_offset_start + (guide_end - guide_start)
    else:
        return None, None, None, None, None, None

    return seq, codon_offset_start, codon_offset_end, cds_strand, extract_start, extract_end


def _splice_aware_codon_offset(interval_trees, transcript_id, extract_start,
                               extract_end, guide_start, guide_end,
                               cds_strand):
    """Return a codon_offset_start that works for spliced coding_seq.

    For a spliced extract (introns removed), the simple formula
    ``guide_start - extract_start`` does not give the right offset in
    ``coding_seq`` because bases are skipped between ``extract_start``
    and the guide's first CDS position.

    The codon_offset value used by :func:`analyze_guide` is the offset
    of the guide's **first base** in ``coding_seq``. When the guide
    starts in the UTR/intron (outside coding_seq), this is
    conceptual — it's the position where GUIDE[1] *would* be if the
    extract included the upstream bases.

    We compute it by counting how many CDS bases lie between the
    extract's first CDS position and the guide's first CDS position.
    """
    tree = interval_trees.get(_TX_CODONS_KEY, {}).get(transcript_id)
    if tree is None:
        return guide_start - extract_start

    exons = tree["exons"]
    strand = tree["strand"]
    pieces = [(begin, end) for begin, end, _phase in exons]
    if strand == "-":
        pieces = list(reversed(pieces))

    # Walk the spliced CDS in mRNA order; track the coding_seq offset
    # at each base. The guide's first CDS base sets the anchor.
    guide_first_cds_pos = None
    cum_offset = 0
    for begin, end in pieces:
        if end <= extract_start:
            continue
        if begin >= extract_end:
            break
        iv_start = max(begin, extract_start)
        iv_end = min(end, extract_end)
        # If the guide starts (or extends) into this piece, anchor.
        for g_pos in range(guide_start, guide_end):
            if iv_start <= g_pos < iv_end:
                guide_first_cds_pos = cum_offset + (g_pos - iv_start)
                break
        if guide_first_cds_pos is not None:
            break
        cum_offset += iv_end - iv_start

    if guide_first_cds_pos is None:
        return guide_start - extract_start

    # The codon_offset formula in analyze_guide assumes the GUIDE is
    # entirely inside coding_seq. When the GUIDE has bases outside
    # coding_seq (e.g. in the UTR or in an intron that we spliced
    # out), we need to shift the conceptual GUIDE[1] position. We
    # count how many GUIDE bases lie outside any CDS exon piece of
    # this transcript (these bases are missing from coding_seq) and
    # subtract that from guide_first_cds_pos to land at the
    # conceptual GUIDE[1] position in coding_seq.
    out = 0
    for g_pos in range(guide_start, guide_end):
        in_exon = False
        for begin, end in pieces:
            if begin <= g_pos < end:
                in_exon = True
                break
        if not in_exon:
            out += 1
        else:
            break
    return guide_first_cds_pos - out


def _codon_aligned_bounds(transcript_id, guide_start, guide_end,
                          cds_start, cds_end, interval_trees):
    """Return codon-aligned (extract_start, extract_end) in genomic coords.

    Uses the transcript's cumulative mRNA position so that the codon
    frame reflects the full CDS, not just the local exon. The window
    is **not** clamped to ``[cds_start, cds_end]`` — those are the
    bounds of the single exon that overlaps the guide, and the codon-
    aligned extract frequently has to cross exon boundaries (e.g. the
    codon containing the edit may share bases with the previous/next
    exon). The caller is expected to splice from the genome and the
    CDS interval tree to skip introns.

    Guide positions outside the CDS (e.g. when the guide extends into
    the UTR) are clamped to the nearest CDS base before the codon
    lookup, so that ``cumulative_mrna_pos`` can resolve.

    Returns ``None`` if the guide falls entirely outside any CDS exon.
    """
    # When the transcript index isn't available, fall back to single-exon
    # alignment. This is the legacy behaviour and remains correct only
    # when the overlapping exon has phase 0.
    if interval_trees is None or transcript_id is None:
        return _legacy_exon_aligned_bounds(
            guide_start, guide_end, cds_start, cds_end,
        )

    # Clamp guide positions to CDS so guides that extend into the UTR
    # still resolve to a position where cumulative_mrna_pos works.
    clamped_a = max(min(guide_start, cds_end - 1), cds_start)
    clamped_b = max(min(guide_end - 1, cds_end - 1), cds_start)

    cum_a = cumulative_mrna_pos(transcript_id, clamped_a, interval_trees)
    cum_b = cumulative_mrna_pos(transcript_id, clamped_b, interval_trees)
    if cum_a is None or cum_b is None:
        return None

    cum_min, cum_max = (cum_a, cum_b) if cum_a <= cum_b else (cum_b, cum_a)
    aligned_cum_start = (cum_min // 3) * 3
    aligned_cum_end = (cum_max // 3 + 1) * 3

    # Translate back to genomic coords; the window may span multiple exons.
    start_geno = mrna_pos_to_genomic(
        transcript_id, aligned_cum_start, interval_trees,
    )
    end_geno = mrna_pos_to_genomic(
        transcript_id, aligned_cum_end - 1, interval_trees,
    )
    if start_geno is None or end_geno is None:
        return None
    # Genomic coords are always returned low→high regardless of strand;
    # mRNA order is reflected inside the sequence itself.
    low = min(start_geno, end_geno)
    high = max(start_geno, end_geno)
    # Allow the extract to extend past the local exon's CDS bounds so
    # the codon-aligned window can span exon junctions. Splicing will
    # skip introns at extraction time.
    extract_start = low
    extract_end = high + 1

    if extract_end <= extract_start:
        return None
    return extract_start, extract_end


def _legacy_exon_aligned_bounds(guide_start, guide_end, cds_start, cds_end):
    """Legacy codon alignment anchored to a single exon.

    Retained as a fallback when the GTF-derived transcript index is not
    available. Only correct when the overlapping exon has phase 0.
    """
    offset_from_cds_start = guide_start - cds_start
    needed_upstream = offset_from_cds_start % 3
    extract_start = max(guide_start - needed_upstream, cds_start)
    offset_from_cds_end = guide_end - extract_start
    needed_downstream = 3 - (offset_from_cds_end % 3)
    extract_end = min(guide_end + needed_downstream, cds_end)
    if extract_end <= extract_start:
        return None
    return extract_start, extract_end


def _splice_genome(genome, interval_trees, chrom, extract_start, extract_end,
                   transcript_id=None):
    """Return the genomic ``+``-strand slice ``[extract_start, extract_end)``.

    If the slice crosses an exon junction of the relevant transcript,
    only the CDS exon bases are kept (introns are skipped). The pieces
    are joined in **genomic** order, so the returned string is always
    on the + strand — callers that need the mRNA (e.g. for translation)
    must reverse-complement the result themselves.

    When *transcript_id* is given, only intervals belonging to that
    transcript are kept. This is important because the chromosome's
    interval tree may contain CDS intervals from multiple genes that
    overlap the same genomic window; splicing across transcripts would
    produce nonsense sequence.

    The CDS exon list is taken from the chromosome's interval tree
    stored on *interval_trees*; if no CDS overlap is found the slice
    is returned verbatim (legacy behaviour).

    Returns ``None`` if the chromosome is missing or the slice is empty.
    """
    if chrom not in genome:
        return None
    if extract_end <= extract_start:
        return None

    tree = interval_trees.get(chrom) if interval_trees is not None else None
    if tree is None:
        return genome[chrom][extract_start:extract_end]

    overlapping = sorted(tree[extract_start:extract_end], key=lambda i: i.begin)
    if not overlapping:
        return None

    pieces = []
    for interval in overlapping:
        if transcript_id is not None:
            data = interval.data
            if not data or data[1] != transcript_id:
                continue
        iv_start = max(interval.begin, extract_start)
        iv_end = min(interval.end, extract_end)
        if iv_end <= iv_start:
            continue
        pieces.append(genome[chrom][iv_start:iv_end])
    return "".join(pieces) if pieces else None
