"""Sequence manipulation utilities."""

from Bio.Seq import Seq


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence.

    Handles both upper and lower case; leaves non-ACGT chars untouched.
    """
    comp = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }
    return "".join(comp.get(base, base) for base in reversed(seq))


def apply_reverse_complement_if_needed(guide_seq, pam_seq, should_rc):
    """Reverse complement guide and PAM if antisense strand is used.

    Args:
        guide_seq: Guide DNA sequence string.
        pam_seq: PAM DNA sequence string.
        should_rc: Bool; True if reverse complement is needed.

    Returns:
        tuple: (final_guide_seq, final_pam_seq)
    """
    if should_rc:
        final_guide = str(Seq(guide_seq).reverse_complement())
        final_pam = str(Seq(pam_seq).reverse_complement())
    else:
        final_guide = guide_seq
        final_pam = pam_seq
    return final_guide, final_pam


def insert_mutation(seq, position, base, strand, guide_length):
    """Insert a mutation into a guide sequence.

    Args:
        seq: Original guide sequence string.
        position: 1-based editing position from 5' end.
        base: Base to insert.
        strand: '+' or '-'.
        guide_length: Length of the guide.

    Returns:
        Mutated sequence string.
    """
    seq_list = list(seq)
    if strand == "-":
        seq_list[guide_length - position] = base
    else:
        seq_list[position - 1] = base
    return "".join(seq_list)
