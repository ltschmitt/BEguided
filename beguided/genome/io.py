"""Genome I/O operations."""

from Bio import SeqIO


def load_genome(genome_file):
    """Read genome file and return dictionary of chromosomes.

    Args:
        genome_file: Path to FASTA file.

    Returns:
        dict: chromosome_id -> sequence string.
    """
    genome_dict = {}
    with open(genome_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            chrom_id = str(record.id).split()[0]
            genome_dict[chrom_id] = str(record.seq)
    return genome_dict
