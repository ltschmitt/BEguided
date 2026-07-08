"""Shared test fixtures for BEguided."""

import os
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def sample_genome_dict():
    """Small test genome dictionary."""
    return {
        "chr1": "ACGAAGACTACGAGTAGCCGTTGGATCGTACCTTCGATCTAGCGATCATCAGAGATATAGGGAACGGTAGCTAGGTAGCTAGCTAGCTAGGGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        "chr2": "AGCAGATAGCCCCATTCTGGACTAGCGGACACCGATCGACTTGCCTGGTAGCCGGTATCTGGGTACTGGGATTCCGAGGGACTGCTCCGAAGATTCGCTGCAGAACAGGAGCCGATAGC",
    }


@pytest.fixture
def valid_preset_tsv():
    """Content of a valid base editor profiles TSV file."""
    return (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\thigh_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
        "saABE8e\tABE\t5\t13\t5\t13\tNNGRRT\tsaCas9\tNature2020\n"
        "ABE9\tABE\t4\t7\t5\t6\tNGG\tSpCas9\tChen2023\n"
    )


@pytest.fixture
def duplicate_preset_tsv():
    """Content with duplicate base editor names."""
    return (
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\thigh_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
        "ABE8e\tABE\t4\t7\t5\t6\tNGG\tSpCas9\tChen2023\n"
    )


@pytest.fixture
def abe8e_preset_config(valid_preset_tsv):
    """ABE8e base editor configuration (sense strand A->G)."""
    from beguided.editing.presets import load_base_editor_profiles
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(valid_preset_tsv)
        temp_path = f.name
    profiles = load_base_editor_profiles(temp_path)
    os.unlink(temp_path)
    return profiles[("A", "G")]["ABE8e"]


@pytest.fixture
def abe8e_antisense_config(valid_preset_tsv):
    """ABE8e base editor configuration (antisense strand T->C)."""
    from beguided.editing.presets import load_base_editor_profiles
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(valid_preset_tsv)
        temp_path = f.name
    profiles = load_base_editor_profiles(temp_path)
    os.unlink(temp_path)
    return profiles[("T", "C")]["ABE8e"]


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a temporary FASTA file for testing."""
    fasta_content = ">chr1\nACGAAGACTACGAGTAGCCGTTGGATCGTACCTTCGATCTAGCGATCATCAGAGATATAGGGAACGGTAGCTAGGTAGCTAGCTAGCTAGGGGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n>chr2\nAGCAGATAGCCCCATTCTGGACTAGCGGACACCGATCGACTTGCCTGGTAGCCGGTATCTGGGTACTGGGATTCCGAGGGACTGCTCCGAAGATTCGCTGCAGAACAGGAGCCGATAGC\n"
    fasta_file = tmp_path / "test_genome.fa"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)


@pytest.fixture
def sample_targets_tsv(tmp_path):
    """Create a temporary targets TSV file."""
    tsv_content = "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\nchr1\t75\tA\tG\nchr2\t30\tT\tC\nchr2\t80\tA\tG\n"
    tsv_file = tmp_path / "test_targets.tsv"
    tsv_file.write_text(tsv_content)
    return str(tsv_file)


@pytest.fixture
def synthetic_beguided_output():
    """Create synthetic beguided output DataFrame."""
    import pandas as pd
    return pd.DataFrame({
        "chromosome": ["chr1", "chr1", "chr2", "chr2"],
        "position": [100, 200, 300, 400],
        "target_base": ["A", "T", "A", "T"],
        "desired_base": ["G", "C", "G", "C"],
        "genome_base": ["G", "C", "G", "C"],
        "strand_used": ["+", "-", "+", "-"],
        "editing_position": [5, 6, 7, 8],
        "bystander_positions": ["3,7", "", "5", ""],
        "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT",
                               "GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC"],
        "pam_sequence": ["AGG", "TGG", "GGG", "CGG"],
        "base_editor": ["ABE8e", "ABE8e", "saABE8e", "saABE8e"]
    })


@pytest.fixture
def synthetic_bystander_output():
    """Create synthetic bystander_scanner output DataFrame."""
    import pandas as pd
    return pd.DataFrame({
        "chromosome": ["chr1", "chr1", "chr2", "chr2", "chr3"],
        "position": [100, 200, 300, 400, 500],
        "target_base": ["A", "T", "A", "T", "A"],
        "desired_base": ["G", "C", "G", "C", "G"],
        "genome_base": ["G", "C", "G", "C", "G"],
        "strand_used": ["+", "-", "+", "-", "+"],
        "editing_position": [5, 6, 7, 8, 9],
        "bystander_positions": ["3,7", "", "5", "", ""],
        "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT",
                               "GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC",
                               "NNNNNNNNNNNNNNNNNNNN"],
        "pam_sequence": ["AGG", "TGG", "GGG", "CGG", "AGG"],
        "base_editor": ["ABE8e", "ABE8e", "saABE8e", "saABE8e", "ABE8e"],
        "bystander_effect": ["missense", "no_bystanders", "silent", "no_bystanders", "no_bystanders"],
        "simulated_outcome": ["AAAAAAAGAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT",
                              "GGGGGAGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC",
                              "NNNNNNNNNNNNNNNNNNNN"],
        "aa_wanted": ["K", "F", "G", "P", "X"],
        "aa_bystander": ["R", "F", "G", "P", "X"]
    })


@pytest.fixture
def synthetic_beguidedoff_output():
    """Create synthetic BEguidedOFF output DataFrame."""
    import pandas as pd
    return pd.DataFrame({
        "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT",
                               "GGGGGGGGGGGGGGGGGGGG"],
        "base_editor": ["ABE8e", "ABE8e", "saABE8e"],
        "ot_rank": [1, 1, 1],
        "ot_chromosome": ["chr1", "chr1", "chr2"],
        "ot_position": [1000, 2000, 3000],
        "ot_strand": ["+", "-", "+"],
        "ot_sequence": ["AAAAAAAAAAAAAAAAAATA", "TTTTTTTTTTTTTTTTTTTT",
                        "GGGGGGGGGGGGGGGGGGGT"],
        "ot_pam": ["AGG", "TGG", "GGG"],
        "n_mismatches": [1, 0, 1],
        "mismatch_positions": ["19", "", "20"],
        "score": [1.2, 0.0, 1.2],
        "n_offtargets_below_threshold": [5, 1, 3],
    })


# BEguidedOFF npz fixture helper (inline in tests to avoid FASTA coupling)
