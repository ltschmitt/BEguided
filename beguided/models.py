"""Shared data models for BEguided."""
from dataclasses import dataclass, field
from typing import Dict, Optional
import re


@dataclass
class PresetConfig:
    """Configuration for a single base editor preset on a specific strand."""

    name: str
    pam_pattern: str
    pam_regex: re.Pattern
    pam_length: int
    guide_length: int
    high_edit_start: int
    high_edit_end: int
    general_edit_start: int
    general_edit_end: int
    strand: str
    reverse_complement: bool
    target_base: str
    desired_base: str


@dataclass
class GuideResult:
    """Result dictionary for a generated guide."""

    chromosome: str
    position: int
    target_base: str
    desired_base: str
    guide_rna_sequence: str
    genome_base: str
    pam_sequence: str
    editing_position: int
    strand_used: str
    bystander_positions: str = ""
    be_preset: str = ""
    genomic_context_upstream: str = ""
    genomic_context_downstream: str = ""
    guide_score: Optional[float] = None
    extra_columns: Dict = field(default_factory=dict, repr=False)


@dataclass
class EditingOutcome:
    """Outcome of simulating edits on a coding sequence."""

    effect: str
    simulated_outcome: str
    aa_wanted: str
    aa_bystander: str


CORE_COLUMN_ORDER = [
    "chromosome",
    "position",
    "target_base",
    "desired_base",
    "genome_base",
    "gene_edit_desired",
    "be_preset",
    "guide_rna_sequence",
    "pam_sequence",
    "strand_used",
    "editing_position",
    "bystander_positions",
    "bystander_effect",
    "guide_score",
]
