"""Shared data models for BEguided."""
from dataclasses import dataclass, field
from typing import Dict, Optional
import re


@dataclass
class BaseEditorConfig:
    """Configuration for a single base editor profile on a specific strand."""

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
    base_editor: str = ""
    genomic_context_upstream: str = ""
    genomic_context_downstream: str = ""
    guide_score: Optional[float] = None
    gc_content: Optional[float] = None
    homopolymer: Optional[bool] = None
    trinuc_edit_score: Optional[float] = None
    extra_columns: Dict = field(default_factory=dict, repr=False)


@dataclass
class EditingOutcome:
    """Outcome of simulating edits on a coding sequence."""

    effect: str
    simulated_outcome: str
    aa_wanted: str
    aa_bystander: str


# Marker strings placed in ``guide_rna_sequence`` for targets that end up in
# the results without a usable guide. Two of them are produced by the pipeline
# for the ``keep_targets_without_guides`` config option (added at the very
# last pipeline step); the third is produced by the web UI's dynamic result
# filters when a filter removes every guide of a target. The trailing
# ``**`` prefix matches the pipeline's existing "invalid guide" convention
# (see ``run_pipeline`` step-1 logging).
NO_GUIDE_FOUND = "**No guide possible**"
BASE_CHANGE_INCOMPATIBLE = "**Base change incompatible**"
GUIDES_REMOVED_BY_FILTER = "**Guides removed by filter**"
PLACEHOLDER_MARKERS = frozenset(
    {NO_GUIDE_FOUND, BASE_CHANGE_INCOMPATIBLE, GUIDES_REMOVED_BY_FILTER}
)


CORE_COLUMN_ORDER = [
    "chromosome",
    "position",
    "target_base",
    "desired_base",
    "genome_base",
    "gene_edit_desired",
    "base_editor",
    "guide_rna_sequence",
    "pam_sequence",
    "strand_used",
    "editing_position",
    "bystander_positions",
    "bystander_effect",
    "guide_score",
    "gc_content",
    "homopolymer",
    "trinuc_edit_score",
]
