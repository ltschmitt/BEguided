"""Web application configuration and preset helpers."""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import pandas as pd

from beguided.editing.presets import load_presets


def _env_bool(name: str, default: bool) -> bool:
    """Return the boolean value of an env var, falling back to ``default``."""
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() not in ("0", "false", "no", "off", "")

_DEFAULT_PRESETS = str(Path(__file__).resolve().parents[2] / "be_presets.tsv")
_DEFAULT_GENOME = str(
    Path(__file__).resolve().parents[2]
    / "genome"
    / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
)
_DEFAULT_GTF = str(
    Path(__file__).resolve().parents[2] / "genome" / "Homo_sapiens.GRCh38.112.gtf"
)


@dataclass
class WebConfig:
    """Configuration for the BEguided web server."""

    presets_file: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_PRESETS_FILE", _DEFAULT_PRESETS)
    )
    genome_fasta: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_GENOME_FASTA", _DEFAULT_GENOME)
    )
    gtf_path: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_GTF_PATH", _DEFAULT_GTF)
    )
    beguidedoff_genome: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_BEGUIDEDOFF_GENOME",
            str(Path(__file__).resolve().parents[2] / "mini_genome.npz"),
        )
    )
    count_threshold_score: float = 15.0
    remove_guides_score: float = 5.0
    mismatch_tolerance: int = 5
    remove_guides_n_offtargets: int = 100
    max_input_targets: int = 1000
    preload_genome: bool = field(
        default_factory=lambda: _env_bool("BEGUIDED_PRELOAD_GENOME", True)
    )
    genome_label: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_GENOME_LABEL", "").strip()
    )
    preload_gtf: bool = field(
        default_factory=lambda: _env_bool("BEGUIDED_PRELOAD_GTF", True)
    )
    # Scheduler / safety knobs
    db_path: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_DB_PATH", "/var/lib/beguided/jobs.db")
    )
    queue_max: int = 8
    cooldown_seconds: float = 60.0
    max_upload_mb: int = 10

    @property
    def has_beguidedoff(self) -> bool:
        """Return True if BEguidedOFF appears to be installed/configured."""
        try:
            import beguidedoff  # noqa: F401
            return True
        except ImportError:
            return False


# Global singleton – replaced at app startup if CLI args override values
_config: WebConfig = WebConfig()


def set_config(cfg: WebConfig) -> None:
    """Override the global web config singleton."""
    global _config
    _config = cfg


def get_config() -> WebConfig:
    """Return the current global web config."""
    return _config


def load_preset_groups(presets_file: str | None = None):
    """Load presets and return grouped options for the dropdown.

    Args:
        presets_file: Path to presets TSV. Uses config default if None.

    Returns:
        dict: Mapping label -> list of preset names.
    """
    if presets_file is None:
        presets_file = _config.presets_file
    df = pd.read_csv(presets_file, sep="\t")

    abe_names = sorted(df[df["editing"] == "A-to-G"]["name"].tolist())
    cbe_names = sorted(df[df["editing"] == "C-to-T"]["name"].tolist())
    all_names = sorted(df["name"].tolist())

    groups = {
        "all": all_names,
        "all ABEs": abe_names,
        "all CBEs": cbe_names,
    }
    for name in all_names:
        groups[name] = [name]

    return groups


def get_preset_display_names(presets_file: str | None = None):
    """Return ordered list of dropdown labels."""
    groups = load_preset_groups(presets_file)
    return ["all", "all ABEs", "all CBEs"] + sorted(
        set(groups.keys()) - {"all", "all ABEs", "all CBEs"}
    )


def get_editing_type_for_presets(
    preset_names: List[str], presets_file: str | None = None
) -> dict:
    """Return a mapping of editing type -> preset names.

    Args:
        preset_names: List of preset names to look up.
        presets_file: Path to presets TSV. Uses config default if None.

    Returns:
        dict: {"A-to-G": [...], "C-to-T": [...]}
    """
    if presets_file is None:
        presets_file = _config.presets_file
    df = pd.read_csv(presets_file, sep="\t")
    df = df[df["name"].isin(preset_names)]
    result = {}
    for editing, group in df.groupby("editing"):
        result[editing] = sorted(group["name"].tolist())
    return result
