"""Web application configuration and base editor helpers."""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import pandas as pd

from beguided.editing.presets import load_base_editor_profiles


def _env_bool(name: str, default: bool) -> bool:
    """Return the boolean value of an env var, falling back to ``default``."""
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() not in ("0", "false", "no", "off", "")


def _csv_env(name: str) -> List[str]:
    """Read a comma-separated env var into a list of stripped non-empty strings.

    Returns ``[]`` if the env var is unset or empty. Used for the
    ``BEGUIDED_TRUSTED_PROXY_IPS`` knob so the operator can override the
    defaults from the shell environment (or via ``docker-compose.override.yml``)
    without touching the code.
    """
    raw = os.getenv(name, "").strip()
    if not raw:
        return []
    return [s.strip() for s in raw.split(",") if s.strip()]



_DEFAULT_PRESETS = str(Path(__file__).resolve().parents[2] / "configs" / "base_editor_profiles.tsv")
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

    base_editor_profiles_file: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_BASE_EDITOR_PROFILES_FILE", _DEFAULT_PRESETS
        )
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
    gc_content_filter: str = "off"
    filter_homopolymers: bool = False
    trinuc_edit_score_min: float = None  # None = disabled (default)
    preload_genome: bool = field(
        default_factory=lambda: _env_bool("BEGUIDED_PRELOAD_GENOME", True)
    )
    genome_label: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_GENOME_LABEL", "").strip()
    )
    preload_gtf: bool = field(
        default_factory=lambda: _env_bool("BEGUIDED_PRELOAD_GTF", True)
    )
    default_base_editor_profiles_file: str = ""
    max_profile_upload_mb: int = 1

    # Scheduler / safety knobs
    db_path: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_DB_PATH", "/var/lib/beguided/jobs.db")
    )
    queue_max: int = field(
        default_factory=lambda: int(os.getenv("BEGUIDED_QUEUE_MAX", "8"))
    )
    cooldown_seconds: float = field(
        default_factory=lambda: float(os.getenv("BEGUIDED_COOLDOWN_SECONDS", "60"))
    )
    max_upload_mb: int = field(
        default_factory=lambda: int(os.getenv("BEGUIDED_MAX_UPLOAD_MB", "10"))
    )

    # Compliance / hardening knobs (see compliance_safety.md)
    job_retention_hours: int = 24
    trusted_proxy_ips: List[str] = field(
        default_factory=lambda: _csv_env("BEGUIDED_TRUSTED_PROXY_IPS") or ["127.0.0.1", "::1"]
    )
    log_ip_salt: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_LOG_IP_SALT", "").strip()
    )
    admin_token: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_ADMIN_TOKEN", "").strip()
    )
    legal_config_path: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_LEGAL_CONFIG", "").strip()
    )

    # Anonymous web-usage statistics (see beguided/web/usage_stats.py).
    # Enabled by default; the default path lives on the same persisted
    # volume as jobs.db so it survives container rebuilds. The data
    # contains no personal data, so it may be retained indefinitely.
    usage_stats_path: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_USAGE_STATS_PATH", "/var/lib/beguided/usage_stats.csv"
        )
    )
    usage_stats_enabled: bool = field(
        default_factory=lambda: _env_bool("BEGUIDED_USAGE_STATS_ENABLED", True)
    )

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
    if not cfg.default_base_editor_profiles_file:
        cfg.default_base_editor_profiles_file = cfg.base_editor_profiles_file
    _config = cfg


def get_config() -> WebConfig:
    """Return the current global web config."""
    return _config


def load_base_editor_groups(base_editor_profiles_file: str | None = None):
    """Load base editor profiles and return grouped options for the dropdown.

    Args:
        base_editor_profiles_file: Path to base editor profiles TSV. Uses config default if None.

    Returns:
        dict: Mapping label -> list of base editor names.
    """
    if base_editor_profiles_file is None:
        base_editor_profiles_file = _config.base_editor_profiles_file
    df = pd.read_csv(base_editor_profiles_file, sep="\t")

    abe_names = sorted(df[df["editing"] == "ABE"]["name"].tolist())
    cbe_names = sorted(df[df["editing"] == "CBE"]["name"].tolist())
    all_names = sorted(df["name"].tolist())

    groups = {
        "all": all_names,
        "all ABEs": abe_names,
        "all CBEs": cbe_names,
    }
    for name in all_names:
        groups[name] = [name]

    return groups


def get_base_editor_display_names(base_editor_profiles_file: str | None = None):
    """Return ordered list of dropdown labels."""
    groups = load_base_editor_groups(base_editor_profiles_file)
    return ["all", "all ABEs", "all CBEs"] + sorted(
        set(groups.keys()) - {"all", "all ABEs", "all CBEs"}
    )


def get_editing_type_for_base_editors(
    base_editor_names: List[str], base_editor_profiles_file: str | None = None
) -> dict:
    """Return a mapping of editing type -> base editor names.

    Args:
        base_editor_names: List of base editor names to look up.
        base_editor_profiles_file: Path to base editor profiles TSV. Uses config default if None.

    Returns:
        dict: {"ABE": [...], "CBE": [...]}
    """
    if base_editor_profiles_file is None:
        base_editor_profiles_file = _config.base_editor_profiles_file
    df = pd.read_csv(base_editor_profiles_file, sep="\t")
    df = df[df["name"].isin(base_editor_names)]
    result = {}
    for editing, group in df.groupby("editing"):
        result[editing] = sorted(group["name"].tolist())
    return result
