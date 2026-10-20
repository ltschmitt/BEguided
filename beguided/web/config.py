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
    / "Homo_sapiens.GRCh38.116.fa"
)
_DEFAULT_GTF = str(
    Path(__file__).resolve().parents[2] / "genome" / "Homo_sapiens.GRCh38.116.gtf"
)
# Mouse reference (GCA_000001635.9 = GRCm39) downloaded by scripts/download_genome.sh
# from the EBI ensemblorganisms tree. The EBI GTF uses ``Ensembl_canonical``
# tags (there are no MANE_Select tags for mouse).
_DEFAULT_MOUSE_GENOME = str(
    Path(__file__).resolve().parents[2]
    / "genome"
    / "Mus_musculus.GRCm39.2026_04.fa"
)
_DEFAULT_MOUSE_GTF = str(
    Path(__file__).resolve().parents[2] / "genome" / "Mus_musculus.GRCm39.2026_04.gtf"
)


@dataclass(frozen=True)
class GenomeOption:
    """A selectable reference genome (FASTA + GTF) for the web UI.

    ``gtf_tag`` is the CDS attribute tag marking the canonical transcript in
    the GTF (``"MANE_Select"`` for human Ensembl, ``"Ensembl_canonical"`` for
    the EBI mouse annotation). ``species`` is the scientific name shown in the
    UI selector (e.g. ``"Homo sapiens"``).
    """

    key: str
    label: str
    genome_fasta: str
    gtf_path: str
    gtf_tag: str = "MANE_Select"
    species: str = ""


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
    mouse_genome_fasta: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_MOUSE_GENOME_FASTA", _DEFAULT_MOUSE_GENOME
        )
    )
    mouse_gtf_path: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_MOUSE_GTF_PATH", _DEFAULT_MOUSE_GTF
        )
    )
    mouse_genome_label: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_MOUSE_GENOME_LABEL", "GRCm39"
        ).strip()
    )
    default_genome: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_DEFAULT_GENOME", "human").strip().lower()
    )
    beguidedoff_genome: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_BEGUIDEDOFF_GENOME",
            str(Path(__file__).resolve().parents[2]
                / "tests" / "fixtures" / "mini_genome.npz"),
        )
    )
    # Search mode for BEguidedOFF: "fast" (CPU FM-index filter; the default —
    # the web VM has no GPU) or "extensive" (GPU brute force).
    beguidedoff_mode: str = field(
        default_factory=lambda: os.getenv(
            "BEGUIDED_BEGUIDEDOFF_MODE", "fast"
        ).strip().lower()
    )
    # FM-index (.fmi) for fast mode. Empty -> derived from the .npz path
    # (same basename, .fmi suffix) when the job config is assembled.
    beguidedoff_index: str = field(
        default_factory=lambda: os.getenv("BEGUIDED_BEGUIDEDOFF_INDEX", "").strip()
    )
    count_threshold_score: float = 15.0
    remove_guides_score: float = 5.0
    mismatch_tolerance: int = 5
    remove_guides_n_offtargets: int = 100
    max_input_targets: int = field(
        default_factory=lambda: int(os.getenv("BEGUIDED_MAX_INPUT_TARGETS", "500"))
    )
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
    workers: int = field(
        default_factory=lambda: max(1, int(os.getenv("BEGUIDED_WORKERS", "1")))
    )
    max_upload_mb: int = field(
        default_factory=lambda: int(os.getenv("BEGUIDED_MAX_UPLOAD_MB", "10"))
    )

    # Compliance / hardening knobs (see compliance_safety.md)
    job_retention_hours: int = field(
        default_factory=lambda: int(os.getenv("BEGUIDED_JOB_RETENTION_HOURS", "24"))
    )
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
    def genomes(self) -> "dict[str, GenomeOption]":
        """Return the selectable genomes, ordered (human first, then mouse).

        An entry is included whenever its FASTA path is non-empty. The
        ``human`` entry is built from the legacy ``genome_fasta`` /
        ``gtf_path`` / ``genome_label`` fields so existing code and tests that
        set those directly keep working.
        """
        out: "dict[str, GenomeOption]" = {}
        if self.genome_fasta:
            out["human"] = GenomeOption(
                key="human",
                label=self.genome_label or "GRCh38",
                genome_fasta=self.genome_fasta,
                gtf_path=self.gtf_path,
                gtf_tag="MANE_Select",
                species="Homo sapiens",
            )
        if self.mouse_genome_fasta:
            out["mouse"] = GenomeOption(
                key="mouse",
                label=self.mouse_genome_label or "GRCm39",
                genome_fasta=self.mouse_genome_fasta,
                gtf_path=self.mouse_gtf_path,
                gtf_tag="Ensembl_canonical",
                species="Mus musculus",
            )
        return out

    def genome_option(self, key: "str | None" = None) -> "GenomeOption | None":
        """Resolve ``key`` to a :class:`GenomeOption`.

        With no (or an unknown) key, the configured ``default_genome`` is
        used, falling back to the first configured genome.
        """
        genomes = self.genomes
        if key and key in genomes:
            return genomes[key]
        if self.default_genome in genomes:
            return genomes[self.default_genome]
        return next(iter(genomes.values()), None)

    @property
    def beguidedoff_index_resolved(self) -> str:
        """FM-index path for fast mode.

        Explicit ``beguidedoff_index`` wins; otherwise the ``.npz`` sibling
        path (same basename, ``.fmi`` suffix) is used. Empty when no genome
        is configured at all.
        """
        if self.beguidedoff_index:
            return self.beguidedoff_index
        if self.beguidedoff_genome:
            return str(Path(self.beguidedoff_genome).with_suffix(".fmi"))
        return ""

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
