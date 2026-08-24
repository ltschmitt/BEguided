"""Legal-page content loaded from an external YAML file.

The legal pages (Imprint / Privacy / Terms of Use) are required for any
public-facing German web presence (TMG/DDG §5, GDPR Art. 13/14).
Because they contain personal data (names, addresses, DPO contacts) they
must not be committed to the public repository.

The operator mounts a ``legal.yaml`` file into the container at a path
of their choice (default ``/etc/beguided/legal.yaml``) and the
contents are loaded once at startup. If the file is missing, every
field falls back to an empty string and a warning is logged — the
server still starts so local development and CI tests work without
any configuration. If the file is present but malformed (invalid
YAML, wrong top-level type), startup fails fast because serving a
half-parsed legal page is worse than refusing to start.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml

logger = logging.getLogger("beguided.web.legal_config")


@dataclass
class LegalConfig:
    """Legal-page content. Every field is a string (or list of strings).

    All fields default to empty so a missing ``legal.yaml`` still produces
    a 200 response — the pages just render with blanks, and the
    startup-time warning is the operator's signal to fix it.
    """

    institution_name: str = ""
    address_lines: List[str] = field(default_factory=list)
    vertretungsberechtigter: str = ""
    contact_email: str = ""
    contact_phone: str = ""
    hosting_provider: str = ""
    content_responsible_mstv: str = ""
    dpo_name: str = ""
    dpo_email: str = ""
    dpo_address: List[str] = field(default_factory=list)
    competent_dpa: str = ""
    log_retention_days: int = 0
    last_updated: str = ""

    # Higher-level institution reference (template uses
    # "Technische Universität Dresden" throughout the imprint).
    hosting_legal_name: str = ""
    hosting_legal_postal: List[str] = field(default_factory=list)
    hosting_legal_email: str = ""
    # Head of institution for "Technical implementation" (template maps
    # this to vertretungsberechtigter when empty).
    head_of_institution: str = ""
    # Public URL of the service (replaces the template's
    # <Please enter your URL here> placeholder).
    website_url: str = ""
    # Data Protection Officer phone/fax (template defaults below).
    dpo_phone: str = ""
    dpo_fax: str = ""
    # Competent supervisory authority head/address/email/phone/web
    # (template defaults below for the Saxon DPA).
    competent_dpa_head: str = ""
    competent_dpa_address: List[str] = field(default_factory=list)
    competent_dpa_email: str = ""
    competent_dpa_phone: str = ""
    competent_dpa_web: str = ""
    # Lead phrase: "The Legal Notice of <legal_notice_host> applies with
    # the following amendments:". When `legal_notice_host_url` is set,
    # the host name is rendered as a hyperlink to that URL; otherwise
    # it is plain text. No fallback when either field is empty.
    legal_notice_host: str = ""
    legal_notice_host_url: str = ""
    # Name of the person to contact for content questions (rendered in
    # the Responsibilities section of the imprint).
    content_contact_name: str = ""
    # Department / institution subdivision (rendered under Technical
    # implementation).
    institution_department: str = ""


# Fields the loader accepts. Anything else in the YAML is ignored so the
# file can carry comments or future fields without breaking older builds.
_KNOWN_SCALAR_FIELDS = (
    "institution_name",
    "vertretungsberechtigter",
    "contact_email",
    "contact_phone",
    "hosting_provider",
    "content_responsible_mstv",
    "dpo_name",
    "dpo_email",
    "competent_dpa",
    "log_retention_days",
    "last_updated",
    "hosting_legal_name",
    "hosting_legal_email",
    "head_of_institution",
    "website_url",
    "dpo_phone",
    "dpo_fax",
    "competent_dpa_head",
    "competent_dpa_email",
    "competent_dpa_phone",
    "competent_dpa_web",
    "legal_notice_host",
    "legal_notice_host_url",
    "content_contact_name",
    "institution_department",
)
_KNOWN_LIST_FIELDS = (
    "address_lines",
    "dpo_address",
    "hosting_legal_postal",
    "competent_dpa_address",
)


def load_legal_config(path: Optional[Path]) -> LegalConfig:
    """Load legal content from ``path`` or return an empty config.

    - ``path is None`` → return ``LegalConfig()`` silently. Used when the
      operator has not configured ``BEGUIDED_LEGAL_CONFIG`` at all (tests,
      local dev with no legal pages).
    - ``path`` points at a missing file OR a directory → return
      ``LegalConfig()`` and log a WARNING. The operator probably forgot
      to mount ``legal.yaml``, or the bind mount target became a directory
      because the host file did not exist when docker compose started.
    - ``path`` exists and parses cleanly → return the populated config
      and log an INFO line so the startup banner confirms the source.
    - ``path`` exists but YAML is malformed or the top-level is not a
      mapping → raise ``ValueError`` so ``create_app()`` aborts startup.
    """
    if path is None or str(path).strip() == "":
        logger.info("No legal config path configured; serving empty legal pages.")
        return LegalConfig()

    p = Path(path).expanduser()
    if not p.exists() or p.is_dir():
        logger.warning(
            "Legal config %s not found or is not a file; serving empty "
            "legal pages. Create the file or update BEGUIDED_LEGAL_CONFIG.",
            p,
        )
        return LegalConfig()

    try:
        raw = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ValueError(
            f"Failed to parse legal config {p}: {exc}. "
            "Fix the YAML or unset BEGUIDED_LEGAL_CONFIG to fall back to "
            "empty values."
        ) from exc

    if raw is None:
        logger.info("Legal config %s is empty; serving empty legal pages.", p)
        return LegalConfig()
    if not isinstance(raw, dict):
        raise ValueError(
            f"Legal config {p} must contain a YAML mapping at the top "
            f"level, got {type(raw).__name__}."
        )

    cfg = LegalConfig()
    for key in _KNOWN_SCALAR_FIELDS:
        if key in raw:
            value = raw[key]
            if key == "log_retention_days":
                try:
                    cfg.log_retention_days = int(value)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Legal config {p}: log_retention_days must be an "
                        f"integer, got {value!r}."
                    ) from exc
            else:
                if value is None:
                    cfg.__setattr__(key, "")
                elif not isinstance(value, str):
                    raise ValueError(
                        f"Legal config {p}: {key} must be a string, "
                        f"got {type(value).__name__}."
                    )
                else:
                    cfg.__setattr__(key, value)
    for key in _KNOWN_LIST_FIELDS:
        if key in raw:
            value = raw[key]
            if value is None:
                cfg.__setattr__(key, [])
            elif isinstance(value, str):
                # Allow a single string for convenience: it becomes a
                # one-element list so the templates can iterate uniformly.
                cfg.__setattr__(key, [value])
            elif isinstance(value, list):
                if not all(isinstance(item, str) for item in value):
                    raise ValueError(
                        f"Legal config {p}: {key} must be a list of "
                        f"strings."
                    )
                cfg.__setattr__(key, list(value))
            else:
                raise ValueError(
                    f"Legal config {p}: {key} must be a string or list "
                    f"of strings, got {type(value).__name__}."
                )

    logger.info(
        "Loaded legal config from %s (institution=%r, dpo=%r).",
        p,
        cfg.institution_name or "<unset>",
        cfg.dpo_name or "<unset>",
    )
    return cfg
