"""Legal-page content loaded from an external YAML file.

The legal pages (Imprint / Privacy / Terms of Use / Accessibility) are
required for any public-facing German web presence (TMG/DDG §5, GDPR
Art. 13/14, BfWebG/BITV 2.0). Because they contain personal data
(names, addresses, DPO contacts) they must not be committed to the
public repository.

The operator mounts a ``legal.yaml`` file into the container at a path
of their choice (default ``/etc/beguided/legal.yaml``) and the
contents are loaded once at startup. If the file is missing, every
field falls back to an empty value and a warning is logged — the
server still starts so local development and CI tests work without
any configuration. If the file is present but malformed (invalid
YAML, wrong top-level type, unknown variable reference) startup fails
fast because serving a half-parsed legal page is worse than refusing
to start.

Variable reuse
--------------

Every string value is rendered as a Jinja template against the *other
values in the file*, so a fact defined once (e.g. ``institution_name``)
can be placed in many paragraphs:

.. code-block:: yaml

   institution_name: "UCC, Medical Systems Biology, TU Dresden"
   controller_block:
     - "{{ institution_name }}"
     - "Fetscherstr. 74"

The file is resolved to a fixed point (iterating a few passes), so
fields may reference each other in either order. A reference to a
variable that is not defined anywhere in the file raises ``ValueError``
at startup (fail-fast on typos).

Blocks
------

Full-paragraph "block" fields (``controller_block``, ``dpo_block``,
``dpa_block``, ``access_contact_block``, ``enforcement_body_block``)
accept a single string or a list of strings; either way the value is
rendered end-to-end (so Jinja ``{% for %}`` loops over other list
fields work) and split into lines. The templates render the lines
escaped, and tokens starting with ``mailto:`` / ``http(s)://`` are
turned into hyperlinks (mailto links show the address without the
``mailto:`` prefix).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import jinja2
import yaml
from jinja2.sandbox import SandboxedEnvironment

logger = logging.getLogger("beguided.web.legal_config")

# How many fixed-point passes the variable resolver makes before giving
# up. One pass resolves every reference whose target is defined earlier
# in the file; a second pass resolves the remainder. Anything that still
# changes after five passes is treated as non-convergent (circular).
_MAX_RESOLVE_PASSES = 5

_TEMPLATE_MARKERS = ("{{", "{%", "{#")


@dataclass
class LegalConfig:
    """Legal-page content.

    Base facts (institution, address, contact) are defined once and
    reused by the templates and by the block paragraphs; the block
    fields carry the full multi-line paragraphs for the privacy policy
    and the accessibility statement. All fields default to empty so a
    missing ``legal.yaml`` still produces a 200 response — the pages
    just render with blanks, and the startup-time warning is the
    operator's signal to fix it.
    """

    # Base facts — single source of truth, referenced elsewhere.
    institution_name: str = ""
    address_lines: List[str] = field(default_factory=list)
    contact_name: str = ""
    contact_email: str = ""
    contact_phone: str = ""
    website_url: str = ""
    head_of_institution: str = ""
    # Person responsible for the content of the service (imprint, per
    # § 18 MStV); usually the PI / head of the institute.
    content_responsible: str = ""
    # Person who technically implemented the web service (imprint).
    technical_implementer: str = ""
    technical_implementer_email: str = ""
    hosting_legal_name: str = ""
    institution_department: str = ""
    # Lead phrase: "The Legal Notice of <legal_notice_host> applies with
    # the following amendments:". When `legal_notice_host_url` is set,
    # the host name is rendered as a hyperlink to that URL; otherwise
    # it is plain text. No fallback when either field is empty.
    legal_notice_host: str = ""
    legal_notice_host_url: str = ""
    # Free-form dates.
    last_updated: str = ""
    accessibility_prepared: str = ""
    accessibility_reviewed: str = ""
    # How long reverse-proxy / web-server logs are kept (separate from
    # the job-database retention enforced by the app itself).
    log_retention_days: int = 0

    # Full-paragraph blocks (lists of rendered lines).
    controller_block: List[str] = field(default_factory=list)
    dpo_block: List[str] = field(default_factory=list)
    dpa_block: List[str] = field(default_factory=list)
    access_contact_block: List[str] = field(default_factory=list)
    enforcement_body_block: List[str] = field(default_factory=list)


# String scalars the loader accepts. Anything else in the YAML is
# ignored so the file can carry comments or future fields without
# breaking older builds.
_KNOWN_SCALAR_FIELDS = (
    "institution_name",
    "contact_name",
    "contact_email",
    "contact_phone",
    "website_url",
    "head_of_institution",
    "content_responsible",
    "technical_implementer",
    "technical_implementer_email",
    "hosting_legal_name",
    "institution_department",
    "legal_notice_host",
    "legal_notice_host_url",
    "last_updated",
    "accessibility_prepared",
    "accessibility_reviewed",
)
# Block fields: a string or a list of strings, resolved as one unit and
# stored as a list of lines.
_KNOWN_BLOCK_FIELDS = (
    "controller_block",
    "dpo_block",
    "dpa_block",
    "access_contact_block",
    "enforcement_body_block",
)
_KNOWN_LIST_FIELDS = (
    "address_lines",
)
# Fields removed by the block-based schema. Still listed here so a
# half-migrated file logs a loud warning instead of silently serving
# blank sections.
_LEGACY_FIELDS = (
    "vertretungsberechtigter",
    "content_responsible_mstv",
    "hosting_provider",
    "content_contact_name",
    "dpo_name",
    "dpo_email",
    "dpo_address",
    "dpo_phone",
    "dpo_fax",
    "competent_dpa",
    "competent_dpa_head",
    "competent_dpa_address",
    "competent_dpa_email",
    "competent_dpa_phone",
    "competent_dpa_web",
    "hosting_legal_postal",
    "hosting_legal_email",
)


def _has_template(text: str) -> bool:
    return any(marker in text for marker in _TEMPLATE_MARKERS)


def _resolve_variables(
    values: Dict[str, object], p: Path
) -> None:
    """Render every string value as a Jinja template against the other
    values, iterating to a fixed point.

    ``values`` is mutated in place. List fields stay lists (so they can
    be iterated from templates); their string items and the scalar
    fields are resolved independently. A reference to an undefined
    variable raises ``ValueError`` so a typo fails startup instead of
    rendering a blank legal paragraph.
    """
    env = SandboxedEnvironment(
        undefined=jinja2.StrictUndefined,
        autoescape=False,
        keep_trailing_newline=True,
    )
    for _ in range(_MAX_RESOLVE_PASSES):
        changed = False
        context = dict(values)
        for key, value in list(values.items()):
            if isinstance(value, str):
                if not _has_template(value):
                    continue
                rendered = _render(env, value, key, context, p)
                if rendered != value:
                    values[key] = rendered
                    changed = True
            elif isinstance(value, list):
                for idx, item in enumerate(value):
                    if not isinstance(item, str) or not _has_template(item):
                        continue
                    rendered = _render(env, item, f"{key}[{idx}]", context, p)
                    if rendered != item:
                        value[idx] = rendered
                        changed = True
        if not changed:
            return
    raise ValueError(
        f"Legal config {p}: variable references do not converge after "
        f"{_MAX_RESOLVE_PASSES} passes (circular reference?)."
    )


def _render(
    env: SandboxedEnvironment,
    text: str,
    key: str,
    context: Dict[str, object],
    p: Path,
) -> str:
    """Render one value, translating an unknown-variable reference into a
    startup error that names the field."""
    try:
        return env.from_string(text).render(context)
    except jinja2.exceptions.UndefinedError as exc:
        raise ValueError(
            f"Legal config {p}: {key} references an undefined variable "
            f"({exc}). Check the spelling against the other fields in "
            f"the file."
        ) from exc


def _to_lines(value: object) -> List[str]:
    """Collapse a rendered string (possibly multi-line) into clean lines."""
    return [line for line in (line.strip() for line in str(value).splitlines()) if line]


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
    - ``path`` exists but YAML is malformed, the top-level is not a
      mapping, a variable reference is unknown, or a field has the wrong
      type → raise ``ValueError`` so ``create_app()`` aborts startup.
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

    legacy = sorted(set(raw) & set(_LEGACY_FIELDS))
    if legacy:
        logger.warning(
            "Legal config %s contains legacy field(s) %s that are no "
            "longer used by the block-based schema; ignoring them. See "
            "configs/web-legal.yaml for the current fields.",
            p,
            ", ".join(legacy),
        )

    values: Dict[str, object] = {}

    # log_retention_days: the only non-string field, validated up front.
    if "log_retention_days" in raw:
        value = raw["log_retention_days"]
        if value is None:
            value = 0
        try:
            values["log_retention_days"] = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Legal config {p}: log_retention_days must be an "
                f"integer, got {value!r}."
            ) from exc

    for key in _KNOWN_SCALAR_FIELDS:
        if key not in raw:
            continue
        value = raw[key]
        if value is None:
            values[key] = ""
        elif not isinstance(value, str):
            raise ValueError(
                f"Legal config {p}: {key} must be a string, "
                f"got {type(value).__name__}."
            )
        else:
            values[key] = value

    for key in _KNOWN_BLOCK_FIELDS:
        if key not in raw:
            continue
        value = raw[key]
        if value is None:
            values[key] = ""
        elif isinstance(value, str):
            values[key] = value
        elif isinstance(value, list):
            if not all(isinstance(item, str) for item in value):
                raise ValueError(
                    f"Legal config {p}: {key} must be a string or a "
                    f"list of strings."
                )
            # Join before resolution so a single Jinja pass can run
            # loops over other list fields across the whole paragraph.
            values[key] = "\n".join(value)
        else:
            raise ValueError(
                f"Legal config {p}: {key} must be a string or list of "
                f"strings, got {type(value).__name__}."
            )

    for key in _KNOWN_LIST_FIELDS:
        if key not in raw:
            continue
        value = raw[key]
        if value is None:
            values[key] = []
        elif isinstance(value, str):
            # Allow a single string for convenience.
            values[key] = [value]
        elif isinstance(value, list):
            if not all(isinstance(item, str) for item in value):
                raise ValueError(
                    f"Legal config {p}: {key} must be a string or list "
                    f"of strings."
                )
            values[key] = list(value)
        else:
            raise ValueError(
                f"Legal config {p}: {key} must be a string or list "
                f"of strings, got {type(value).__name__}."
            )

    _resolve_variables(values, p)

    # A resolved list item may itself expand to multiple lines
    # (e.g. "{% for line in address_lines %}..."), so re-split.
    for key in _KNOWN_LIST_FIELDS:
        if key in values:
            values[key] = [
                line
                for item in values[key]  # type: ignore[str-index]
                for line in _to_lines(item)
            ]

    cfg = LegalConfig()
    for key in _KNOWN_SCALAR_FIELDS:
        if key in values:
            cfg.__setattr__(key, values[key])
    for key in _KNOWN_BLOCK_FIELDS:
        if key in values:
            cfg.__setattr__(key, _to_lines(values[key]))
    for key in _KNOWN_LIST_FIELDS:
        if key in values:
            cfg.__setattr__(key, values[key])
    if "log_retention_days" in values:
        cfg.log_retention_days = values["log_retention_days"]  # type: ignore[arg-type]

    logger.info(
        "Loaded legal config from %s (institution=%r, controller block=%s).",
        p,
        cfg.institution_name or "<unset>",
        "set" if cfg.controller_block else "<unset>",
    )
    return cfg
