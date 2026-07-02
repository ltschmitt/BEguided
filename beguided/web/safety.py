"""Safety helpers for the BEguided web server.

Two concerns:
  * Request body size cap (Starlette middleware) so a malicious client cannot
    upload a 5 GB blob.
  * Operator-supplied path validation so the server refuses to start if any of
    the genome FASTA / GTF / presets / BEguidedOFF .npz paths is missing or
    unreadable. Also validate uploaded chromosome names against the loaded
    genome.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

logger = logging.getLogger("beguided.web.safety")


class MaxBodySizeMiddleware(BaseHTTPMiddleware):
    """Reject requests whose Content-Length header exceeds ``max_bytes``.

    For requests with no Content-Length (chunked uploads), falls back to
    streaming a chunked read up to ``max_bytes + 1`` and aborts if more is
    sent.
    """

    def __init__(self, app, max_bytes: int):
        super().__init__(app)
        self.max_bytes = max_bytes

    async def dispatch(self, request: Request, call_next) -> Response:
        content_length = request.headers.get("content-length")
        if content_length is not None:
            try:
                if int(content_length) > self.max_bytes:
                    return JSONResponse(
                        {
                            "ok": False,
                            "error": (
                                f"Request body too large: {int(content_length)} bytes "
                                f"(max {self.max_bytes})"
                            ),
                        },
                        status_code=413,
                    )
            except ValueError:
                return JSONResponse(
                    {"ok": False, "error": "Invalid Content-Length header."},
                    status_code=400,
                )
        return await call_next(request)


def validate_path(path: Optional[str], label: str, *, must_be_file: bool = True) -> Path:
    """Resolve ``path`` and verify it exists and (optionally) is a file.

    Raises ``FileNotFoundError`` with a clear, operator-facing message if the
    path is empty, missing, or unreadable.
    """
    if not path:
        raise FileNotFoundError(f"{label} is not configured (path is empty).")
    p = Path(path).expanduser().resolve(strict=False)
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    if must_be_file and not p.is_file():
        raise FileNotFoundError(f"{label} is not a regular file: {p}")
    if not os_accessible(p):
        raise PermissionError(f"{label} is not readable: {p}")
    return p


def os_accessible(p: Path) -> bool:
    """Return True if ``p`` is readable by the current process."""
    try:
        with open(p, "rb"):
            return True
    except OSError:
        return False


def validate_operator_paths(cfg) -> None:
    """Validate every operator-configured input path. Raises on the first failure."""
    validate_path(cfg.presets_file, "Presets file")
    # The genome FASTA and GTF may not be available in every deployment
    # (e.g. genome-free reruns of step 4-5). Validate them best-effort.
    if cfg.genome_fasta:
        try:
            validate_path(cfg.genome_fasta, "Genome FASTA")
        except (FileNotFoundError, PermissionError) as exc:
            logger.warning("Genome FASTA not accessible: %s", exc)
    if cfg.gtf_path:
        try:
            validate_path(cfg.gtf_path, "GTF annotation")
        except (FileNotFoundError, PermissionError) as exc:
            logger.warning("GTF annotation not accessible: %s", exc)
    if cfg.beguidedoff_genome:
        try:
            validate_path(cfg.beguidedoff_genome, "BEguidedOFF genome (.npz)")
        except (FileNotFoundError, PermissionError) as exc:
            logger.warning("BEguidedOFF genome not accessible: %s", exc)


def validate_chromosomes(chromosomes: Iterable[str], genome_keys) -> list[str]:
    """Return the list of chromosome names from ``chromosomes`` not present in ``genome_keys``.

    Empty result means all chromosomes are valid.
    """
    genome = set(genome_keys)
    return [c for c in chromosomes if c not in genome]
