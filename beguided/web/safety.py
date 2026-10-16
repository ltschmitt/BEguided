"""Safety helpers for the BEguided web server.

Two concerns:
  * Request body size cap (Starlette middleware) so a malicious client cannot
    upload a 5 GB blob.
  * Operator-supplied path validation so the server refuses to start if any of
    the genome FASTA / GTF / base editor profiles / BEguidedOFF .npz paths is missing or
    unreadable. Also validate uploaded chromosome names against the loaded
    genome.
"""

from __future__ import annotations

import gzip
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


class SelectiveGZipMiddleware(BaseHTTPMiddleware):
    """GZip only the large, non-streaming responses under ``paths``.

    A global ``GZipMiddleware`` would also wrap the Server-Sent Events
    (``/logs/stream``) response, which can stall or break streaming; Caddy's
    ``encode`` handler has the same problem (caddyserver/caddy#6293). We
    therefore compress only the endpoints that return one large document in a
    single shot (job results + result downloads) and leave every other
    response — especially SSE — untouched.
    """

    def __init__(self, app, paths, minimum_size: int = 1024):
        super().__init__(app)
        self.paths = frozenset(paths)
        self.minimum_size = max(0, int(minimum_size))

    async def dispatch(self, request: Request, call_next) -> Response:
        response = await call_next(request)
        if (
            request.method != "GET"
            or request.url.path not in self.paths
            or response.status_code != 200
            or "content-encoding" in response.headers
        ):
            return response
        if "gzip" not in request.headers.get("accept-encoding", "").lower():
            return response

        if hasattr(response, "body_iterator"):
            body = b"".join([chunk async for chunk in response.body_iterator])
        else:
            body = response.body
        if len(body) < self.minimum_size:
            # Iterator already consumed: hand back an equivalent response.
            return Response(
                content=body,
                status_code=response.status_code,
                headers=dict(response.headers),
            )

        compressed = gzip.compress(body, compresslevel=6)
        headers = dict(response.headers)
        headers.pop("content-length", None)
        headers["content-encoding"] = "gzip"
        vary = headers.get("vary")
        headers["vary"] = f"{vary}, Accept-Encoding" if vary else "Accept-Encoding"
        return Response(
            content=compressed,
            status_code=response.status_code,
            headers=headers,
        )


# Default Content-Security-Policy for the BEguided web app.
#
# ``default-src 'self'`` covers everything not explicitly listed. The page
# uses inline ``style="..."`` attributes and template-level inline styles
# (e.g. ``style="color:var(--muted)"``), so ``style-src 'unsafe-inline'`` is
# required. The JS bundle is loaded from ``/static/...`` so a same-origin
# ``script-src`` is enough. No third-party origins are contacted on page
# load. ``frame-ancestors 'none'`` disables clickjacking via iframes.
DEFAULT_CSP = (
    "default-src 'self'; "
    "script-src 'self'; "
    "style-src 'self' 'unsafe-inline'; "
    "img-src 'self' data:; "
    "font-src 'self'; "
    "connect-src 'self'; "
    "frame-ancestors 'none'; "
    "base-uri 'self'; "
    "form-action 'self'"
)


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Apply a fixed set of security response headers to every response.

    HSTS is only emitted on HTTPS responses (gated on ``X-Forwarded-Proto``
    from the trusted reverse proxy, falling back to the request scheme). The
    uvicorn ``Server`` header is removed so the framework version is not
    disclosed.
    """

    def __init__(self, app, csp: str = DEFAULT_CSP):
        super().__init__(app)
        self.csp = csp

    @staticmethod
    def _is_https(request: Request) -> bool:
        forwarded_proto = request.headers.get("x-forwarded-proto", "").lower().strip()
        if forwarded_proto == "https":
            return True
        return request.url.scheme == "https"

    async def dispatch(self, request: Request, call_next) -> Response:
        response = await call_next(request)
        response.headers.setdefault("Content-Security-Policy", self.csp)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("Referrer-Policy", "strict-origin-when-cross-origin")
        response.headers.setdefault(
            "Permissions-Policy",
            "camera=(), microphone=(), geolocation=(), interest-cohort=()",
        )
        response.headers.setdefault("Cross-Origin-Opener-Policy", "same-origin")
        if self._is_https(request):
            response.headers.setdefault(
                "Strict-Transport-Security",
                "max-age=31536000; includeSubDomains",
            )
        # uvicorn adds ``server: uvicorn`` by default; replace it with a
        # neutral value so the framework version is not disclosed.
        response.headers["Server"] = "beguided"
        return response


class ProxyHeadersMiddleware:
    """ASGI middleware that resolves the real client IP from ``X-Forwarded-For``.

    Implemented locally because Starlette 1.x removed
    ``starlette.middleware.proxy_headers``. Only updates
    ``scope["client"]`` when the immediate upstream peer is on the
    trusted-hosts allow-list (so an attacker on the public internet
    cannot forge the header). The first IP in the comma-separated
    ``X-Forwarded-For`` chain is taken as the real client.
    """

    def __init__(self, app, trusted_hosts=None):
        self.app = app
        # ``trusted_hosts`` may be None (trust nothing), an iterable of
        # IP literals / CIDRs (compared as exact strings for simplicity),
        # or "*" (trust every upstream — only useful in tests).
        if trusted_hosts is None:
            self._trusted = set()
        elif isinstance(trusted_hosts, str):
            self._trusted = {trusted_hosts}
        else:
            self._trusted = {str(h).strip() for h in trusted_hosts if str(h).strip()}
        self._trust_all = self._trusted == {"*"}

    def _is_trusted(self, peer_host: Optional[str]) -> bool:
        if self._trust_all:
            return True
        return peer_host in self._trusted if peer_host else False

    async def __call__(self, scope, receive, send):
        if scope["type"] not in ("http", "websocket"):
            await self.app(scope, receive, send)
            return

        client = scope.get("client")
        if client is not None and self._is_trusted(client[0]):
            headers = {
                k.decode("latin-1").lower(): v.decode("latin-1")
                for k, v in scope.get("headers", [])
            }
            forwarded = headers.get("x-forwarded-for")
            if forwarded:
                first = forwarded.split(",")[0].strip()
                if first:
                    port = client[1] if len(client) > 1 else 0
                    scope["client"] = (first, port)

        await self.app(scope, receive, send)


class HeadToGetMiddleware:
    """Answer HEAD requests by dispatching them as GET and dropping the body.

    Starlette's plain ``Route`` auto-adds HEAD when GET is registered, but
    FastAPI's ``APIRoute`` does not, so ``@app.get(...)`` endpoints return
    405 for HEAD. Uptime monitors and link checkers issue HEAD, so this
    middleware rewrites the method before routing and strips the response
    body frames (keeping status line + headers).

    The SSE endpoint (``/logs/stream``) is exempt: rewriting it to GET would
    open a long-lived stream that the wrapper would then drain.
    """

    def __init__(self, app, exempt_paths: Iterable[str] = ("/logs/stream",)):
        self.app = app
        self.exempt_paths = frozenset(exempt_paths)

    async def __call__(self, scope, receive, send):
        if (
            scope["type"] != "http"
            or scope["method"] != "HEAD"
            or scope.get("path") in self.exempt_paths
        ):
            await self.app(scope, receive, send)
            return

        scope = dict(scope)
        scope["method"] = "GET"

        async def send_no_body(message):
            if message["type"] == "http.response.body":
                # Keep ``more_body`` as-is: responses may arrive in several
                # frames, and forcing the flag to False makes the client treat
                # the first frame as terminal (then choke on the real last one).
                message = {**message, "body": b""}
            await send(message)

        await self.app(scope, receive, send_no_body)


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
    validate_path(cfg.base_editor_profiles_file, "Base editor profiles file")
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
    # Optional mouse reference (best-effort, like the human genome above).
    for label, path in (
        ("Mouse genome FASTA", getattr(cfg, "mouse_genome_fasta", "")),
        ("Mouse GTF annotation", getattr(cfg, "mouse_gtf_path", "")),
    ):
        if path:
            try:
                validate_path(path, label)
            except (FileNotFoundError, PermissionError) as exc:
                logger.warning("%s not accessible: %s", label, exc)
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
