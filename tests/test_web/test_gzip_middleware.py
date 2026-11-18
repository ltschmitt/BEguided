"""Tests for the selective GZip middleware.

Only the large single-shot result endpoints may be compressed; the SSE
(``text/event-stream``) endpoint must never be, or the live log stream stalls.
"""

from __future__ import annotations

import asyncio
import gzip

import pytest

# The default pixi environment has no starlette/fastapi (web feature only);
# skip cleanly there instead of erroring at collection time — same pattern as
# test_head_requests.py.
pytest.importorskip("starlette")

from starlette.requests import Request  # noqa: E402
from starlette.responses import Response  # noqa: E402

from beguided.web.safety import SelectiveGZipMiddleware  # noqa: E402

PATHS = ("/results", "/download/results")


def _request(path: str, accept_encoding: str | None = "gzip") -> Request:
    headers = []
    if accept_encoding:
        headers.append((b"accept-encoding", accept_encoding.encode()))
    scope = {
        "type": "http",
        "method": "GET",
        "path": path,
        "headers": headers,
        "query_string": b"",
        "scheme": "http",
        "server": ("test", 80),
        "client": ("1.2.3.4", 1234),
    }
    return Request(scope)


def _dispatch(path: str, body: bytes, accept_encoding: str | None = "gzip"):
    mw = SelectiveGZipMiddleware(lambda scope, receive, send: None, paths=PATHS)

    async def call_next(_request):
        return Response(body, media_type="application/json")

    return asyncio.run(mw.dispatch(_request(path, accept_encoding), call_next))


def test_large_result_is_gzipped():
    body = b"x" * 5000
    resp = _dispatch("/results", body)
    assert resp.headers["content-encoding"] == "gzip"
    assert "Accept-Encoding" in resp.headers["vary"]
    assert gzip.decompress(resp.body) == body


def test_small_result_is_not_gzipped():
    body = b"x" * 10
    resp = _dispatch("/results", body)
    assert "content-encoding" not in resp.headers
    assert resp.body == body


def test_non_result_path_is_not_gzipped():
    body = b"x" * 5000
    resp = _dispatch("/logs/stream", body)
    assert "content-encoding" not in resp.headers
    assert resp.body == body


def test_no_accept_encoding_is_not_gzipped():
    body = b"x" * 5000
    resp = _dispatch("/results", body, accept_encoding=None)
    assert "content-encoding" not in resp.headers
    assert resp.body == body
