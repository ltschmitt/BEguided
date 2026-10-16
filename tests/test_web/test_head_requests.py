"""Tests for the HEAD → GET middleware.

FastAPI's ``APIRoute`` does not auto-add HEAD to ``@app.get`` routes (unlike
Starlette's plain ``Route``), so ``curl -I`` / uptime monitors used to get 405.
``HeadToGetMiddleware`` rewrites HEAD to GET and strips the body frames.
"""

from __future__ import annotations

import asyncio

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("starlette")

from beguided.web.safety import HeadToGetMiddleware  # noqa: E402


class _RecordingApp:
    """Minimal ASGI app: 200 + body for GET, 405 for anything else."""

    def __init__(self) -> None:
        self.methods: list[str] = []

    async def __call__(self, scope, receive, send):
        self.methods.append(scope["method"])
        if scope["method"] == "GET":
            await send(
                {
                    "type": "http.response.start",
                    "status": 200,
                    "headers": [(b"content-type", b"text/plain")],
                }
            )
            # Two frames, to exercise multi-frame body handling.
            await send(
                {"type": "http.response.body", "body": b"he", "more_body": True}
            )
            await send({"type": "http.response.body", "body": b"llo"})
        else:
            await send({"type": "http.response.start", "status": 405, "headers": []})
            await send({"type": "http.response.body", "body": b""})


def _run(mw, method: str, path: str = "/"):
    sent: list[dict] = []

    async def send(message):
        sent.append(message)

    async def receive():
        return {"type": "http.request"}

    scope = {"type": "http", "method": method, "path": path, "headers": []}
    asyncio.run(mw(scope, receive, send))
    return sent


def test_head_is_rewritten_to_get_with_empty_body():
    inner = _RecordingApp()
    sent = _run(HeadToGetMiddleware(inner), "HEAD", "/")

    assert inner.methods == ["GET"], "inner app must see GET"
    start = next(m for m in sent if m["type"] == "http.response.start")
    assert start["status"] == 200
    assert (b"content-type", b"text/plain") in start["headers"], "headers preserved"
    body = [m for m in sent if m["type"] == "http.response.body"]
    assert body and all(m["body"] == b"" for m in body), "body frames stripped"
    # Framing must be preserved (first frame more_body=True, last False).
    assert body[0]["more_body"] is True
    assert body[-1].get("more_body", False) is False


def test_get_is_untouched():
    inner = _RecordingApp()
    sent = _run(HeadToGetMiddleware(inner), "GET", "/")
    assert inner.methods == ["GET"]
    bodies = b"".join(m["body"] for m in sent if m["type"] == "http.response.body")
    assert bodies == b"hello"


def test_sse_path_is_exempt_from_rewrite():
    """A HEAD on the SSE endpoint must not open a drained GET stream."""
    inner = _RecordingApp()
    _run(HeadToGetMiddleware(inner), "HEAD", "/logs/stream")
    assert inner.methods == ["HEAD"]


def test_head_through_real_fastapi_app():
    """End-to-end: an ``@app.get`` route answers HEAD after the middleware."""
    from fastapi import FastAPI
    from fastapi.responses import HTMLResponse
    from fastapi.testclient import TestClient

    def make_app(with_middleware: bool) -> FastAPI:
        app = FastAPI()

        @app.get("/", response_class=HTMLResponse)
        def index() -> str:
            return "<html>hi</html>"

        if with_middleware:
            app.add_middleware(HeadToGetMiddleware)
        return app

    # Without the middleware FastAPI returns 405 for HEAD.
    assert TestClient(make_app(False)).head("/").status_code == 405

    client = TestClient(make_app(True))
    resp = client.head("/")
    assert resp.status_code == 200
    assert resp.content == b""
    assert client.get("/").content == b"<html>hi</html>"
