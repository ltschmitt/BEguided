"""Unit tests for ``beguided.web.genome_cache``."""

from __future__ import annotations

import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest


@pytest.fixture(autouse=True)
def _reset_module_state():
    """Ensure module state is clean for every test."""
    from beguided.web import genome_cache

    genome_cache._reset_for_test()
    yield
    genome_cache._reset_for_test()


def _write_fasta(path: Path, chroms: dict[str, str]) -> None:
    content = "".join(f">{name}\n{seq}\n" for name, seq in chroms.items())
    path.write_text(content)


def _wait_until(key: str, status: str, timeout: float = 2.0) -> dict:
    from beguided.web import genome_cache

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        snap = genome_cache.get_status(key)
        if snap["status"] == status:
            return snap
        time.sleep(0.02)
    raise AssertionError(f"key {key!r} never reached status {status!r}")


def test_preload_async_sets_ready_state(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT", "chr2": "TTTT"})

    genome_cache.preload_async("human", str(fasta))
    status = _wait_until("human", "ready")

    assert status["chromosome_count"] == 2
    assert status["path"] == str(fasta)
    # No label supplied and the FASTA filename doesn't match the Ensembl
    # pattern, so derive_genome_label falls back to "GRCh38".
    assert status["label"] == "GRCh38"

    genome = genome_cache.wait_until_ready("human", timeout=1.0)
    assert genome == {"chr1": "ACGT", "chr2": "TTTT"}


def test_preload_async_uses_explicit_label(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    genome_cache.preload_async("human", str(fasta), label="GRCh38.112")
    status = _wait_until("human", "ready")

    assert status["label"] == "GRCh38.112"


def test_preload_async_is_idempotent(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    genome_cache.preload_async("human", str(fasta))
    # Call a second time with the same path; must not raise.
    genome_cache.preload_async("human", str(fasta))
    _wait_until("human", "ready")
    assert genome_cache.is_ready("human")


def test_preload_async_rejects_different_path(tmp_path):
    from beguided.web import genome_cache

    a = tmp_path / "a.fa"
    b = tmp_path / "b.fa"
    _write_fasta(a, {"chr1": "AAAA"})
    _write_fasta(b, {"chr1": "CCCC"})

    genome_cache.preload_async("human", str(a))
    with pytest.raises(RuntimeError, match="already loaded"):
        genome_cache.preload_async("human", str(b))


def test_multiple_keys_have_independent_state(tmp_path):
    from beguided.web import genome_cache

    human = tmp_path / "human.fa"
    mouse = tmp_path / "mouse.fa"
    _write_fasta(human, {"1": "ACGT"})
    _write_fasta(mouse, {"1": "ACGT", "2": "TTTT"})

    genome_cache.preload_async("human", str(human), label="GRCh38")
    genome_cache.preload_async("mouse", str(mouse), label="GRCm39")

    _wait_until("human", "ready")
    _wait_until("mouse", "ready")

    aggregate = genome_cache.get_status()
    assert aggregate["status"] == "ready"
    assert aggregate["chromosome_count"] == 3
    assert {g["key"] for g in aggregate["genomes"]} == {"human", "mouse"}

    # With several registered genomes an unkeyed wait is ambiguous.
    with pytest.raises(RuntimeError, match="requires a key"):
        genome_cache.wait_until_ready(timeout=0.1)


def test_fail_soft_marks_error_without_exiting(tmp_path):
    """fail_fast=False (non-default genome) marks the key as error."""
    from beguided.web import genome_cache

    fasta = tmp_path / "broken.fa"
    _write_fasta(fasta, {"1": "ACGT"})

    with patch(
        "beguided.web.genome_cache.load_genome",
        side_effect=RuntimeError("boom"),
    ):
        genome_cache.preload_async("mouse", str(fasta), fail_fast=False)
        status = _wait_until("mouse", "error")

    assert "boom" in status["error"]
    with pytest.raises(RuntimeError, match="boom"):
        genome_cache.wait_until_ready("mouse", timeout=0.1)


def test_wait_until_ready_raises_after_error(tmp_path):
    """Simulate a load failure and verify waiters surface the error."""
    from beguided.web import genome_cache

    fasta = tmp_path / "broken.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    entry = genome_cache._get_or_create("human")
    with entry.cond:
        entry.status = "loading"
        entry.path = str(fasta)

    def _simulate_failure():
        with entry.cond:
            entry.status = "error"
            entry.error = "RuntimeError: boom"
            entry.finished_at = time.time()
            entry.cond.notify_all()

    results = {}

    def waiter():
        try:
            results["genome"] = genome_cache.wait_until_ready("human", timeout=2.0)
            results["ok"] = True
        except Exception as exc:
            results["ok"] = False
            results["error"] = exc

    t = threading.Thread(target=waiter, daemon=True)
    t.start()
    time.sleep(0.05)
    _simulate_failure()
    t.join(timeout=2.0)

    assert results.get("ok") is False, results
    assert "boom" in str(results["error"])
    assert genome_cache.get_status("human")["status"] == "error"


def test_wait_until_ready_raises_when_not_started():
    from beguided.web import genome_cache

    with pytest.raises(RuntimeError, match="not been started"):
        genome_cache.wait_until_ready(timeout=0.1)
    with pytest.raises(RuntimeError, match="not been preloaded"):
        genome_cache.wait_until_ready("human", timeout=0.1)


def test_wait_until_ready_times_out(tmp_path):
    """If preload never finishes within timeout, raise TimeoutError."""
    from beguided.web import genome_cache

    fasta = tmp_path / "slow.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    sentinel = threading.Event()

    def slow_load(_path):
        sentinel.wait(timeout=5.0)
        return {"chr1": "ACGT"}

    with patch("beguided.web.genome_cache.load_genome", side_effect=slow_load):
        genome_cache.preload_async("human", str(fasta))
        with pytest.raises(TimeoutError):
            genome_cache.wait_until_ready("human", timeout=0.2)

        # Release the loader so the daemon thread can exit cleanly.
        sentinel.set()


def test_wait_blocks_then_returns(tmp_path):
    """A caller waiting on the condition wakes up once preload finishes."""
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    barrier = threading.Event()
    started = threading.Event()
    results: dict = {}

    def waiter():
        try:
            started.set()
            results["genome"] = genome_cache.wait_until_ready("human", timeout=5.0)
            results["ok"] = True
        except Exception as exc:  # pragma: no cover - test failure path
            results["ok"] = False
            results["error"] = exc
        finally:
            barrier.set()

    # Put the key into "loading" state first so the waiter parks on the
    # condition, then drive the real load synchronously.
    entry = genome_cache._get_or_create("human")
    with entry.cond:
        entry.status = "loading"
        entry.path = str(fasta)

    t = threading.Thread(target=waiter, daemon=True)
    t.start()
    assert started.wait(timeout=2.0), "waiter never started"
    # Give the waiter time to actually call wait().
    time.sleep(0.1)

    genome_cache._load_blocking(entry, False)

    assert barrier.wait(timeout=3.0), "waiter never woke up"
    assert results.get("ok") is True, results
    assert results["genome"] == {"chr1": "ACGT"}
