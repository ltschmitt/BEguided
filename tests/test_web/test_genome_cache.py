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


def test_preload_async_sets_ready_state(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT", "chr2": "TTTT"})

    genome_cache.preload_async(str(fasta))
    # Wait up to 2s for the daemon thread to finish.
    deadline = time.monotonic() + 2.0
    while time.monotonic() < deadline:
        if genome_cache.is_ready():
            break
        time.sleep(0.02)

    status = genome_cache.get_status()
    assert status["status"] == "ready"
    assert status["chromosome_count"] == 2
    assert status["path"] == str(fasta)
    # No label supplied and the FASTA filename doesn't match the Ensembl
    # pattern, so derive_genome_label falls back to "GRCh38".
    assert status["label"] == "GRCh38"

    genome = genome_cache.wait_until_ready(timeout=1.0)
    assert genome == {"chr1": "ACGT", "chr2": "TTTT"}


def test_preload_async_uses_explicit_label(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    genome_cache.preload_async(str(fasta), label="GRCh38.112")
    deadline = time.monotonic() + 2.0
    while not genome_cache.is_ready():
        time.sleep(0.02)

    assert genome_cache.get_status()["label"] == "GRCh38.112"


def test_preload_async_is_idempotent(tmp_path):
    from beguided.web import genome_cache

    fasta = tmp_path / "tiny.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    genome_cache.preload_async(str(fasta))
    # Call a second time with the same path; must not raise.
    genome_cache.preload_async(str(fasta))
    deadline = time.monotonic() + 2.0
    while not genome_cache.is_ready():
        time.sleep(0.02)
    assert genome_cache.is_ready()


def test_preload_async_rejects_different_path(tmp_path):
    from beguided.web import genome_cache

    a = tmp_path / "a.fa"
    b = tmp_path / "b.fa"
    _write_fasta(a, {"chr1": "AAAA"})
    _write_fasta(b, {"chr1": "CCCC"})

    genome_cache.preload_async(str(a))
    with pytest.raises(RuntimeError, match="already loaded"):
        genome_cache.preload_async(str(b))
    genome_cache._reset_for_test()


def test_wait_until_ready_raises_after_error(tmp_path):
    """Simulate a load failure: ``load_genome`` raises and ``os._exit`` is
    intercepted so the test process keeps running."""
    from beguided.web import genome_cache

    fasta = tmp_path / "broken.fa"
    _write_fasta(fasta, {"chr1": "ACGT"})

    # Mark the preload as in-flight so ``wait_until_ready`` will block.
    # We then fail the load by raising from a synchronous helper rather than
    # driving the daemon thread (which calls os._exit and would kill pytest).
    genome_cache._state["status"] = "loading"
    genome_cache._state["path"] = str(fasta)

    with patch(
        "beguided.web.genome_cache.load_genome",
        side_effect=RuntimeError("boom"),
    ):
        # Simulate the daemon thread's failure path: it sets state, notifies
        # waiters, and exits the process. We replay that here synchronously.
        def _simulate_failure():
            with genome_cache._cond:
                genome_cache._state["status"] = "error"
                genome_cache._state["error"] = "RuntimeError: boom"
                genome_cache._state["finished_at"] = time.time()
                genome_cache._cond.notify_all()

        # Start the waiter first.
        results = {}

        def waiter():
            try:
                results["genome"] = genome_cache.wait_until_ready(timeout=2.0)
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
        assert genome_cache.get_status()["status"] == "error"


def test_wait_until_ready_raises_when_not_started():
    from beguided.web import genome_cache

    with pytest.raises(RuntimeError, match="not been started"):
        genome_cache.wait_until_ready(timeout=0.1)


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
        genome_cache.preload_async(str(fasta))
        with pytest.raises(TimeoutError):
            genome_cache.wait_until_ready(timeout=0.2)

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
            results["genome"] = genome_cache.wait_until_ready(timeout=5.0)
            results["ok"] = True
        except Exception as exc:  # pragma: no cover - test failure path
            results["ok"] = False
            results["error"] = exc
        finally:
            barrier.set()

    # Put the cache into "loading" state first so the waiter parks on the
    # condition. Then trigger the actual load in a daemon thread.
    genome_cache._state["status"] = "loading"
    genome_cache._state["path"] = str(fasta)

    t = threading.Thread(target=waiter, daemon=True)
    t.start()
    assert started.wait(timeout=2.0), "waiter never started"
    # Give the waiter time to actually call wait().
    time.sleep(0.1)

    # Trigger a load that updates state to ready (and notifies waiters).
    genome_cache._load_blocking(str(fasta))

    assert barrier.wait(timeout=3.0), "waiter never woke up"
    assert results.get("ok") is True, results
    assert results["genome"] == {"chr1": "ACGT"}
