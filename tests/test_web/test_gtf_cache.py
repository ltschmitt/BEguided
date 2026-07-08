"""Unit tests for ``beguided.web.gtf_cache``."""

from __future__ import annotations

import threading
import time
from unittest.mock import patch

import pytest


@pytest.fixture(autouse=True)
def _reset_module_state():
    """Ensure module state is clean for every test."""
    from beguided.web import gtf_cache

    gtf_cache._reset_for_test()
    yield
    gtf_cache._reset_for_test()


def _write_gtf(path, chroms: dict[str, list[tuple[int, int]]]) -> None:
    """Write a minimal GTF with MANE_Select CDS records.

    Args:
        path: Destination file path.
        chroms: Mapping of chromosome name -> list of (start, end) CDS intervals
            (1-based, inclusive). All transcripts use strand "+" and
            transcript_id "T1".
    """
    lines = [
        '##description: test',
        '##provider: beguided-tests',
    ]
    for chrom, intervals in chroms.items():
        for start, end in intervals:
            lines.append(
                "\t".join([
                    chrom,
                    "test",
                    "CDS",
                    str(start),
                    str(end),
                    ".",
                    "+",
                    "0",
                    f'transcript_id "T1"; tag "MANE_Select";',
                ])
            )
    path.write_text("\n".join(lines) + "\n")


def test_preload_async_sets_ready_state(tmp_path):
    from beguided.web import gtf_cache

    gtf = tmp_path / "tiny.gtf"
    _write_gtf(gtf, {"chr1": [(100, 200), (300, 400)], "chr2": [(10, 50)]})

    gtf_cache.preload_async(str(gtf))
    # Wait up to 2s for the daemon thread to finish.
    deadline = time.monotonic() + 2.0
    while time.monotonic() < deadline:
        if gtf_cache.is_ready():
            break
        time.sleep(0.02)

    status = gtf_cache.get_status()
    assert status["status"] == "ready"
    assert status["chromosome_count"] == 2
    assert status["path"] == str(gtf)

    trees = gtf_cache.wait_until_ready(timeout=1.0)
    assert set(trees.keys()) == {"chr1", "chr2", "__tx_codons__"}
    # Each interval is converted to 0-based half-open [start, end).
    assert list(trees["chr1"][99:200])  # noqa: E501 - probe the tree
    assert list(trees["chr1"][299:400])
    assert list(trees["chr2"][9:50])


def test_preload_async_is_idempotent(tmp_path):
    from beguided.web import gtf_cache

    gtf = tmp_path / "tiny.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    gtf_cache.preload_async(str(gtf))
    # Call a second time with the same path; must not raise.
    gtf_cache.preload_async(str(gtf))
    deadline = time.monotonic() + 2.0
    while not gtf_cache.is_ready():
        time.sleep(0.02)
    assert gtf_cache.is_ready()


def test_preload_async_rejects_different_path(tmp_path):
    from beguided.web import gtf_cache

    a = tmp_path / "a.gtf"
    b = tmp_path / "b.gtf"
    _write_gtf(a, {"chr1": [(1, 10)]})
    _write_gtf(b, {"chr1": [(1, 10)]})

    gtf_cache.preload_async(str(a))
    with pytest.raises(RuntimeError, match="already loaded"):
        gtf_cache.preload_async(str(b))
    gtf_cache._reset_for_test()


def test_wait_until_ready_raises_after_error(tmp_path):
    """Simulate a load failure: ``load_gtf_to_interval_tree`` raises and
    ``os._exit`` is intercepted so the test process keeps running."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "broken.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    # Mark the preload as in-flight so ``wait_until_ready`` will block.
    # We then fail the load by raising from a synchronous helper rather than
    # driving the daemon thread (which calls os._exit and would kill pytest).
    gtf_cache._state["status"] = "loading"
    gtf_cache._state["path"] = str(gtf)

    with patch(
        "beguided.web.gtf_cache.load_gtf_to_interval_tree",
        side_effect=RuntimeError("boom"),
    ):
        # Simulate the daemon thread's failure path: it sets state, notifies
        # waiters, and exits the process. We replay that here synchronously.
        def _simulate_failure():
            with gtf_cache._cond:
                gtf_cache._state["status"] = "error"
                gtf_cache._state["error"] = "RuntimeError: boom"
                gtf_cache._state["finished_at"] = time.time()
                gtf_cache._cond.notify_all()

        # Start the waiter first.
        results = {}

        def waiter():
            try:
                results["trees"] = gtf_cache.wait_until_ready(timeout=2.0)
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
        assert gtf_cache.get_status()["status"] == "error"


def test_wait_until_ready_raises_when_not_started():
    from beguided.web import gtf_cache

    with pytest.raises(RuntimeError, match="not been started"):
        gtf_cache.wait_until_ready(timeout=0.1)


def test_wait_until_ready_times_out(tmp_path):
    """If preload never finishes within timeout, raise TimeoutError."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "slow.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    sentinel = threading.Event()

    def slow_load(_path):
        sentinel.wait(timeout=5.0)
        return {"chr1": []}

    with patch("beguided.web.gtf_cache.load_gtf_to_interval_tree", side_effect=slow_load):
        gtf_cache.preload_async(str(gtf))
        with pytest.raises(TimeoutError):
            gtf_cache.wait_until_ready(timeout=0.2)

        # Release the loader so the daemon thread can exit cleanly.
        sentinel.set()


def test_wait_blocks_then_returns(tmp_path):
    """A caller waiting on the condition wakes up once preload finishes."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "tiny.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    barrier = threading.Event()
    started = threading.Event()
    results: dict = {}

    def waiter():
        try:
            started.set()
            results["trees"] = gtf_cache.wait_until_ready(timeout=5.0)
            results["ok"] = True
        except Exception as exc:  # pragma: no cover - test failure path
            results["ok"] = False
            results["error"] = exc
        finally:
            barrier.set()

    # Put the cache into "loading" state first so the waiter parks on the
    # condition. Then trigger the actual load in a daemon thread.
    gtf_cache._state["status"] = "loading"
    gtf_cache._state["path"] = str(gtf)

    t = threading.Thread(target=waiter, daemon=True)
    t.start()
    assert started.wait(timeout=2.0), "waiter never started"
    # Give the waiter time to actually call wait().
    time.sleep(0.1)

    # Trigger a load that updates state to ready (and notifies waiters).
    gtf_cache._load_blocking(str(gtf))

    assert barrier.wait(timeout=3.0), "waiter never woke up"
    assert results.get("ok") is True, results
    assert "chr1" in results["trees"]
