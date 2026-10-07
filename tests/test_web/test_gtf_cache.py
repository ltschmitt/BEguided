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


def _write_gtf(
    path,
    chroms: dict[str, list[tuple[int, int]]],
    tag: str = "MANE_Select",
) -> None:
    """Write a minimal GTF with canonical CDS records.

    Args:
        path: Destination file path.
        chroms: Mapping of chromosome name -> list of (start, end) CDS intervals
            (1-based, inclusive). All transcripts use strand "+" and
            transcript_id "T1".
        tag: CDS attribute tag to emit on every record.
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
                    f'transcript_id "T1"; tag "{tag}";',
                ])
            )
    path.write_text("\n".join(lines) + "\n")


def _wait_until(key: str, status: str, timeout: float = 2.0) -> dict:
    from beguided.web import gtf_cache

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        snap = gtf_cache.get_status(key)
        if snap["status"] == status:
            return snap
        time.sleep(0.02)
    raise AssertionError(f"key {key!r} never reached status {status!r}")


def test_derive_gtf_label():
    from beguided.web.gtf_cache import derive_gtf_label

    assert derive_gtf_label("/x/Homo_sapiens.GRCh38.112.gtf") == "GRCh38.112"
    assert (
        derive_gtf_label("/x/Homo_sapiens.GRCh38.chr_patch_hapl_scaff.gtf")
        == "GRCh38"
    )
    assert (
        derive_gtf_label("/x/Mus_musculus.GRCm39.2026_04.gtf")
        == "GRCm39.2026_04"
    )
    assert derive_gtf_label("/x/anything.gtf") == "GRCh38"


def test_preload_async_sets_ready_state(tmp_path):
    from beguided.web import gtf_cache

    gtf = tmp_path / "Homo_sapiens.GRCh38.112.gtf"
    _write_gtf(gtf, {"chr1": [(100, 200), (300, 400)], "chr2": [(10, 50)]})

    gtf_cache.preload_async("human", str(gtf))
    status = _wait_until("human", "ready")

    assert status["chromosome_count"] == 2
    assert status["path"] == str(gtf)
    assert status["label"] == "GRCh38.112"
    assert status["select_tag"] == "MANE_Select"

    trees = gtf_cache.wait_until_ready("human", timeout=1.0)
    assert set(trees.keys()) == {"chr1", "chr2", "__tx_codons__"}
    # Each interval is converted to 0-based half-open [start, end).
    assert list(trees["chr1"][99:200])  # noqa: E501 - probe the tree
    assert list(trees["chr1"][299:400])
    assert list(trees["chr2"][9:50])


def test_preload_async_honours_select_tag(tmp_path):
    """The mouse annotation uses Ensembl_canonical rather than MANE_Select."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "Mus_musculus.GRCm39.2026_04.gtf"
    _write_gtf(gtf, {"1": [(1, 30)]}, tag="Ensembl_canonical")

    gtf_cache.preload_async(
        "mouse", str(gtf), select_tag="Ensembl_canonical"
    )
    status = _wait_until("mouse", "ready")
    assert status["select_tag"] == "Ensembl_canonical"

    trees = gtf_cache.wait_until_ready("mouse", timeout=1.0)
    assert list(trees["1"][0:30])


def test_preload_async_is_idempotent(tmp_path):
    from beguided.web import gtf_cache

    gtf = tmp_path / "tiny.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    gtf_cache.preload_async("human", str(gtf))
    # Call a second time with the same path; must not raise.
    gtf_cache.preload_async("human", str(gtf))
    _wait_until("human", "ready")
    assert gtf_cache.is_ready("human")


def test_preload_async_rejects_different_path(tmp_path):
    from beguided.web import gtf_cache

    a = tmp_path / "a.gtf"
    b = tmp_path / "b.gtf"
    _write_gtf(a, {"chr1": [(1, 10)]})
    _write_gtf(b, {"chr1": [(1, 10)]})

    gtf_cache.preload_async("human", str(a))
    with pytest.raises(RuntimeError, match="already loaded"):
        gtf_cache.preload_async("human", str(b))


def test_wait_until_ready_raises_after_error(tmp_path):
    """Simulate a load failure and verify waiters surface the error."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "broken.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    entry = gtf_cache._get_or_create("human")
    with entry.cond:
        entry.status = "loading"
        entry.path = str(gtf)

    def _simulate_failure():
        with entry.cond:
            entry.status = "error"
            entry.error = "RuntimeError: boom"
            entry.finished_at = time.time()
            entry.cond.notify_all()

    results = {}

    def waiter():
        try:
            results["trees"] = gtf_cache.wait_until_ready("human", timeout=2.0)
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
    assert gtf_cache.get_status("human")["status"] == "error"


def test_wait_until_ready_raises_when_not_started():
    from beguided.web import gtf_cache

    with pytest.raises(RuntimeError, match="not been started"):
        gtf_cache.wait_until_ready(timeout=0.1)
    with pytest.raises(RuntimeError, match="not been preloaded"):
        gtf_cache.wait_until_ready("human", timeout=0.1)


def test_wait_until_ready_times_out(tmp_path):
    """If preload never finishes within timeout, raise TimeoutError."""
    from beguided.web import gtf_cache

    gtf = tmp_path / "slow.gtf"
    _write_gtf(gtf, {"chr1": [(1, 10)]})

    sentinel = threading.Event()

    def slow_load(*_args, **_kwargs):
        sentinel.wait(timeout=5.0)
        return {"chr1": []}

    with patch("beguided.web.gtf_cache.load_gtf_to_interval_tree", side_effect=slow_load):
        gtf_cache.preload_async("human", str(gtf))
        with pytest.raises(TimeoutError):
            gtf_cache.wait_until_ready("human", timeout=0.2)

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
            results["trees"] = gtf_cache.wait_until_ready("human", timeout=5.0)
            results["ok"] = True
        except Exception as exc:  # pragma: no cover - test failure path
            results["ok"] = False
            results["error"] = exc
        finally:
            barrier.set()

    # Put the key into "loading" state first so the waiter parks on the
    # condition, then drive the real load synchronously.
    entry = gtf_cache._get_or_create("human")
    with entry.cond:
        entry.status = "loading"
        entry.path = str(gtf)

    t = threading.Thread(target=waiter, daemon=True)
    t.start()
    assert started.wait(timeout=2.0), "waiter never started"
    # Give the waiter time to actually call wait().
    time.sleep(0.1)

    gtf_cache._load_blocking(entry, False)

    assert barrier.wait(timeout=3.0), "waiter never woke up"
    assert results.get("ok") is True, results
    assert "chr1" in results["trees"]
