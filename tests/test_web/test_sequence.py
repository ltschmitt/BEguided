"""Tests for the sequence input parser and the /sequence/parse endpoint."""

from __future__ import annotations

import base64

import pandas as pd
import pytest

pytest.importorskip("fastapi")

from beguided.web.sequence_input import (
    DEFAULT_CHROMOSOME,
    MAX_SEQUENCE_LENGTH,
    MIN_SEQUENCE_LENGTH,
    SequenceParseError,
    parse_sequence_text,
    records_to_fasta,
    records_to_targets_dataframe,
)


# A canonical 40-nt STRIPPED sequence (what the parser will produce).
# Position 17 (1-based) is 'A' and is the default edit site.
CANONICAL_STRIPPED = "ATAGAGTATAGACCCAA" + "CCTATCCAGTCTACCTAGGGAATC"
CANONICAL_EDIT_POS = 17  # The 'A' that the bracket will become.


def _wrap(stripped: str, edit_base: str | None = None, desired_base: str = "G",
          header: str | None = None, edit_position: int | None = None) -> str:
    """Build user-facing input: replace the character at ``edit_position``
    (1-based) of ``stripped`` with a ``[X>Y]`` bracket.

    The default ``edit_base`` is the character already at ``edit_position``
    in ``stripped`` (so the bracket correctly represents the existing
    base). The caller can override either parameter to construct specific
    edit directions.

    Optionally prepends a FASTA header.
    """
    pos = edit_position or CANONICAL_EDIT_POS
    base = edit_base if edit_base is not None else stripped[pos - 1]
    bracketed = (
        stripped[:pos - 1]
        + f"[{base}>{desired_base}]"
        + stripped[pos:]
    )
    if header:
        return f">{header}\n{bracketed}\n"
    return bracketed + "\n"


# ---------------------------------------------------------------------------
# parse_sequence_text: pure-parser tests
# ---------------------------------------------------------------------------


class TestParseSequenceText:
    def test_single_sequence_no_header(self):
        text = _wrap(CANONICAL_STRIPPED)
        records = parse_sequence_text(text)
        assert len(records) == 1
        r = records[0]
        assert r.chromosome == DEFAULT_CHROMOSOME
        assert r.sequence == CANONICAL_STRIPPED
        assert r.edit_position == CANONICAL_EDIT_POS
        assert r.target_base == "A"
        assert r.desired_base == "G"
        assert r.genome_base == "A"

    def test_bracket_becomes_target_base_at_position(self):
        # The character at edit_position in the stripped sequence must be
        # the bracket's target_base, regardless of what surrounded the
        # bracket in the user input.
        text = "ATAG[T>C]GAC" + "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        # Pad to >= MIN_SEQUENCE_LENGTH
        text = text + "GTACGTACGTACGTACGTACGTACGTACGT"
        records = parse_sequence_text(text + "\n")
        assert records[0].sequence[records[0].edit_position - 1] == "T"
        assert records[0].target_base == "T"
        assert records[0].desired_base == "C"

    def test_text_after_bracket_shifts_one_position(self):
        # The character originally right after `]` should appear at
        # edit_position + 1 in the stripped sequence.
        text = "ATAG[T>C]GAC" + "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        records = parse_sequence_text(text + "\n")
        assert records[0].sequence[records[0].edit_position] == "G"
        assert records[0].sequence[records[0].edit_position + 1] == "A"
        assert records[0].sequence[records[0].edit_position + 2] == "C"

    def test_single_sequence_with_header(self):
        text = _wrap(CANONICAL_STRIPPED, header="my_target")
        records = parse_sequence_text(text)
        assert len(records) == 1
        assert records[0].chromosome == "my_target"

    def test_header_takes_first_whitespace_token(self):
        text = _wrap(CANONICAL_STRIPPED, header="chr1 extra annotation here")
        records = parse_sequence_text(text)
        assert records[0].chromosome == "chr1"

    def test_multi_sequence_fasta(self):
        # Position 17 in CANONICAL_STRIPPED is 'A' → default [A>G].
        # Position 13 is 'C' → use [C>T] for the second record.
        text = (
            _wrap(CANONICAL_STRIPPED, header="t1")
            + _wrap(
                CANONICAL_STRIPPED, header="t2",
                edit_base="C", desired_base="T",
                edit_position=13,
            )
        )
        records = parse_sequence_text(text)
        assert [r.chromosome for r in records] == ["t1", "t2"]
        assert [r.desired_base for r in records] == ["G", "T"]

    def test_lowercase_bracket_normalized(self):
        text = _wrap(CANONICAL_STRIPPED, edit_base="a", desired_base="g")
        records = parse_sequence_text(text)
        assert records[0].target_base == "A"
        assert records[0].desired_base == "G"

    def test_bracket_with_whitespace(self):
        seq = (
            CANONICAL_STRIPPED[:CANONICAL_EDIT_POS - 1]
            + "[ A > G ]"
            + CANONICAL_STRIPPED[CANONICAL_EDIT_POS:]
        )
        records = parse_sequence_text(seq)
        assert records[0].target_base == "A"
        assert records[0].desired_base == "G"

    def test_all_four_directions_accepted(self):
        # A>G, T>C, C>T, G>A are the four supported editing directions
        # (sense + antisense complements of ABE and CBE).
        for edit in ("A>G", "T>C", "C>T", "G>A"):
            # Build a stripped sequence where position 17 carries
            # ``edit[0]`` (the bracket target_base).
            target_base = edit[0]
            desired_base = edit[2]
            stripped = (
                CANONICAL_STRIPPED[:CANONICAL_EDIT_POS - 1]
                + target_base
                + CANONICAL_STRIPPED[CANONICAL_EDIT_POS:]
            )
            text = _wrap(stripped, edit_base=target_base, desired_base=desired_base)
            records = parse_sequence_text(text)
            assert records[0].target_base == target_base
            assert records[0].desired_base == desired_base
            # And the stripped sequence matches what we built.
            assert records[0].sequence == stripped

    def test_missing_bracket_rejected(self):
        text = CANONICAL_STRIPPED + "\n"
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "bracket" in str(exc.value).lower()

    def test_multiple_brackets_rejected(self):
        seq = (
            CANONICAL_STRIPPED[:10]
            + "[A>G]"
            + CANONICAL_STRIPPED[10:30]
            + "[C>T]"
            + CANONICAL_STRIPPED[30:]
        )
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(seq)
        assert "multiple" in str(exc.value).lower()

    def test_malformed_bracket_rejected(self):
        bad_brackets = ["[A-G]", "[AG]", "[N>G]", "[A>]", "[>G]"]
        for bad in bad_brackets:
            text = (
                CANONICAL_STRIPPED[:CANONICAL_EDIT_POS - 1]
                + bad
                + CANONICAL_STRIPPED[CANONICAL_EDIT_POS:]
            )
            with pytest.raises(SequenceParseError):
                parse_sequence_text(text)

    def test_unsupported_direction_rejected(self):
        # A>T is not supported by any base editor profile's sense or antisense direction.
        text = (
            CANONICAL_STRIPPED[:CANONICAL_EDIT_POS - 1]
            + "[A>T]"
            + CANONICAL_STRIPPED[CANONICAL_EDIT_POS:]
        )
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "not supported" in str(exc.value).lower()

    def test_empty_input_rejected(self):
        with pytest.raises(SequenceParseError):
            parse_sequence_text("")
        with pytest.raises(SequenceParseError):
            parse_sequence_text("   \n  \n")

    def test_sequence_too_short_rejected(self):
        # Stripped length below MIN_SEQUENCE_LENGTH. Build a 10-nt stripped
        # sequence and wrap it with a bracket — the resulting stripped
        # length is still 10 (the bracket itself is replaced by one base).
        stripped = "A" * 10
        text = stripped[:5] + "[A>G]" + stripped[5:] + "\n"
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "too short" in str(exc.value).lower()

    def test_sequence_too_long_rejected(self):
        # Stripped length = MAX_SEQUENCE_LENGTH + 1, plus the bracket
        # itself.
        stripped = "A" * (MAX_SEQUENCE_LENGTH + 1)
        text = stripped[:10] + "[A>G]" + stripped[10:] + "\n"
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "too long" in str(exc.value).lower()

    def test_non_acgt_characters_rejected(self):
        # 'N' is not allowed in the stripped sequence.
        stripped = "ATAGANTATAGACCCAA" + "CCTATCCAGTCTACCTAGGGAATC"
        text = _wrap(stripped) + "\n"
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "non-acgt" in str(exc.value).lower()

    def test_empty_header_in_fasta_rejected(self):
        text = f">\n{CANONICAL_STRIPPED}\n"
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert "header is empty" in str(exc.value).lower()

    def test_record_index_set_on_error(self):
        # Second record has a bad bracket; record_index should be 1.
        text = (
            _wrap(CANONICAL_STRIPPED, header="t1")
            + f">bad\n{CANONICAL_STRIPPED}\n"  # no bracket
        )
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert exc.value.record_index == 1

    def test_duplicate_header_rejected(self):
        text = (
            _wrap(CANONICAL_STRIPPED, header="t1")
            + _wrap(CANONICAL_STRIPPED, header="t1")
        )
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert exc.value.record_index == 1
        assert "duplicate" in str(exc.value).lower()
        assert "'t1'" in str(exc.value)

    def test_header_collapsing_to_same_name_rejected(self):
        text = (
            _wrap(CANONICAL_STRIPPED, header="custom")
            + _wrap(CANONICAL_STRIPPED, header="custom extra")
        )
        with pytest.raises(SequenceParseError) as exc:
            parse_sequence_text(text)
        assert exc.value.record_index == 1
        assert "duplicate" in str(exc.value).lower()

    def test_edit_position_one_based(self):
        # Position 1 means the bracket is at the start.
        seq = "A" + "G" * 30
        text = f"[A>G]{seq}\n"
        records = parse_sequence_text(text)
        assert records[0].edit_position == 1


# ---------------------------------------------------------------------------
# records_to_targets_dataframe + records_to_fasta
# ---------------------------------------------------------------------------


class TestRecordsSerialization:
    def test_targets_dataframe_has_core_columns(self):
        records = parse_sequence_text(_wrap(CANONICAL_STRIPPED))
        df = records_to_targets_dataframe(records)
        assert list(df["chromosome"]) == [DEFAULT_CHROMOSOME]
        assert list(df["position"]) == [CANONICAL_EDIT_POS]
        assert list(df["target_base"]) == ["A"]
        assert list(df["desired_base"]) == ["G"]
        assert "sequence" in df.columns
        assert "sequence_preview" in df.columns

    def test_fasta_round_trip(self):
        records = parse_sequence_text(
            _wrap(CANONICAL_STRIPPED, header="t1")
            + _wrap(
                CANONICAL_STRIPPED, header="t2",
                edit_base="C", desired_base="T",
                edit_position=13,
            )
        )
        fasta = records_to_fasta(records)
        assert fasta.startswith(">t1\n")
        assert ">t2\n" in fasta
        # The bracket-stripped sequence must not contain any '[' or ']'.
        for line in fasta.splitlines():
            if line.startswith(">"):
                continue
            assert "[" not in line
            assert "]" not in line


# ---------------------------------------------------------------------------
# /sequence/parse + /run: FastAPI-level integration
# ---------------------------------------------------------------------------


@pytest.fixture
def cfg_and_db(tmp_path, monkeypatch):
    """Minimal WebConfig + scheduler for FastAPI integration tests."""
    from beguided.web import config as cfg_mod
    from beguided.web import scheduler as sched_mod
    from beguided.web.config import WebConfig, set_config
    from beguided.web.scheduler import Scheduler, install

    db = tmp_path / "jobs.db"
    presets = tmp_path / "presets.tsv"
    presets.write_text(
        "name\tediting\tgeneral_edit_start\tgeneral_edit_end\t"
        "high_edit_start\thigh_edit_end\tpam\tcas\tsource\n"
        "ABE8e\tABE\t2\t13\t3\t8\tNGG\tSpCas9\tNature2022\n"
    )
    set_config(WebConfig(
        base_editor_profiles_file=str(presets),
        genome_fasta="/nonexistent",
        gtf_path="/nonexistent",
        beguidedoff_genome="",
        max_input_targets=100,
        db_path=str(db),
        queue_max=3,
        cooldown_seconds=0,
        max_upload_mb=5,
        usage_stats_enabled=False,
    ))
    sch = Scheduler(db_path=str(db), queue_max=3, cooldown_seconds=0)
    install(sch)
    sch.start()
    yield cfg_mod.get_config(), sch


def _client():
    from fastapi.testclient import TestClient
    from beguided.web.app import create_app
    return TestClient(create_app())


class TestSequenceParseEndpoint:
    def test_single_sequence_returns_targets_and_fasta(self, cfg_and_db):
        res = _client().post(
            "/sequence/parse",
            data={"text": _wrap(CANONICAL_STRIPPED)},
        )
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is True
        assert body["row_count"] == 1
        assert body["raw_b64"]
        assert body["custom_genome_b64"]

        targets_tsv = base64.b64decode(body["raw_b64"]).decode()
        assert "custom\t17\tA\tG" in targets_tsv

        fasta = base64.b64decode(body["custom_genome_b64"]).decode()
        assert fasta.startswith(">custom\n")
        assert "[" not in fasta  # bracket stripped
        assert CANONICAL_STRIPPED in fasta  # stripped sequence present

    def test_multi_sequence_fasta(self, cfg_and_db):
        text = (
            _wrap(CANONICAL_STRIPPED, header="t1")
            + _wrap(CANONICAL_STRIPPED, header="t2")
        )
        res = _client().post("/sequence/parse", data={"text": text})
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is True
        assert body["row_count"] == 2

        preview = body["preview"]
        assert preview[0]["chromosome"] == "t1"
        assert preview[1]["chromosome"] == "t2"

    def test_missing_bracket_reported(self, cfg_and_db):
        res = _client().post(
            "/sequence/parse", data={"text": CANONICAL_STRIPPED + "\n"}
        )
        assert res.status_code == 200
        body = res.json()
        assert body["ok"] is False
        assert "bracket" in body["error"].lower()

    def test_multiple_brackets_reported(self, cfg_and_db):
        seq = (
            CANONICAL_STRIPPED[:10]
            + "[A>G]"
            + CANONICAL_STRIPPED[10:30]
            + "[C>T]"
            + CANONICAL_STRIPPED[30:]
        )
        res = _client().post("/sequence/parse", data={"text": seq + "\n"})
        body = res.json()
        assert body["ok"] is False
        assert "multiple" in body["error"].lower()

    def test_empty_input_reported(self, cfg_and_db):
        res = _client().post("/sequence/parse", data={"text": ""})
        body = res.json()
        assert body["ok"] is False

    def test_preview_contains_sequence_snippet(self, cfg_and_db):
        res = _client().post(
            "/sequence/parse", data={"text": _wrap(CANONICAL_STRIPPED)}
        )
        body = res.json()
        assert body["ok"] is True
        assert "[" in body["preview"][0]["sequence_preview"]
        assert "]" in body["preview"][0]["sequence_preview"]


class TestRunWithSequenceMode:
    def test_run_in_sequence_mode_enqueues_job(self, cfg_and_db, monkeypatch):
        """`/run` with mode=sequence should accept the custom_genome_b64
        payload and enqueue a job that runs only the guidefinder step."""
        # Parse first to get the b64 blobs.
        cfg, sch = cfg_and_db
        parse_res = _client().post(
            "/sequence/parse", data={"text": _wrap(CANONICAL_STRIPPED)}
        )
        body = parse_res.json()
        assert body["ok"] is True

        # Capture the kwargs that tasks.execute_job receives so we can
        # verify mode/custom_genome_b64 made it through.
        captured = {}

        def fake_execute_job(**kwargs):
            captured.update(kwargs)
            return {
                "job_id": kwargs["job_id"],
                "status": "complete",
                "message": "fake ok",
                "records": [{"chromosome": "custom", "position": 17}],
            }

        import beguided.web.scheduler as sched_mod
        sched_mod.execute_job = fake_execute_job

        run_res = _client().post(
            "/run",
            data={
                "raw_b64": body["raw_b64"],
                "base_editor": "ABE8e",
                "mode": "sequence",
                "custom_genome_b64": body["custom_genome_b64"],
            },
        )
        assert run_res.status_code == 200, run_res.text
        run_body = run_res.json()
        assert run_body["ok"] is True

        # Wait for the worker thread to invoke our fake.
        for _ in range(20):
            if captured:
                break
            import time
            time.sleep(0.05)
        assert captured.get("mode") == "sequence"
        assert captured.get("custom_genome_b64") == body["custom_genome_b64"]

    def test_run_sequence_mode_requires_custom_genome_b64(self, cfg_and_db):
        res = _client().post(
            "/run",
            data={
                "raw_b64": base64.b64encode(b"chrom\t1\tA\tG\n").decode(),
                "base_editor": "ABE8e",
                "mode": "sequence",
                "custom_genome_b64": "",
            },
        )
        assert res.status_code == 400
        body = res.json()
        assert body["ok"] is False
        assert "custom_genome_b64" in body["error"]

    def test_run_rejects_unknown_mode(self, cfg_and_db):
        res = _client().post(
            "/run",
            data={
                "raw_b64": base64.b64encode(b"chrom\t1\tA\tG\n").decode(),
                "base_editor": "ABE8e",
                "mode": "wat",
            },
        )
        assert res.status_code == 400

    def test_run_spreadsheet_mode_default(self, cfg_and_db):
        """Without a ``mode`` field the request should be treated as
        spreadsheet mode and not require custom_genome_b64."""
        captured = {}

        def fake_execute_job(**kwargs):
            captured.update(kwargs)
            return {
                "job_id": kwargs["job_id"],
                "status": "complete",
                "message": "fake ok",
                "records": [],
            }

        import beguided.web.scheduler as sched_mod
        sched_mod.execute_job = fake_execute_job

        res = _client().post(
            "/run",
            data={
                "raw_b64": base64.b64encode(b"chr1\t1\tA\tG\n").decode(),
                "base_editor": "ABE8e",
            },
        )
        assert res.status_code == 200
        for _ in range(20):
            if captured:
                break
            import time
            time.sleep(0.05)
        assert captured.get("mode", "spreadsheet") == "spreadsheet"
        assert captured.get("custom_genome_b64") in (None, "")


# ---------------------------------------------------------------------------
# End-to-end: parser output → load_genome → run_guidefinder.
# ---------------------------------------------------------------------------


class TestSequenceModeEndToEnd:
    def test_parsed_genome_loadable_by_load_genome(self, cfg_and_db):
        """The custom FASTA produced by the parser should be loadable by
        ``beguided.genome.io.load_genome`` and round-trip cleanly.
        """
        # Position 13 in CANONICAL_STRIPPED is 'C', so a [C>T] bracket
        # there is valid; position 17 is 'A' for [A>G].
        records = parse_sequence_text(
            _wrap(CANONICAL_STRIPPED, header="t1")
            + _wrap(
                CANONICAL_STRIPPED, header="t2",
                edit_base="C", desired_base="T",
                edit_position=13,
            )
        )
        fasta = records_to_fasta(records)

        import tempfile
        from beguided.genome.io import load_genome

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", delete=False
        ) as fh:
            fh.write(fasta)
            path = fh.name

        try:
            genome = load_genome(path)
        finally:
            import os
            os.unlink(path)

        assert set(genome.keys()) == {"t1", "t2"}
        # No bracket characters should leak into the loaded sequences.
        for seq in genome.values():
            assert "[" not in seq
            assert "]" not in seq

    def test_targets_dataframe_round_trip_via_run_guidefinder(self):
        """The targets DataFrame produced by the parser should be accepted by
        ``run_genome_wide`` without modification (it has the four core
        columns the guidefinder expects).
        """
        from pathlib import Path
        from beguided.guidefinder.core import run_genome_wide
        from beguided.editing.presets import load_base_editor_profiles
        from beguided.genome.io import load_genome
        import tempfile, os

        repo_root = Path(__file__).resolve().parents[2]
        presets_path = repo_root / "configs" / "base_editor_profiles.tsv"
        if not presets_path.exists():
            pytest.skip("configs/base_editor_profiles.tsv not found in repo root")

        records = parse_sequence_text(_wrap(CANONICAL_STRIPPED, header="t1"))
        targets_df = records_to_targets_dataframe(records)
        # The guidefinder only consumes the four core columns.
        core = targets_df[["chromosome", "position", "target_base", "desired_base"]]

        fasta = records_to_fasta(records)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", delete=False
        ) as fh:
            fh.write(fasta)
            path = fh.name
        try:
            genome = load_genome(path)
        finally:
            os.unlink(path)

        profiles = load_base_editor_profiles(str(presets_path))
        base_editor_names = sorted({
            name
            for editor in profiles.values()
            for name in editor.keys()
        })

        results = run_genome_wide(genome, core, profiles, base_editor_names)
        # We don't assert a specific guide; just that the pipeline didn't
        # reject the input outright.
        assert isinstance(results, list)

    def test_end_to_end_user_reporting_a_rich_sequences(self):
        """Regression: the two A/T-rich sequences reported by a user
        produced 29 guides in step 1 but step 3 dropped all of them
        because the filter_mode='no_bystanders' shortcut also acted as a
        bystander filter on bystander_positions, which is empty for nearly
        no guide on A-rich sequences. Sequence-mode jobs must keep all
        guides generated in step 1.
        """
        from pathlib import Path
        from beguided.pipeline import run_pipeline

        repo_root = Path(__file__).resolve().parents[2]
        presets_path = repo_root / "configs" / "base_editor_profiles.tsv"
        if not presets_path.exists():
            pytest.skip("configs/base_editor_profiles.tsv not found in repo root")

        user_text = (
            ">target_C\n"
            "AGAGGGGCATAGCCATGTAGCATC[A>G]CCAATTGGCCATGCGGTGGCCGAGA\n"
            ">target_T\n"
            "GAGATTACCGGTACCCATGGCCCT[T>C]TATTGCCAGGCATTGGAGGGAAG\n"
        )
        records = parse_sequence_text(user_text)
        assert len(records) == 2

        targets_df = records_to_targets_dataframe(records)
        core = targets_df[["chromosome", "position", "target_base", "desired_base"]]

        fasta = records_to_fasta(records)

        import tempfile
        with tempfile.TemporaryDirectory() as work_dir:
            targets_path = Path(work_dir) / "targets.tsv"
            genome_path = Path(work_dir) / "custom.fa"
            targets_path.write_text(core.to_csv(sep="\t", index=False))
            genome_path.write_text(fasta)

            config = {
                "input": {
                    "targets": str(targets_path),
                    "genome": str(genome_path),
                    "presets": str(presets_path),
                    "gff_file": None,
                },
                "output_prefix": work_dir,
                "use_base_editors": "all",
                "keep_bystanders": "all",
                "guide_score_threshold": 0.0,
                "search_offtargets": False,
                "skip_bystander": True,
            }
            run_pipeline(config, start_step=1, end_step=5)
            out = pd.read_csv(Path(work_dir) / "beguided_results.tsv", sep="\t")

        # We must not have dropped all guides. With these two A/T-rich
        # sequences we expect well over zero guides to survive.
        assert len(out) > 0
        # And the two target chromosomes should both be represented.
        assert set(out["chromosome"]) == {"target_C", "target_T"}
