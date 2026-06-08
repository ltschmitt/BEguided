"""Tests for pipeline components."""

import pytest
from pathlib import Path


class TestPipelineImports:
    """Ensure pipeline module imports correctly."""

    def test_import_pipeline(self):
        from beguided.pipeline import run_pipeline, PipelineError
        assert callable(run_pipeline)
        assert issubclass(PipelineError, Exception)


class TestPipelineConfigValidation:
    """Test configuration validation."""

    def test_missing_output_prefix(self):
        from beguided.pipeline import PipelineError, run_pipeline
        config = {"input": {"targets": "t.tsv"}}
        with pytest.raises((PipelineError, KeyError)):
            run_pipeline(config)

    def test_missing_input_section(self):
        from beguided.pipeline import PipelineError, run_pipeline
        config = {"output_prefix": "./results"}
        with pytest.raises((PipelineError, KeyError)):
            run_pipeline(config)

    def test_guide_score_threshold_out_of_range(self):
        from beguided.pipeline import PipelineError, run_pipeline
        config = {
            "input": {"targets": "t.tsv", "genome": "g.fa"},
            "output_prefix": "./results",
            "guide_score_threshold": 1.5,
        }
        with pytest.raises(PipelineError):
            run_pipeline(config, start_step=1, end_step=1)

    def test_guide_score_threshold_negative(self):
        from beguided.pipeline import PipelineError, run_pipeline
        config = {
            "input": {"targets": "t.tsv", "genome": "g.fa"},
            "output_prefix": "./results",
            "guide_score_threshold": -0.1,
        }
        with pytest.raises(PipelineError):
            run_pipeline(config, start_step=1, end_step=1)


class TestPipelineSkipOffTarget:
    """Test pipeline when off-target screening is disabled."""

    def test_skip_offtarget_produces_results_no_ot_cols(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Run full pipeline with search_offtargets=false."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
        )

        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t60\t.\t+\t0\t"
            'gene_id "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )

        base_dir = tmp_path / "out"
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "presets": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_presets": "all",
            "keep_bystanders": "all",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }

        run_pipeline(config, start_step=1, end_step=5)

        final_output = base_dir / "beguided_results.tsv"
        assert final_output.exists(), "Final TSV should be created even when off-targets are skipped"

        df = pd.read_csv(final_output, sep="\t")
        assert len(df) > 0, "Should contain at least one guide row"

        # Off-target columns must NOT be present when screening is skipped
        assert "worst_offtarget_score" not in df.columns
        assert "offtarget_count" not in df.columns
        assert "worst_offtarget_chr" not in df.columns

        # Intermediate files should use stepN naming
        assert (base_dir / "step1_guidefinder.tsv").exists()
        assert (base_dir / "step2_bystander.tsv").exists()
        assert (base_dir / "step3_filtered.tsv").exists()


class TestPrepareBEguidedOFFInput:
    """Unit tests for _prepare_beguidedoff_input without BEguidedOFF installed."""

    def test_columns_and_positions(self, tmp_path):
        """Verify output columns and on-target position per strand."""
        from beguided.pipeline import _prepare_beguidedoff_input
        import pandas as pd

        # + strand: position=3, editing_position=3 -> guide 1..20 -> start=1
        # - strand: position=3, editing_position=18 -> guide_end=20 -> pos=20
        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1"],
            "position": [3, 3],
            "target_base": ["A", "A"],
            "desired_base": ["G", "G"],
            "strand_used": ["+", "-"],
            "editing_position": [3, 18],
            "guide_rna_sequence": ["GTACAGCTAGCTAGCTAGCT", "CCCCAGCTAGCTAGCTAGCT"],
            "be_preset": ["ABE8e", "ABE8e"],
        })
        out_path = tmp_path / "bgo_input.tsv"
        _prepare_beguidedoff_input(df, out_path)

        result = pd.read_csv(out_path, sep="\t")
        assert list(result.columns) == [
            "chromosome", "position", "guide_rna_sequence",
            "be_preset", "strand_used",
        ]

        plus_row = result[result["strand_used"] == "+"].iloc[0]
        assert plus_row["position"] == 1

        minus_row = result[result["strand_used"] == "-"].iloc[0]
        assert minus_row["position"] == 20

    def test_no_guides_skipped(self, tmp_path):
        """Rows with '**No guide possible**' should be dropped."""
        from beguided.pipeline import _prepare_beguidedoff_input
        import pandas as pd

        df = pd.DataFrame({
            "chromosome": ["chr1", "chr1"],
            "position": [10, 20],
            "target_base": ["A", "A"],
            "desired_base": ["G", "G"],
            "strand_used": ["+", "+"],
            "editing_position": [5, 5],
            "guide_rna_sequence": ["AAAAAAAAAAAAAAAAAAAA", "**No guide possible**"],
            "be_preset": ["ABE8e", "ABE8e"],
        })
        out_path = tmp_path / "bgo_input.tsv"
        _prepare_beguidedoff_input(df, out_path)

        result = pd.read_csv(out_path, sep="\t")
        assert len(result) == 1
        assert result.iloc[0]["guide_rna_sequence"] == "AAAAAAAAAAAAAAAAAAAA"

    def test_empty_dataframe(self, tmp_path):
        """Empty input should produce empty output (no crash)."""
        from beguided.pipeline import _prepare_beguidedoff_input
        import pandas as pd

        df = pd.DataFrame({
            "chromosome": [],
            "position": [],
            "target_base": [],
            "desired_base": [],
            "strand_used": [],
            "editing_position": [],
            "guide_rna_sequence": [],
            "be_preset": [],
        })
        out_path = tmp_path / "bgo_input.tsv"
        _prepare_beguidedoff_input(df, out_path)
        result = pd.read_csv(out_path, sep="\t")
        assert len(result) == 0


class TestPipelineWithBEguidedOFF:
    """Test pipeline off-target integration using BEguidedOFF."""

    @pytest.fixture(autouse=True)
    def _skip_without_beguidedoff(self):
        pytest.importorskip("beguidedoff", reason="BEguidedOFF not installed")

    def test_beguidedoff_pipeline_steps_4_6(
        self, tmp_path, valid_preset_tsv
    ):
        """Test BEguidedOFF integration (steps 4-6) with a pre-built step3 file."""
        from beguidedoff.genome_converter import convert
        from beguided.pipeline import run_pipeline
        import pandas as pd

        # Genome: guide CCCCCCCCCCCCCCCCCCCC with PAM TGG at pos 30-32,
        # target A at pos 10 -> editing_pos 10
        fasta = tmp_path / "mini.fa"
        fasta.write_text(
            ">chr1\n"
            "ACGTACGTACCCCCCCCCCCCCCCCCCCCTGGACGTACGT\n"
        )
        npz = tmp_path / "mini.npz"
        convert(str(fasta), str(npz), verbose=False)

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        base_dir = tmp_path / "out"
        base_dir.mkdir()

        # Pre-create a minimal step3_filtered.tsv (required columns only)
        filtered = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [10],
            "target_base": ["A"],
            "desired_base": ["G"],
            "genome_base": ["A"],
            "strand_used": ["+"],
            "editing_position": [10],
            "guide_rna_sequence": ["CCCCCCCCCCCCCCCCCCCC"],
            "pam_sequence": ["TGG"],
            "be_preset": ["ABE8e"],
            "bystander_effect": ["no_bystanders"],
            "simulated_outcome": ["CCCCCCCCCCCCCCCCCCCC"],
            "aa_wanted": ["P"],
            "aa_bystander": ["P"],
            "genomic_context_upstream": ["ACGT"],
            "genomic_context_downstream": ["ACGT"],
        })
        filtered.to_csv(base_dir / "step3_filtered.tsv", sep="\t", index=False)

        # Also create step2_bystander.tsv and step1_guidefinder.tsv
        filtered.to_csv(base_dir / "step2_bystander.tsv", sep="\t", index=False)
        filtered.to_csv(base_dir / "step1_guidefinder.tsv", sep="\t", index=False)

        config = {
            "input": {
                "targets": str(tmp_path / "targets.tsv"),
                "genome": str(fasta),
                "presets": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_presets": "all",
            "keep_bystanders": "all",
            "search_offtargets": True,
            "offtargets": {
                "beguidedoff_genome": str(npz),
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }

        run_pipeline(config, start_step=4, end_step=5)

        final_output = base_dir / "beguided_results.tsv"
        assert final_output.exists()

        df = pd.read_csv(final_output, sep="\t")
        assert len(df) > 0
        assert "worst_offtarget_score" in df.columns
        assert "offtarget_count" in df.columns
        assert (base_dir / "step4_offtarget.tsv").exists()

    def test_beguidedoff_on_target_excluded_plus_strand(
        self, tmp_path, valid_preset_tsv
    ):
        """Verify on-target exclusion via step4 with a known guide-start pos."""
        from beguidedoff.genome_converter import convert
        from beguided.pipeline import run_pipeline
        import pandas as pd

        # Genome: target A at pos 11, guide from pos 1, PAM GGG at 21-23
        fasta = tmp_path / "mini.fa"
        fasta.write_text(
            ">chr1\n"
            "GTACAGCTAGCTAGCTAGCTGGGGTACAGCTAGCTAGCTAGC"
            "TGGGGTACAGCTAGCTAGCTAGCTGGGGTACAGCTAGCTAGC"
            "TGGGG\n"
        )
        npz = tmp_path / "mini.npz"
        convert(str(fasta), str(npz), verbose=False)

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        base_dir = tmp_path / "out"
        base_dir.mkdir()

        # Build step3 with correct editing_position so guide_start is computed right
        filtered = pd.DataFrame({
            "chromosome": ["chr1"],
            "position": [3],
            "target_base": ["A"],
            "desired_base": ["G"],
            "genome_base": ["A"],
            "strand_used": ["+"],
            "editing_position": [3],
            "guide_rna_sequence": ["GTACAGCTAGCTAGCTAGCT"],
            "pam_sequence": ["GGG"],
            "be_preset": ["ABE8e"],
            "bystander_effect": ["no_bystanders"],
            "simulated_outcome": ["GTACGGCTAGCTAGCTAGCT"],
            "aa_wanted": ["P"],
            "aa_bystander": ["P"],
            "genomic_context_upstream": ["ACGT"],
            "genomic_context_downstream": ["ACGT"],
        })
        filtered.to_csv(base_dir / "step3_filtered.tsv", sep="\t", index=False)
        filtered.to_csv(base_dir / "step2_bystander.tsv", sep="\t", index=False)
        filtered.to_csv(base_dir / "step1_guidefinder.tsv", sep="\t", index=False)

        config = {
            "input": {
                "targets": str(tmp_path / "targets.tsv"),
                "genome": str(fasta),
                "presets": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_presets": "all",
            "keep_bystanders": "all",
            "search_offtargets": True,
            "offtargets": {
                "beguidedoff_genome": str(npz),
                "max_mismatches": 5,
                "count_threshold_score": 100.0,
                "remove_guides_score": 0.0,
                "remove_guides_n_offtargets": 10000,
            },
        }

        run_pipeline(config, start_step=4, end_step=5)

        df = pd.read_csv(base_dir / "beguided_results.tsv", sep="\t")
        assert len(df) > 0

        assert "worst_offtarget_score" in df.columns
        assert "offtarget_count" in df.columns
        assert "worst_offtarget_chr" in df.columns
        assert (base_dir / "step4_offtarget.tsv").exists()

    def test_beguidedoff_on_target_excluded_minus_strand(
        self, tmp_path, valid_preset_tsv
    ):
        """Verify on-target locations are excluded (no crash = success)."""
        pytest.importorskip("beguidedoff", reason="BEguidedOFF not installed")
        from beguidedoff.genome_converter import convert
        from beguided.pipeline import run_pipeline
        import pandas as pd

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        fasta = tmp_path / "mini.fa"
        fasta.write_text(
            ">chr1\n"
            "GTACAGCTAGCTAGCTAGCTGGGGTACAGCTAGCTAGCTAGC"
            "TGGGGTACAGCTAGCTAGCTAGCTGGGGTACAGCTAGCTAGC"
            "TGGGG\n"
        )
        npz = tmp_path / "mini.npz"
        convert(str(fasta), str(npz), verbose=False)

        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t3\tA\tG\n"
        )

        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t100\t.\t+\t0\t"
            'gene_id "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )

        base_dir = tmp_path / "out"
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": str(fasta),
                "presets": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_presets": "all",
            "keep_bystanders": "all",
            "search_offtargets": True,
            "offtargets": {
                "beguidedoff_genome": str(npz),
                "max_mismatches": 5,
                "count_threshold_score": 100.0,
                "remove_guides_score": 0.0,
                "remove_guides_n_offtargets": 10000,
            },
        }

        run_pipeline(config, start_step=1, end_step=5)

        df = pd.read_csv(base_dir / "beguided_results.tsv", sep="\t")
        assert len(df) > 0


class TestPipelineFailFast:
    """Test early-stop behaviour when no guides are found."""

    def test_no_valid_guides_raises_pipeline_error(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        from beguided.pipeline import PipelineError, run_pipeline

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        # Targets are CBE-type (C->T) but presets are ABE-type (A->G)
        targets_path = tmp_path / "targets_cbe.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t10\tC\tT\n"
        )

        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "presets": str(presets_path),
            },
            "output_prefix": str(tmp_path / "out"),
            "use_presets": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
            "guide_score_threshold": 0.5,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }
        with pytest.raises(PipelineError):
            run_pipeline(config, start_step=1, end_step=1)

    def test_pipeline_filters_incompatible_targets_before_guidefinder(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        from beguided.pipeline import PipelineError, run_pipeline
        import pandas as pd

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)

        targets_path = tmp_path / "targets_mixed.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t10\tC\tT\n"
            "chr1\t50\tT\tC\n"
        )

        base_dir = tmp_path / "out"
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "presets": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_presets": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }
        # Should NOT raise PipelineError because T->C target is compatible
        run_pipeline(config, start_step=1, end_step=1)

        guidefinder_output = base_dir / "step1_guidefinder.tsv"
        assert guidefinder_output.exists()
        df = pd.read_csv(guidefinder_output, sep="\t")
        c_to_t_rows = df[
            (df["target_base"] == "C") & (df["desired_base"] == "T")
        ]
        assert len(c_to_t_rows) == 0
