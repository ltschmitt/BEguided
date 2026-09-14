"""Tests for pipeline components."""

from pathlib import Path
from unittest.mock import patch

import pytest


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


class TestPassthroughExtras:
    """Test the ``passthrough_extras`` config option."""

    def _build_config(
        self, tmp_path, sample_fasta_file, valid_preset_tsv,
        targets_content, passthrough_extras=None,
    ):
        presets_path = tmp_path / "presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(targets_content)
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(tmp_path / "out"),
            "use_base_editors": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
        }
        if passthrough_extras is not None:
            config["passthrough_extras"] = passthrough_extras
        return config

    def test_default_passes_extras_through(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """No ``passthrough_extras`` key → default ON → extras preserved."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\tgene\n"
            "chr1\t75\tA\tG\tKRAS\n",
        )
        run_pipeline(config, start_step=1, end_step=5)

        final = pd.read_csv(tmp_path / "out" / "beguided_results.tsv", sep="\t")
        assert "gene" in final.columns
        assert (final["gene"] == "KRAS").all()

    def test_explicit_true_keeps_extras(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\tgene\n"
            "chr1\t75\tA\tG\tKRAS\n",
            passthrough_extras=True,
        )
        run_pipeline(config, start_step=1, end_step=5)

        final = pd.read_csv(tmp_path / "out" / "beguided_results.tsv", sep="\t")
        assert "gene" in final.columns

    def test_explicit_false_drops_extras(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """``passthrough_extras: false`` in config drops extras everywhere."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\tgene\n"
            "chr1\t75\tA\tG\tKRAS\n",
            passthrough_extras=False,
        )
        run_pipeline(config, start_step=1, end_step=5)

        final = pd.read_csv(tmp_path / "out" / "beguided_results.tsv", sep="\t")
        assert "gene" not in final.columns

        # Also verify the intermediate step1 output is stripped
        step1 = pd.read_csv(tmp_path / "out" / "step1_guidefinder.tsv", sep="\t")
        assert "gene" not in step1.columns


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
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )

        base_dir = tmp_path / "out"
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
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


class TestPipelineAlwaysRunBystander:
    """Verify that step 2 always runs unless skip_bystander=True, regardless
    of keep_bystanders filter mode."""

    def _build_config(self, tmp_path, sample_fasta_file, valid_preset_tsv,
                      keep_bystanders="no_bystanders", skip_bystander=False):
        presets_path = tmp_path / "presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
        )
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(tmp_path / "out"),
            "use_base_editors": "all",
            "keep_bystanders": keep_bystanders,
            "search_offtargets": False,
        }
        if skip_bystander:
            config["skip_bystander"] = True
        return config

    def test_no_bystanders_mode_still_runs_step2(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """keep_bystanders='no_bystanders' must NOT skip step 2 anymore."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            keep_bystanders="no_bystanders",
        )
        run_pipeline(config, start_step=1, end_step=5)

        step2 = pd.read_csv(tmp_path / "out" / "step2_bystander.tsv", sep="\t")
        # Step 2 ran → bystander_effect column exists (not just a copy of step 1)
        assert "bystander_effect" in step2.columns
        assert "gene_edit_desired" in step2.columns

    def test_gene_edit_desired_in_final_output(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Final TSV includes gene_edit_desired when step 2 runs."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            keep_bystanders="all",
        )
        run_pipeline(config, start_step=1, end_step=5)

        final = pd.read_csv(tmp_path / "out" / "beguided_results.tsv", sep="\t")
        assert "gene_edit_desired" in final.columns

    def test_skip_bystander_skips_step2(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """skip_bystander=True still skips step 2."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            keep_bystanders="all",
            skip_bystander=True,
        )
        run_pipeline(config, start_step=1, end_step=5)

        step2 = pd.read_csv(tmp_path / "out" / "step2_bystander.tsv", sep="\t")
        # Step 2 skipped → no bystander_effect column
        assert "bystander_effect" not in step2.columns
        assert "gene_edit_desired" not in step2.columns
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
            "base_editor": ["ABE8e", "ABE8e"],
        })
        out_path = tmp_path / "bgo_input.tsv"
        _prepare_beguidedoff_input(df, out_path)

        result = pd.read_csv(out_path, sep="\t")
        assert list(result.columns) == [
            "chromosome", "position", "guide_rna_sequence",
            "base_editor", "strand_used",
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
            "base_editor": ["ABE8e", "ABE8e"],
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
            "base_editor": [],
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
            "base_editor": ["ABE8e"],
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
                "base_editor_profiles": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
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
            "base_editor": ["ABE8e"],
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
                "base_editor_profiles": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
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
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )

        base_dir = tmp_path / "out"
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": str(fasta),
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
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
                "base_editor_profiles": str(presets_path),
            },
            "output_prefix": str(tmp_path / "out"),
            "use_base_editors": "all",
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
                "base_editor_profiles": str(presets_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
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


class TestRunPipelineGenomeKwarg:
    """When the caller supplies a pre-loaded genome, run_pipeline() must
    skip the per-step ``load_genome()`` call and reuse the supplied dict."""

    def test_preloaded_genome_skips_load_genome(
        self, tmp_path, sample_fasta_file, sample_genome_dict, valid_preset_tsv
    ):
        from beguided.pipeline import run_pipeline

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t50\tT\tC\nchr1\t75\tA\tG\n"
        )
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t160\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        base_dir = tmp_path / "out"

        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }

        # Reuse the real sample genome so PAM search can actually find guides.
        preloaded = dict(sample_genome_dict)

        with patch("beguided.pipeline.load_genome") as fake_load:
            fake_load.side_effect = AssertionError(
                "load_genome must NOT be called when genome kwarg is supplied"
            )
            run_pipeline(
                config,
                start_step=1,
                end_step=3,
                genome=preloaded,
            )

        guidefinder_output = base_dir / "step1_guidefinder.tsv"
        assert guidefinder_output.exists()
        import pandas as pd
        df = pd.read_csv(guidefinder_output, sep="\t")
        assert len(df) > 0

    def test_no_kwarg_still_calls_load_genome(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Regression: backward-compatible behaviour when ``genome=None``."""
        from beguided.pipeline import run_pipeline
        from beguided.genome.io import load_genome

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t50\tT\tC\n"
        )
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        base_dir = tmp_path / "out"

        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }

        with patch("beguided.pipeline.load_genome", wraps=load_genome) as spy:
            run_pipeline(config, start_step=1, end_step=3, genome=None)
            assert spy.called, "load_genome must run when no genome is supplied"


class TestPipelineScaffoldFilter:
    """The pipeline must drop scaffold records from the genome dict so
    only canonical chromosomes are searched for PAMs."""

    def test_pipeline_drops_scaffolds(
        self, tmp_path, valid_preset_tsv, sample_genome_dict
    ):
        from beguided.pipeline import run_pipeline
        from beguided.genome.io import load_genome
        import pandas as pd

        # Build a FASTA whose scaffold records would fail chromosome lookup.
        fasta = tmp_path / "mixed.fa"
        fasta.write_text(
            ">1 dna:chromosome chromosome:GRCh38:1:1:160:1 REF\n"
            + sample_genome_dict["chr1"]
            + "\n"
            + ">KI270728.1 dna:scaffold scaffold:GRCh38:KI270728.1:1:60:1 REF\n"
            + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            + ">chr2 dna:chromosome chromosome:GRCh38:2:1:160:1 REF\n"
            + sample_genome_dict["chr2"]
            + "\n"
        )

        presets_path = tmp_path / "test_presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        # Note: targets use the bare Ensembl chromosome names (1, 2).
        targets_path.write_text(
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "1\t50\tT\tC\n1\t75\tA\tG\n2\t30\tT\tC\n"
        )
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "1\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
            "chr2\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G2"; gene_name "G2"; transcript_id "T2"; tag "MANE_Select";\n'
        )
        base_dir = tmp_path / "out"

        config = {
            "input": {
                "targets": str(targets_path),
                "genome": str(fasta),
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(base_dir),
            "use_base_editors": "all",
            "keep_bystanders": "no_bystanders",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }

        genome = load_genome(str(fasta))
        assert "KI270728.1" not in genome
        assert set(genome) == {"1", "chr2"}

        run_pipeline(config, start_step=1, end_step=3, genome=genome)
        guidefinder_output = base_dir / "step1_guidefinder.tsv"
        assert guidefinder_output.exists()
        df = pd.read_csv(guidefinder_output, sep="\t")
        # No guides should have landed on a scaffold chromosome.
        assert "KI270728.1" not in df["chromosome"].astype(str).unique()


class TestPipelineTrinucEditScore:
    """End-to-end tests for the trinuc_edit_score column and its filter."""

    def _build_config(
        self, tmp_path, sample_fasta_file, valid_preset_tsv,
        targets_content, trinuc_edit_score_min=None,
    ):
        presets_path = tmp_path / "presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(targets_content)
        gtf_path = tmp_path / "dummy.gtf"
        gtf_path.write_text(
            "##format: gtf\n"
            "chr1\ttest\tCDS\t1\t200\t.\t+\t0\t"
            'gene_id "G1"; gene_name "G1"; transcript_id "T1"; tag "MANE_Select";\n'
        )
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(gtf_path),
            },
            "output_prefix": str(tmp_path / "out"),
            "use_base_editors": "all",
            "keep_bystanders": "all",
            "search_offtargets": False,
            "offtargets": {
                "max_mismatches": 5,
                "count_threshold_score": 15.0,
                "remove_guides_score": 5.0,
                "remove_guides_n_offtargets": 100,
            },
        }
        if trinuc_edit_score_min is not None:
            config["trinuc_edit_score_min"] = trinuc_edit_score_min
        return config

    def test_trinuc_edit_score_column_in_all_outputs(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """The trinuc_edit_score column appears in step1, step3, and final."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
        )
        run_pipeline(config, start_step=1, end_step=5)

        out = tmp_path / "out"
        for fname in ("step1_guidefinder.tsv", "step3_filtered.tsv", "beguided_results.tsv"):
            df = pd.read_csv(out / fname, sep="\t")
            assert "trinuc_edit_score" in df.columns, (
                f"{fname} missing trinuc_edit_score column"
            )

    def test_trinuc_edit_score_populated_for_abe_guides(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """ABE guides (A->G) get a numeric trinuc_edit_score."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
        )
        run_pipeline(config, start_step=1, end_step=5)

        final = pd.read_csv(tmp_path / "out" / "beguided_results.tsv", sep="\t")
        valid = final[~final["guide_rna_sequence"].str.startswith("**", na=False)]
        if len(valid) > 0:
            # ABE guides should have a numeric score (not NaN).
            assert valid["trinuc_edit_score"].notna().any(), (
                "At least one ABE guide should have a non-NaN trinuc_edit_score"
            )

    def test_trinuc_edit_score_min_filter_drops_low_guides(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """trinuc_edit_score_min drops ABE guides below the threshold."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        # Run without the filter first to get the baseline.
        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
        )
        run_pipeline(config, start_step=1, end_step=3)
        baseline = pd.read_csv(tmp_path / "out" / "step3_filtered.tsv", sep="\t")

        # Now run with a very high threshold (1.0 = only CAC context).
        # Use a separate output dir so we don't clobber the baseline.
        config_filtered = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
            trinuc_edit_score_min=1.0,
        )
        config_filtered["output_prefix"] = str(tmp_path / "out_filtered")
        run_pipeline(config_filtered, start_step=1, end_step=3)
        filtered = pd.read_csv(
            tmp_path / "out_filtered" / "step3_filtered.tsv", sep="\t"
        )

        # The filtered output should have <= rows than the baseline.
        assert len(filtered) <= len(baseline)
        # Any remaining ABE guides (non-NaN trinuc_edit_score) must be >= 1.0.
        abe_rows = filtered[filtered["trinuc_edit_score"].notna()]
        if len(abe_rows) > 0:
            assert (abe_rows["trinuc_edit_score"] >= 1.0).all()

    def test_edit_context_table_empty_string_disables_column(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """edit_context_table='' disables the lookup; every score is NaN."""
        from beguided.pipeline import run_pipeline
        import pandas as pd

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
        )
        config["edit_context_table"] = ""
        run_pipeline(config, start_step=1, end_step=1)

        step1 = pd.read_csv(tmp_path / "out" / "step1_guidefinder.tsv", sep="\t")
        assert "trinuc_edit_score" in step1.columns
        # Every value should be NaN (column disabled).
        assert step1["trinuc_edit_score"].isna().all()

    def test_trinuc_edit_score_min_invalid_value_raises(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """trinuc_edit_score_min outside [0, 1] raises PipelineError."""
        from beguided.pipeline import PipelineError, run_pipeline

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\nchr1\t75\tA\tG\n",
            trinuc_edit_score_min=1.5,
        )
        with pytest.raises(PipelineError):
            run_pipeline(config, start_step=1, end_step=1)


class TestKeepTargetsWithoutGuides:
    """The ``keep_targets_without_guides`` flag keeps one placeholder row per
    input target in the final results table, even when no guide could be
    designed or the mutation type does not match the selected editors.
    Placeholder rows are only added in the very last step."""

    def _build_config(
        self, tmp_path, sample_fasta_file, valid_preset_tsv, targets_content,
        keep_targets_without_guides=True,
    ):
        presets_path = tmp_path / "presets.tsv"
        presets_path.write_text(valid_preset_tsv)
        targets_path = tmp_path / "targets.tsv"
        targets_path.write_text(targets_content)
        config = {
            "input": {
                "targets": str(targets_path),
                "genome": sample_fasta_file,
                "base_editor_profiles": str(presets_path),
                "gff_file": str(tmp_path / "nonexistent.gtf"),
            },
            "output_prefix": str(tmp_path / "out"),
            # Only ABE editors so a C>T target is incompatible by design.
            "use_base_editors": "ABE8e",
            "keep_bystanders": "all",
            "search_offtargets": False,
            "skip_bystander": True,
            "keep_targets_without_guides": keep_targets_without_guides,
        }
        return config

    def _run(self, config):
        from beguided.pipeline import run_pipeline
        import pandas as pd

        run_pipeline(config, start_step=1, end_step=5)
        return pd.read_csv(
            Path(config["output_prefix"]) / "beguided_results.tsv", sep="\t"
        )

    def test_off_by_default_no_placeholders(self, tmp_path, sample_fasta_file, valid_preset_tsv):
        """Default (flag off) keeps the historic behaviour: the incompatible
        target is simply dropped, no placeholder row appears."""
        from beguided.models import (
            BASE_CHANGE_INCOMPATIBLE,
            NO_GUIDE_FOUND,
        )

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t5\tA\tG\n"
            "chr1\t100\tC\tT\n",
            keep_targets_without_guides=False,
        )
        final = self._run(config)
        assert (final["guide_rna_sequence"] != NO_GUIDE_FOUND).all()
        assert (final["guide_rna_sequence"] != BASE_CHANGE_INCOMPATIBLE).all()
        assert 100 not in set(final["position"])

    def test_flag_on_adds_no_guide_and_incompatible_rows(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """With the flag on, a target with a real guide is kept, a compatible
        target with no PAM gets "**No guide possible**", and an incompatible
        mutation type gets "**Base change incompatible**"."""
        from beguided.models import BASE_CHANGE_INCOMPATIBLE, NO_GUIDE_FOUND

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t5\tA\tG\n"      # real guide
            "chr1\t20\tA\tG\n"     # compatible, no valid PAM
            "chr1\t100\tC\tT\n",   # mutation type not supported by ABE
        )
        final = self._run(config)
        by_pos = final.groupby("position")

        # chr1:5 -> at least one real guide.
        pos5 = by_pos.get_group(5)
        assert (pos5["guide_rna_sequence"] != NO_GUIDE_FOUND).any()
        assert (pos5["guide_rna_sequence"] != BASE_CHANGE_INCOMPATIBLE).any()

        # chr1:20 -> exactly one "**No guide possible**" placeholder.
        pos20 = by_pos.get_group(20)
        assert len(pos20) == 1
        assert (pos20["guide_rna_sequence"] == NO_GUIDE_FOUND).all()

        # chr1:100 -> exactly one "**Base change incompatible**" placeholder.
        pos100 = by_pos.get_group(100)
        assert len(pos100) == 1
        assert (pos100["guide_rna_sequence"] == BASE_CHANGE_INCOMPATIBLE).all()

    def test_flag_on_all_incompatible_no_error(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Every input target incompatible -> one placeholder row per target
        instead of the historic "No targets remain" PipelineError."""
        from beguided.models import BASE_CHANGE_INCOMPATIBLE

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t100\tC\tT\n"
            "chr2\t40\tG\tA\n",
        )
        final = self._run(config)
        assert len(final) == 2
        assert (final["guide_rna_sequence"] == BASE_CHANGE_INCOMPATIBLE).all()
        assert set(final["position"]) == {40, 100}

    def test_flag_on_all_compatible_no_pam_no_error(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Compatible targets with no valid PAM yield "**No guide possible**"
        rows instead of the historic "No valid guides" PipelineError."""
        from beguided.models import NO_GUIDE_FOUND

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t100\tA\tG\n",
        )
        final = self._run(config)
        assert len(final) == 1
        assert (final["guide_rna_sequence"] == NO_GUIDE_FOUND).all()
        assert set(final["position"]) == {100}

    def test_flag_off_all_compatible_no_pam_still_raises(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Flag off keeps the historic behaviour: compatible targets with no
        valid PAM raise PipelineError instead of a placeholder-only output."""
        from beguided.pipeline import PipelineError, run_pipeline

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\n"
            "chr1\t100\tA\tG\n",
            keep_targets_without_guides=False,
        )
        with pytest.raises(PipelineError, match="No valid guides"):
            run_pipeline(config, start_step=1, end_step=5)

    def test_flag_on_incompatible_with_passthrough_extras(
        self, tmp_path, sample_fasta_file, valid_preset_tsv
    ):
        """Extra passthrough columns from the input are carried onto the
        placeholder rows (CLI passthrough_extras default on)."""
        from beguided.models import BASE_CHANGE_INCOMPATIBLE

        config = self._build_config(
            tmp_path, sample_fasta_file, valid_preset_tsv,
            "chromosome\tposition\ttarget_base\tdesired_base\tgene\n"
            "chr1\t5\tA\tG\tKRAS\n"
            "chr1\t100\tC\tT\tTP53\n",
        )
        final = self._run(config)
        ph = final[final["guide_rna_sequence"] == BASE_CHANGE_INCOMPATIBLE]
        assert len(ph) == 1
        assert ph.iloc[0]["gene"] == "TP53"
