"""Tests for CLI commands."""

import pytest


class TestCLIImports:
    """Ensure CLI module imports correctly."""

    def test_import_main(self):
        from beguided.cli import main
        assert callable(main)

    def test_import_commands(self):
        from beguided.cli import (
            cmd_guidefinder,
            cmd_bystander,
            cmd_filter,
            cmd_run_pipeline,
        )
        assert all(callable(f) for f in [
            cmd_guidefinder, cmd_bystander, cmd_filter, cmd_run_pipeline
        ])
