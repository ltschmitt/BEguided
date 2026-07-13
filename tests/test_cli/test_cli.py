"""Tests for CLI entry points."""

import pytest


class TestCLIImports:
    """Ensure CLI modules import correctly."""

    def test_import_main(self):
        from beguided.cli import main
        assert callable(main)

    def test_import_web_main(self):
        from beguided.web.cli import main
        assert callable(main)
