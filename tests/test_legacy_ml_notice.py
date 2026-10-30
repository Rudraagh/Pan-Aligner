"""The legacy healthy/unhealthy ML path must say it is not a valid classifier evaluation."""

from __future__ import annotations

import contextlib
import importlib
import io
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

LEGACY_MODULES = ("pipeline", "extract_features", "train_ml", "predict_ml", "evaluate")


def _import_or_skip(test: unittest.TestCase, name: str):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as exc:  # scikit-learn is only declared for the legacy path
        test.skipTest(f"{name} needs {exc.name}, which is not installed in this environment")


class LegacyMlNoticeTests(unittest.TestCase):
    def test_module_docstrings_mark_them_legacy(self) -> None:
        for name in LEGACY_MODULES:
            with self.subTest(module=name):
                module = _import_or_skip(self, name)
                self.assertTrue((module.__doc__ or "").startswith("LEGACY"), module.__doc__)

    def test_entry_points_print_the_notice_before_doing_anything(self) -> None:
        from common import LEGACY_ML_NOTICE

        for name in LEGACY_MODULES:
            with self.subTest(module=name):
                module = _import_or_skip(self, name)
                stderr, stdout = io.StringIO(), io.StringIO()
                with mock.patch.object(sys, "argv", [f"{name}.py", "--help"]), contextlib.redirect_stderr(stderr), \
                        contextlib.redirect_stdout(stdout), self.assertRaises(SystemExit):
                    module.main()
                self.assertIn(LEGACY_ML_NOTICE, stderr.getvalue())

    def test_notice_names_the_leakage_and_the_confound(self) -> None:
        from common import LEGACY_ML_NOTICE

        for phrase in ("not a valid classifier", "leak", "genomic position", "docs/13"):
            with self.subTest(phrase=phrase):
                self.assertIn(phrase, LEGACY_ML_NOTICE)


if __name__ == "__main__":
    unittest.main()
