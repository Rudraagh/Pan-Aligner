"""GUI wiring that must not touch tracked docs, and GUI commands for the forwarded main.py modes."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT, ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import gui_app
import main as project_main
import proposed_methodology
import shared_query_eval


def _bare_gui():
    gui = gui_app.PanAlignerGUI.__new__(gui_app.PanAlignerGUI)
    gui.cheatsheet_metrics_text, gui.cheatsheet_local_text, gui.cheatsheet_audit_text = "metrics", "local", "audit"
    gui._load_json = lambda path: None
    return gui


class CoverageAuditTests(unittest.TestCase):
    def test_refresh_shows_the_audit_file_without_rewriting_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            audit = Path(tmp) / "paper_coverage_audit.md"
            audit.write_text("Current audit text\n", encoding="utf-8")
            before = (audit.read_bytes(), audit.stat().st_mtime_ns)
            gui, shown = _bare_gui(), {}
            gui._set_text_widget = lambda widget, content: shown.__setitem__(widget, content)
            with mock.patch.object(gui_app, "PAPER_COVERAGE_REPORT", audit):
                gui._refresh_cheatsheet()
            self.assertEqual((audit.read_bytes(), audit.stat().st_mtime_ns), before)
            self.assertEqual(shown["audit"], "Current audit text\n")

    def test_missing_audit_is_reported_and_not_created(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            audit = Path(tmp) / "paper_coverage_audit.md"
            gui, shown = _bare_gui(), {}
            gui._set_text_widget = lambda widget, content: shown.__setitem__(widget, content)
            with mock.patch.object(gui_app, "PAPER_COVERAGE_REPORT", audit), \
                    mock.patch.object(gui_app, "messagebox") as box, mock.patch.object(gui_app, "open_path") as opener:
                gui._refresh_cheatsheet()
                gui._open_paper_coverage_report()
            self.assertFalse(audit.exists())
            self.assertIn("not found", shown["audit"])
            opener.assert_not_called()
            box.showinfo.assert_called_once()

    def test_gui_has_no_embedded_copy_of_the_audit(self) -> None:
        source = (ROOT / "gui_app.py").read_text(encoding="utf-8")
        self.assertNotIn("PBSIM2 and correctness evaluation with paftools are not implemented", source)
        self.assertFalse(hasattr(gui_app.PanAlignerGUI, "_write_paper_coverage_report"))


class LocalReproductionSummaryTests(unittest.TestCase):
    READINESS = {"paftools_runnable": True, "exact_paper_style_paftools_correctness_eval_possible": False,
                 "why_not_possible": "Current dataset headers do not encode truth interval coordinates."}

    def test_points_to_the_truth_evaluation_when_it_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = Path(tmp) / "REPORT.md"
            report.write_text("# report\n", encoding="utf-8")
            with mock.patch.object(gui_app, "TRUTH_EVAL_REPORT", report):
                text = _bare_gui()._build_local_reproduction_summary(None, self.READINESS)
        self.assertIn("Truth-labelled evaluation (locus scale)", text)
        self.assertIn(str(report), text)
        self.assertIn("labelled FASTAs", text)  # readiness now scoped to the original inputs
        self.assertNotIn("Why exact paftools scoring is left out", text)

    def test_says_when_the_truth_evaluation_is_missing(self) -> None:
        with mock.patch.object(gui_app, "TRUTH_EVAL_REPORT", Path(tempfile.gettempdir()) / "no_such_truth_report.md"):
            text = _bare_gui()._build_local_reproduction_summary(None, self.READINESS)
        self.assertIn("not run yet", text)


class ForwardedModeCommandTests(unittest.TestCase):
    def test_commands_put_the_forwarded_mode_first(self) -> None:
        expected = {
            "proposed-final": ["--proposed", "--final-proposed-evaluation", "--backend", "gpu"],
            "proposed-demo": ["--proposed", "--pg-scunk-stress-demo", "--backend", "gpu"],
            "shared-query": ["--shared-query-eval", "--backend", "gpu"],
        }
        for mode, arguments in expected.items():
            with self.subTest(mode=mode):
                self.assertEqual(gui_app.build_forwarded_command(mode, "gpu", python="py"), ["py", str(gui_app.MAIN_SCRIPT), *arguments])

    def test_invalid_backend_or_mode_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            gui_app.build_forwarded_command("proposed-final", "tpu")
        with self.assertRaises(ValueError):
            gui_app.build_forwarded_command("truth-eval", "cpu")

    def test_every_gui_command_reaches_its_tool_through_main(self) -> None:
        targets = {"proposed-final": proposed_methodology, "proposed-demo": proposed_methodology, "shared-query": shared_query_eval}
        for mode, module in targets.items():
            command = gui_app.build_forwarded_command(mode, "cpu")
            with self.subTest(mode=mode), mock.patch.object(module, "main") as target:
                project_main.main(command[2:])
                target.assert_called_once_with(command[3:])

    def test_gui_mode_list_offers_the_forwarded_modes(self) -> None:
        values = {value for _, value in gui_app.PIPELINE_MODES}
        self.assertLessEqual(set(gui_app.FORWARDED_GUI_MODES), values)
        self.assertLessEqual({"--full-pipeline", "--theory-only", "--evaluate"}, values)


if __name__ == "__main__":
    unittest.main()
