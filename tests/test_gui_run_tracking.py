import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import gui_app


class GuiRunTrackingTests(unittest.TestCase):
    def test_successful_run_records_reported_artifacts_without_file_ordering(self):
        output_dir = ROOT / "outputs" / "alignments" / "custom_query_real_psen2_presentation_test"
        prediction = output_dir / "prediction.json"
        plot = output_dir / "custom_query_scores.png"
        if not prediction.exists() or not plot.exists():
            self.skipTest("Custom-query demonstration artifacts are not present.")

        tracker = gui_app.PanAlignerGUI.__new__(gui_app.PanAlignerGUI)
        tracker._run_output_snapshot = {}
        after = {str(path.resolve()): (1, path.stat().st_size) for path in output_dir.iterdir() if path.is_file()}
        tracker._snapshot_output_files = lambda: after
        manifest_path = ROOT / "tests" / "_tmp" / "latest_gui_run_test.json"

        summary = {
            "mode": "custom-query-analysis",
            "selected_gene": "PSEN2",
            "analysis_status": "ALIGNMENT_DETECTED",
            "prediction_json": str(prediction),
            "score_plot": str(plot),
        }
        with patch.object(gui_app, "LATEST_GUI_RUN_PATH", manifest_path), patch.object(Path, "write_text") as write_text:
            tracker._record_completed_run(summary)

        record = json.loads(write_text.call_args.args[0])
        self.assertEqual(record["mode"], "custom-query-analysis")
        self.assertEqual(record["selected_gene"], "PSEN2")
        self.assertEqual(record["status"], "ALIGNMENT_DETECTED")
        self.assertIn(str(plot.resolve()), record["artifacts"])
        self.assertIn(str(prediction.resolve()), record["artifacts"])
