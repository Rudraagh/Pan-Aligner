from __future__ import annotations

import json
import shutil
import sys
import unittest
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.stress_experiment import run_pg_scunk_stress_experiment


class PGScunkStressExperimentTests(unittest.TestCase):
    def test_demo_is_deterministic_and_exercises_real_feedback(self) -> None:
        first_dir = ROOT / "tests" / "_tmp" / f"step5_demo_{uuid.uuid4().hex}_first"
        second_dir = ROOT / "tests" / "_tmp" / f"step5_demo_{uuid.uuid4().hex}_second"
        try:
            first = run_pg_scunk_stress_experiment(output_dir=first_dir)
            second = run_pg_scunk_stress_experiment(output_dir=second_dir)
            first_metadata = json.loads(Path(first["metadata_json"]).read_text(encoding="utf-8"))
            second_metadata = json.loads(Path(second["metadata_json"]).read_text(encoding="utf-8"))
            self.assertTrue(first_metadata["initial_below_target"])
            self.assertGreaterEqual(first_metadata["unique_configurations_evaluated"], 2)
            self.assertGreater(first_metadata["best_quality"], first_metadata["initial_quality"])
            self.assertEqual(first_metadata["best_quality"], second_metadata["best_quality"])
            self.assertEqual(first_metadata["best_configuration"], second_metadata["best_configuration"])
            self.assertTrue(Path(first["pg_scunk_iterations_csv"]).exists())
            self.assertTrue(Path(first["pg_scunk_summary_json"]).exists())
            self.assertTrue(Path(first["pg_scunk_quality_plot"]).exists())
        finally:
            shutil.rmtree(first_dir, ignore_errors=True)
            shutil.rmtree(second_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
