from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.models import AssemblyInput
from proposed.quality_feedback import run_pg_scunk_feedback_loop


class ProposedQualityFeedbackTests(unittest.TestCase):
    def _assemblies(self) -> list[AssemblyInput]:
        return [
            AssemblyInput("hap1", "TOY", "HEALTHY", "ACGTACGTACGT"),
            AssemblyInput("hap2", "TOY", "HEALTHY", "ACGTACGTTCGT"),
            AssemblyInput("hap3", "TOY", "UNHEALTHY", "ACGTTCGTACGT"),
            AssemblyInput("hap4", "TOY", "UNHEALTHY", "ACGTTCGTTCGT"),
        ]

    def test_feedback_loop_returns_round_summaries_and_best_graph(self) -> None:
        result = run_pg_scunk_feedback_loop(
            self._assemblies(),
            base_parameters={"k": 5, "window_size": 4, "min_node_support": 3, "max_anchors": 3, "min_anchor_spacing": 2.0},
            backend="cpu",
            quality_threshold=0.90,
            max_rounds=2,
        )
        self.assertGreaterEqual(len(result["round_summaries"]), 1)
        self.assertIn("overall_quality_score", result["best_quality"])
        self.assertIsNotNone(result["graph"])


if __name__ == "__main__":
    unittest.main()
