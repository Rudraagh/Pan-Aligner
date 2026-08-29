from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.minimizer_engine import compute_minimizers


class MinimizerEngineTests(unittest.TestCase):
    def test_minimizer_engine_returns_deterministic_records(self) -> None:
        payload = compute_minimizers("ACGTACGTACGT", k=3, window_size=4, assembly_id="hap1", backend="cpu")
        self.assertEqual(payload["backend_used"], "cpu")
        self.assertGreaterEqual(payload["minimizer_count"], 1)
        first = payload["records"][0]
        self.assertEqual(first["source_assembly"], "hap1")
        self.assertEqual(len(first["minimizer"]), 3)

    def test_gpu_request_falls_back_honestly_when_backend_missing(self) -> None:
        payload = compute_minimizers("ACGTACGTACGT", k=3, window_size=4, assembly_id="hap1", backend="gpu")
        self.assertIn(payload["backend_used"], {"cpu", "gpu"})
        if payload["backend_used"] == "cpu":
            self.assertIn("fell back", payload["backend_note"])


if __name__ == "__main__":
    unittest.main()
