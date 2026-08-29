from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.haplotype_scoring import pairwise_haplotype_consistency
from proposed.minimizer_engine import compute_minimizers
from proposed.models import MinimizerRecord


def _records(sequence: str, assembly_id: str) -> list[MinimizerRecord]:
    payload = compute_minimizers(sequence, k=3, window_size=3, assembly_id=assembly_id, backend="cpu")
    return [MinimizerRecord(**item) for item in payload["records"]]


class HaplotypeScoringTests(unittest.TestCase):
    def test_related_sequences_score_higher_than_divergent_sequences(self) -> None:
        similar = pairwise_haplotype_consistency(_records("ACGTACGTAC", "a"), _records("ACGTACGTTC", "b"))
        divergent = pairwise_haplotype_consistency(_records("ACGTACGTAC", "a"), _records("TTTTGGGGCC", "c"))
        self.assertGreater(similar, divergent)
        self.assertGreaterEqual(similar, 0.0)
        self.assertLessEqual(similar, 1.0)


if __name__ == "__main__":
    unittest.main()
