"""Known-answer checks of how the reported metrics are defined.

These pin down properties the evaluation reports rely on, including the ones
that make some proposed consistency terms constant by construction.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from parse_gaf import AlignmentResult
from proposed.graph_constructor import build_proposed_pangenome_graph
from proposed.models import AssemblyInput


def alignment(query_start, query_end, block, matches=90, query_length=100):
    return AlignmentResult("q", query_length, query_start, query_end, "+", ">s1", 1000, 0, 100, matches, block, 60)


class AlignmentMetricTests(unittest.TestCase):
    def test_coverage_is_the_aligned_query_fraction_even_with_deletions(self) -> None:
        # 100 bp query fully aligned with a 5 bp deletion: block length 105.
        self.assertEqual(alignment(0, 100, 105).coverage, 1.0)
        self.assertEqual(alignment(10, 60, 52).coverage, 0.5)
        self.assertEqual(alignment(0, 0, 0, query_length=0).coverage, 0.0)

    def test_identity_is_matches_over_alignment_columns(self) -> None:
        self.assertAlmostEqual(alignment(0, 100, 105, matches=95).identity, 95 / 105)


class ProposedConsistencyDefinitionTests(unittest.TestCase):
    """Edges are created from observed paths, so two transition terms are always 1."""

    def setUp(self) -> None:
        self.graph = build_proposed_pangenome_graph(
            [
                AssemblyInput("a", "G", "HEALTHY", "ACGTTGCATGACCTAGGATCCAGTTACG"),
                AssemblyInput("b", "G", "HEALTHY", "ACGTTGCATGACGTAGGATCCAGTTACG"),
                AssemblyInput("c", "G", "UNHEALTHY", "TTTTACGTTGCATGACCTAGGATCC"),
            ],
            k=4, window_size=3, backend="cpu",
        )

    def test_path_transition_fraction_is_always_one(self) -> None:
        self.assertTrue(all(path.transition_fraction == 1.0 for path in self.graph.paths))

    def test_node_local_path_consistency_is_always_one(self) -> None:
        self.assertTrue(all(node.local_path_consistency == 1.0 for node in self.graph.nodes))

    def test_scores_and_ratios_stay_in_unit_interval(self) -> None:
        quality = self.graph.quality
        for key in ("mean_node_consistency", "mean_edge_consistency", "mean_path_consistency", "graph_path_consistency",
                    "mean_node_support_ratio", "anchor_coverage_ratio", "branching_node_ratio", "overall_quality_score"):
            with self.subTest(metric=key):
                self.assertGreaterEqual(quality[key], 0.0)
                self.assertLessEqual(quality[key], 1.0)
        # path score = 0.40*1 + 0.30*node + 0.30*edge, so it can never fall below 0.40
        self.assertTrue(all(path.consistency_score >= 0.40 for path in self.graph.paths))


if __name__ == "__main__":
    unittest.main()
