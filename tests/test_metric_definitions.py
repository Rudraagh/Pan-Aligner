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


class ProposedQualityScoreTests(unittest.TestCase):
    """Overall quality uses graph structure only: no anchor term, no min_node_support effect."""

    ASSEMBLIES = [
        AssemblyInput("a", "G", "HEALTHY", "ACGTTGCATGACCTAGGATCCAGTTACG"),
        AssemblyInput("b", "G", "HEALTHY", "ACGTTGCATGACGTAGGATCCAGTTACG"),
        AssemblyInput("c", "G", "UNHEALTHY", "TTTTACGTTGCATGACCTAGGATCC"),
    ]

    def build(self, **overrides):
        return build_proposed_pangenome_graph(self.ASSEMBLIES, **{"k": 4, "window_size": 3, "backend": "cpu", **overrides})

    def test_overall_quality_formula_has_no_anchor_term(self) -> None:
        q = self.build().quality
        expected = (
            0.25 * q["mean_node_consistency"] + 0.20 * q["mean_edge_consistency"] + 0.15 * q["mean_path_consistency"]
            + 0.15 * q["mean_node_support_ratio"] + 0.10 * (1.0 - min(1.0, q["branching_node_ratio"]))
        ) / 0.85
        self.assertAlmostEqual(q["overall_quality_score"], expected, places=12)

    def test_anchor_settings_change_reported_coverage_but_not_quality(self) -> None:
        # Uniqueness-only weights with one anchor pick a private node, covering 1 of 3 haplotypes.
        few = self.build(max_anchors=1, anchor_alpha=1.0, anchor_beta=0.0, anchor_gamma=0.0)
        many = self.build(max_anchors=12, min_anchor_spacing=1.0)
        self.assertAlmostEqual(few.quality["anchor_coverage_ratio"], 1 / 3)
        self.assertEqual(many.quality["anchor_coverage_ratio"], 1.0)
        self.assertEqual(few.quality["overall_quality_score"], many.quality["overall_quality_score"])

    def test_min_node_support_changes_labels_but_not_quality(self) -> None:
        low, high = self.build(min_node_support=1), self.build(min_node_support=3)
        self.assertNotEqual(low.metadata["rescued_node_count"], high.metadata["rescued_node_count"])
        self.assertEqual(low.quality, high.quality)


if __name__ == "__main__":
    unittest.main()
