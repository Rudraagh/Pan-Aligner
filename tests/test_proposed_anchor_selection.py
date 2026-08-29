from __future__ import annotations

import sys
import shutil
import unittest
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.anchor_selection import (
    compute_anchor_score,
    score_anchor_candidates,
    select_informative_anchors,
    validate_anchor_weights,
    write_anchor_reports,
)
from proposed.graph_constructor import build_proposed_pangenome_graph
from proposed.models import AnchorCandidate, AssemblyInput


class AnchorSelectionTests(unittest.TestCase):
    def _weight_candidates(self) -> list[AnchorCandidate]:
        return [
            AnchorCandidate("a", "AAA", 1, 1.0, 0.50, 0.0, ["h1"], node_id="n_a", uniqueness=0.90, depth=0.20),
            AnchorCandidate("b", "CCC", 4, 4.0, 0.50, 0.0, ["h1", "h2"], node_id="n_b", uniqueness=0.50, depth=0.95),
            AnchorCandidate("c", "GGG", 2, 8.0, 0.95, 0.0, ["h1"], node_id="n_c", uniqueness=0.50, depth=0.50),
        ]

    def test_spacing_filter_prefers_strong_non_overlapping_anchors(self) -> None:
        anchors = [
            AnchorCandidate("a1", "AAA", 4, 2.0, 0.9, 0.92, ["h1", "h2"]),
            AnchorCandidate("a2", "CCC", 4, 2.5, 0.8, 0.88, ["h1", "h2"]),
            AnchorCandidate("a3", "GGG", 3, 8.0, 0.85, 0.87, ["h1", "h2"]),
        ]
        selected = select_informative_anchors(anchors, max_anchors=3, min_spacing=2.0)
        self.assertEqual([anchor.anchor_id for anchor in selected], ["a1", "a3"])

    def test_manual_formula_matches_expected_score(self) -> None:
        score = compute_anchor_score(0.8, 0.6, 0.9, 0.4, 0.3, 0.3)
        self.assertAlmostEqual(score, 0.4 * 0.8 + 0.3 * 0.6 + 0.3 * 0.9)

    def test_each_weight_term_changes_the_ranking(self) -> None:
        uniqueness_first = score_anchor_candidates(self._weight_candidates(), alpha=1.0, beta=0.0, gamma=0.0)
        depth_first = score_anchor_candidates(self._weight_candidates(), alpha=0.0, beta=1.0, gamma=0.0)
        consistency_first = score_anchor_candidates(self._weight_candidates(), alpha=0.0, beta=0.0, gamma=1.0)
        self.assertEqual(max(uniqueness_first, key=lambda item: item.final_score).anchor_id, "a")
        self.assertEqual(max(depth_first, key=lambda item: item.final_score).anchor_id, "b")
        self.assertEqual(max(consistency_first, key=lambda item: item.final_score).anchor_id, "c")

    def test_weights_and_graph_consistency_are_recorded_in_candidates(self) -> None:
        graph = build_proposed_pangenome_graph(
            [
                AssemblyInput("h1", "TOY", "HEALTHY", "ACGTACGTAC"),
                AssemblyInput("h2", "TOY", "UNHEALTHY", "ACGTTCGTAC"),
            ],
            k=3,
            window_size=2,
            min_node_support=2,
            anchor_alpha=0.4,
            anchor_beta=0.3,
            anchor_gamma=0.3,
            backend="cpu",
        )
        node_scores = {node.node_id: node.consistency_score for node in graph.nodes}
        self.assertTrue(graph.anchor_candidates)
        self.assertTrue(all(candidate.consistency_score == node_scores[candidate.node_id] for candidate in graph.anchor_candidates))
        self.assertEqual(graph.metadata["anchor_scoring"]["weights"], {"alpha": 0.4, "beta": 0.3, "gamma": 0.3})
        self.assertTrue(all(candidate.selected for candidate in graph.anchors))

    def test_graph_anchor_selection_is_deterministic(self) -> None:
        assemblies = [
            AssemblyInput("h1", "TOY", "HEALTHY", "ACGTACGTAC"),
            AssemblyInput("h2", "TOY", "UNHEALTHY", "ACGTTCGTAC"),
        ]
        first = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        second = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        self.assertEqual(
            [(anchor.node_id, anchor.final_score) for anchor in first.anchors],
            [(anchor.node_id, anchor.final_score) for anchor in second.anchors],
        )

    def test_anchor_reports_include_components_and_selection(self) -> None:
        candidates = score_anchor_candidates(self._weight_candidates(), alpha=0.4, beta=0.3, gamma=0.3)
        select_informative_anchors(candidates, max_anchors=2, min_spacing=0.0)
        temp_dir = ROOT / "tests" / "_tmp" / f"anchors_{uuid.uuid4().hex}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            reports = write_anchor_reports(candidates, temp_dir)
            text = Path(reports["anchor_scores_csv"]).read_text(encoding="utf-8")
            self.assertIn("anchor_id,node_id,haplotype_count,U,D,C,alpha,beta,gamma,final_score,selected", text)
            self.assertTrue(Path(reports["anchor_summary_json"]).exists())
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def test_weights_must_sum_to_one(self) -> None:
        with self.assertRaises(ValueError):
            validate_anchor_weights(0.5, 0.5, 0.5)


if __name__ == "__main__":
    unittest.main()
