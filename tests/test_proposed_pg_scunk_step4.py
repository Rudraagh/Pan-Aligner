from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.models import AssemblyInput
from proposed.quality_feedback import (
    canonical_configuration,
    diagnose_quality,
    generate_feedback_candidates,
    run_pg_scunk_feedback_loop,
)


def _graph(quality: float, marker: str) -> SimpleNamespace:
    return SimpleNamespace(
        marker=marker,
        paths=[object(), object()],
        quality={
            "overall_quality_score": quality,
            "anchor_coverage_ratio": 0.4,
            "mean_node_consistency": 0.7,
            "mean_edge_consistency": 0.8,
            "mean_path_consistency": 0.8,
            "mean_node_support_ratio": 0.7,
            "branching_node_ratio": 0.1,
            "node_count": 4,
            "edge_count": 3,
            "anchor_count": 2,
        },
    )


class PGScunkStep4Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.assemblies = [AssemblyInput("h1", "TOY", "HEALTHY", "ACGTACGTAC"), AssemblyInput("h2", "TOY", "UNHEALTHY", "ACGTTCGTAC")]
        self.parameters = {"k": 5, "window_size": 3, "min_node_support": 1, "max_anchors": 4, "min_anchor_spacing": 2.0, "anchor_alpha": 0.35, "anchor_beta": 0.35, "anchor_gamma": 0.30}

    def _run(self, qualities: list[float], **kwargs):
        components = kwargs.pop("components", {})
        graphs = [_graph(value, f"g{index}") for index, value in enumerate(qualities)]
        for graph in graphs:
            graph.quality.update(components)
        quality_threshold = kwargs.pop("quality_threshold", 0.99)
        base_parameters = kwargs.pop("base_parameters", self.parameters)
        with patch("proposed.quality_feedback.build_proposed_pangenome_graph", side_effect=graphs) as builder:
            result = run_pg_scunk_feedback_loop(self.assemblies, base_parameters=base_parameters, backend="cpu", quality_threshold=quality_threshold, min_anchor_coverage=kwargs.pop("min_anchor_coverage", 0.3), max_rounds=kwargs.pop("max_rounds", 2), **kwargs)
        return result, builder

    def test_initial_threshold_is_immediate_success(self) -> None:
        result, builder = self._run([0.95], quality_threshold=0.9)
        self.assertTrue(result["threshold_reached"])
        self.assertEqual(result["unique_configurations_evaluated"], 1)
        self.assertEqual(builder.call_count, 1)

    def test_improvement_updates_best_so_far(self) -> None:
        result, _ = self._run([0.30, 0.75, 0.40, 0.35])
        self.assertEqual(result["best_overall_quality"], 0.75)
        self.assertEqual(result["best_graph"].marker, "g1")
        self.assertEqual(result["best_evaluation"]["candidate_id"], "c002")

    def test_multiple_rounds_evaluate_at_least_three_unique_configurations(self) -> None:
        result, _ = self._run([0.30, 0.31, 0.32, 0.33], max_rounds=2)
        self.assertGreaterEqual(result["unique_configurations_evaluated"], 3)
        keys = [item["configuration_key"] for item in result["evaluations"]]
        self.assertEqual(len(keys), len(set(keys)))

    def test_duplicate_candidate_is_filtered(self) -> None:
        fixed = {**self.parameters, "k": 3, "window_size": 2}
        visited = {canonical_configuration(fixed)}
        diagnosis = {"action": "retune_path_resolution", "reason": "test"}
        candidates = generate_feedback_candidates(fixed, diagnosis, visited)
        self.assertEqual(candidates, [])

    def test_worse_later_candidate_cannot_replace_best(self) -> None:
        result, _ = self._run([0.80, 0.20, 0.30, 0.10])
        self.assertEqual(result["best_overall_quality"], 0.80)
        self.assertEqual(result["best_graph"].marker, "g0")
        self.assertEqual(result["final_quality"], 0.10)
        self.assertFalse(result["threshold_reached"])

    def test_stagnation_stops_after_patience(self) -> None:
        result, _ = self._run([0.40, 0.40, 0.40, 0.40], stagnation_patience=1)
        self.assertEqual(result["stagnation_reason"], "stagnation_patience_exceeded")

    def test_exhaustion_stops_when_no_unseen_configuration_exists(self) -> None:
        fixed = {**self.parameters, "k": 3, "window_size": 2, "min_node_support": 1}
        result, _ = self._run([0.30], base_parameters=fixed, max_rounds=4, components={"anchor_coverage_ratio": 0.9, "mean_node_consistency": 0.9, "mean_node_support_ratio": 0.9, "mean_edge_consistency": 0.9, "mean_path_consistency": 0.1, "branching_node_ratio": 0.0})
        self.assertEqual(result["stagnation_reason"], "no_unseen_configurations")

    def test_diagnosis_actions_are_distinct_and_deterministic(self) -> None:
        anchor = diagnose_quality({"anchor_coverage_ratio": 0.1, "mean_node_consistency": 0.9, "mean_node_support_ratio": 0.9, "mean_edge_consistency": 0.9, "mean_path_consistency": 0.9, "branching_node_ratio": 0.0})
        node = diagnose_quality({"anchor_coverage_ratio": 0.9, "mean_node_consistency": 0.1, "mean_node_support_ratio": 0.9, "mean_edge_consistency": 0.9, "mean_path_consistency": 0.9, "branching_node_ratio": 0.0})
        self.assertEqual(anchor["action"], "retune_anchor_selection")
        self.assertEqual(node["action"], "retune_minimizer_nodes")

    def test_same_mocked_input_has_deterministic_trajectory(self) -> None:
        first, _ = self._run([0.30, 0.40, 0.35, 0.36])
        second, _ = self._run([0.30, 0.40, 0.35, 0.36])
        self.assertEqual(first["evaluations"], second["evaluations"])
        self.assertEqual(first["best_configuration"], second["best_configuration"])

    def test_every_controller_parameter_reaches_graph_builder(self) -> None:
        _, builder = self._run([0.95], quality_threshold=0.9)
        kwargs = builder.call_args.kwargs
        self.assertEqual({"k", "window_size", "min_node_support", "max_anchors", "min_anchor_spacing", "anchor_alpha", "anchor_beta", "anchor_gamma"}, set(kwargs).intersection({"k", "window_size", "min_node_support", "max_anchors", "min_anchor_spacing", "anchor_alpha", "anchor_beta", "anchor_gamma"}))


if __name__ == "__main__":
    unittest.main()
