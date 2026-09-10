from __future__ import annotations

import shutil
import sys
import unittest
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.graph_consistency import calculate_graph_path_consistency, write_graph_consistency_reports
from proposed.graph_constructor import build_proposed_pangenome_graph
from proposed.models import AssemblyInput, HaplotypePath, ProposedGraph, ProposedGraphEdge, ProposedGraphNode


def _graph(nodes: list[ProposedGraphNode], edges: list[ProposedGraphEdge], paths: list[HaplotypePath]) -> ProposedGraph:
    return ProposedGraph(
        gene="TOY",
        assembly_count=len(paths),
        parameters={},
        backend_used="cpu",
        nodes=nodes,
        edges=edges,
        anchors=[],
        quality={},
        paths=paths,
        metadata={},
    )


def _node(node_id: str, support: int, *, rescued: bool = False) -> ProposedGraphNode:
    return ProposedGraphNode(node_id, node_id, support, 0.0, 0.0, is_rescued=rescued)


def _edge(source: str, target: str, support: int = 0) -> ProposedGraphEdge:
    return ProposedGraphEdge(source, target, support, 0.0)


class GraphPathConsistencyTests(unittest.TestCase):
    def test_identical_paths_have_high_consistency(self) -> None:
        graph = _graph(
            [_node("n1", 3), _node("n2", 3), _node("n3", 3)],
            [_edge("n1", "n2", 3), _edge("n2", "n3", 3)],
            [
                HaplotypePath("h1", "HEALTHY", ["n1", "n2", "n3"], [0, 1, 2]),
                HaplotypePath("h2", "HEALTHY", ["n1", "n2", "n3"], [0, 1, 2]),
                HaplotypePath("h3", "UNHEALTHY", ["n1", "n2", "n3"], [0, 1, 2]),
            ],
        )
        summary = calculate_graph_path_consistency(graph)
        self.assertEqual(summary["mean_node_consistency"], 1.0)
        self.assertEqual(summary["mean_edge_consistency"], 1.0)
        self.assertEqual(summary["mean_path_consistency"], 1.0)

    def test_branching_context_lowers_shared_node_consistency(self) -> None:
        coherent = _graph(
            [_node("p", 2), _node("shared", 2), _node("s", 2)],
            [_edge("p", "shared", 2), _edge("shared", "s", 2)],
            [
                HaplotypePath("h1", "HEALTHY", ["p", "shared", "s"], [0, 1, 2]),
                HaplotypePath("h2", "UNHEALTHY", ["p", "shared", "s"], [0, 1, 2]),
            ],
        )
        branched = _graph(
            [_node("p1", 1), _node("p2", 1), _node("shared", 2), _node("s1", 1), _node("s2", 1)],
            [_edge("p1", "shared"), _edge("p2", "shared"), _edge("shared", "s1"), _edge("shared", "s2")],
            [
                HaplotypePath("h1", "HEALTHY", ["p1", "shared", "s1"], [0, 1, 2]),
                HaplotypePath("h2", "UNHEALTHY", ["p2", "shared", "s2"], [0, 1, 2]),
            ],
        )
        calculate_graph_path_consistency(coherent)
        calculate_graph_path_consistency(branched)
        coherent_shared = next(node for node in coherent.nodes if node.node_id == "shared")
        branched_shared = next(node for node in branched.nodes if node.node_id == "shared")
        self.assertEqual(coherent_shared.consistency_score, 1.0)
        self.assertLess(branched_shared.consistency_score, coherent_shared.consistency_score)
        self.assertEqual(branched_shared.consistency_score, 0.75)

    def test_private_paths_remain_self_consistent(self) -> None:
        graph = _graph(
            [_node("a", 1, rescued=True), _node("b", 1, rescued=True), _node("c", 1, rescued=True), _node("d", 1, rescued=True)],
            [_edge("a", "b"), _edge("c", "d")],
            [
                HaplotypePath("h1", "HEALTHY", ["a", "b"], [0, 1]),
                HaplotypePath("h2", "UNHEALTHY", ["c", "d"], [0, 1]),
            ],
        )
        summary = calculate_graph_path_consistency(graph)
        self.assertTrue(all(node.consistency_score >= 0.875 for node in graph.nodes))
        self.assertTrue(all(path.consistency_score > 0.80 for path in graph.paths))
        self.assertLess(summary["mean_edge_consistency"], 1.0)

    def test_untraversed_cross_haplotype_edge_is_low_consistency(self) -> None:
        graph = _graph(
            [_node("a", 1), _node("b", 1), _node("c", 1), _node("d", 1)],
            [_edge("a", "b"), _edge("c", "d"), _edge("a", "c")],
            [
                HaplotypePath("h1", "HEALTHY", ["a", "b"], [0, 1]),
                HaplotypePath("h2", "UNHEALTHY", ["c", "d"], [0, 1]),
            ],
        )
        calculate_graph_path_consistency(graph)
        cross_edge = next(edge for edge in graph.edges if (edge.source, edge.target) == ("a", "c"))
        self.assertEqual(cross_edge.consistency_score, 0.0)

    def test_rescued_path_and_reports_are_preserved(self) -> None:
        graph = build_proposed_pangenome_graph(
            [
                AssemblyInput("shared", "TOY", "HEALTHY", "ACGTACGTAC"),
                AssemblyInput("private", "TOY", "UNHEALTHY", "TTTGGGCCCA"),
            ],
            k=3,
            window_size=2,
            min_node_support=2,
            backend="cpu",
        )
        private_path = next(path for path in graph.paths if path.assembly_id == "private")
        self.assertGreater(private_path.rescued_node_count, 0)
        self.assertEqual(private_path.valid_transition_count, len(private_path.node_ids) - 1)
        temp_dir = ROOT / "tests" / "_tmp" / f"consistency_{uuid.uuid4().hex}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            reports = write_graph_consistency_reports(graph, temp_dir)
            self.assertEqual(set(reports), {
                "node_consistency_csv", "edge_consistency_csv", "path_consistency_csv", "graph_consistency_json",
            })
            self.assertTrue(all(Path(path).exists() for path in reports.values()))
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def test_repeated_runs_are_deterministic(self) -> None:
        assemblies = [
            AssemblyInput("h1", "TOY", "HEALTHY", "ACGTACGTAC"),
            AssemblyInput("h2", "TOY", "UNHEALTHY", "ACGTTCGTAC"),
        ]
        first = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        second = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        self.assertEqual(first.metadata["graph_consistency"], second.metadata["graph_consistency"])
        self.assertEqual(
            [node.consistency_score for node in first.nodes],
            [node.consistency_score for node in second.nodes],
        )


if __name__ == "__main__":
    unittest.main()
