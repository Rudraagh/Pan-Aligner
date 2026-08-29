from __future__ import annotations

import sys
import json
import shutil
import unittest
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.graph_constructor import build_proposed_pangenome_graph, validate_proposed_graph, write_proposed_gfa, write_proposed_graph_bundle
from proposed.models import AssemblyInput


class ProposedGraphConstructorTests(unittest.TestCase):
    def _assemblies(self) -> list[AssemblyInput]:
        return [
            AssemblyInput("hap1", "TOY", "HEALTHY", "ACGTACGTACGT"),
            AssemblyInput("hap2", "TOY", "HEALTHY", "ACGTACGTTCGT"),
            AssemblyInput("hap3", "TOY", "UNHEALTHY", "ACGTTCGTACGT"),
            AssemblyInput("hap4", "TOY", "UNHEALTHY", "ACGTTCGTTCGT"),
        ]

    def test_graph_builds_nodes_edges_and_anchors(self) -> None:
        graph = build_proposed_pangenome_graph(
            self._assemblies(),
            k=3,
            window_size=3,
            min_node_support=2,
            max_anchors=4,
            min_anchor_spacing=1.5,
            backend="cpu",
        )
        self.assertGreater(len(graph.nodes), 0)
        self.assertGreater(len(graph.edges), 0)
        self.assertGreater(len(graph.anchors), 0)
        self.assertEqual(graph.gene, "TOY")
        self.assertIn("overall_quality_score", graph.quality)
        self.assertEqual(len(graph.paths), len(self._assemblies()))
        self.assertFalse(validate_proposed_graph(graph))

    def test_private_minimizers_are_rescued_into_their_haplotype_path(self) -> None:
        assemblies = [
            AssemblyInput("shared", "TOY", "HEALTHY", "ACGTACGTAC"),
            AssemblyInput("divergent", "TOY", "UNHEALTHY", "TTTGGGCCCA"),
        ]
        graph = build_proposed_pangenome_graph(
            assemblies,
            k=3,
            window_size=2,
            min_node_support=2,
            backend="cpu",
        )
        paths = {path.assembly_id: path for path in graph.paths}
        private_node_ids = {node.node_id for node in graph.nodes if node.support == 1}
        self.assertTrue(private_node_ids)
        self.assertTrue(set(paths["divergent"].node_ids).intersection(private_node_ids))
        self.assertGreater(graph.metadata["rescued_private_node_count"], 0)
        self.assertEqual(graph.metadata["discarded"]["minimizer_records"], 0)

    def test_identical_haplotypes_share_the_same_path_nodes(self) -> None:
        assemblies = [
            AssemblyInput("hap_a", "TOY", "HEALTHY", "ACGTACGTAC"),
            AssemblyInput("hap_b", "TOY", "HEALTHY", "ACGTACGTAC"),
        ]
        graph = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        self.assertEqual(graph.paths[0].node_ids, graph.paths[1].node_ids)
        self.assertEqual(graph.metadata["haplotype_specific_node_count"], 0)
        self.assertEqual(graph.metadata["rescued_node_count"], 0)

    def test_divergent_haplotypes_create_distinct_path_structure(self) -> None:
        assemblies = [
            AssemblyInput("hap_a", "TOY", "HEALTHY", "ACGTACGTAC"),
            AssemblyInput("hap_b", "TOY", "UNHEALTHY", "TTTGGGCCCA"),
        ]
        graph = build_proposed_pangenome_graph(assemblies, k=3, window_size=2, min_node_support=2, backend="cpu")
        self.assertNotEqual(graph.paths[0].node_ids, graph.paths[1].node_ids)
        self.assertFalse(validate_proposed_graph(graph))

    def test_graph_can_be_exported_to_gfa(self) -> None:
        graph = build_proposed_pangenome_graph(
            self._assemblies(),
            k=3,
            window_size=3,
            min_node_support=2,
            max_anchors=4,
            min_anchor_spacing=1.5,
            backend="cpu",
        )
        temp_dir = ROOT / "tests" / "_tmp" / f"graph_export_{uuid.uuid4().hex}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        try:
            bundle = write_proposed_graph_bundle(graph, temp_dir)
            gfa_path = write_proposed_gfa(graph, temp_dir / "toy.gfa")
            text = gfa_path.read_text(encoding="utf-8")
            self.assertIn("\nS\t", text)
            self.assertIn("\nL\t", text)
            p_records = [line.split("\t") for line in text.splitlines() if line.startswith("P\t")]
            self.assertEqual(len(p_records), len(graph.paths))
            known_node_ids = {node.node_id for node in graph.nodes}
            for record, path in zip(p_records, graph.paths):
                self.assertEqual(record[1], path.assembly_id)
                self.assertEqual([segment[:-1] for segment in record[2].split(",")], path.node_ids)
                self.assertTrue(all(segment[:-1] in known_node_ids for segment in record[2].split(",")))
            payload = json.loads(Path(bundle["graph_json"]).read_text(encoding="utf-8"))
            self.assertEqual(payload["metadata"]["path_count"], len(graph.paths))
            self.assertEqual(len(payload["paths"]), len(graph.paths))
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
