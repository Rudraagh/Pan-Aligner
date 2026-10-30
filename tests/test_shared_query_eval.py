import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from shared_query_eval import (
    SequenceGraph,
    WalkStep,
    aggregate,
    minimizer_walk,
    reverse_complement,
    score_proposed_graph,
    score_sequence_graph,
)


class SequenceGraphTests(unittest.TestCase):
    def setUp(self) -> None:
        # Bubble: ACGT -> (G | TT) -> CA
        self.graph = SequenceGraph(
            {"s1": "ACGT", "s2": "G", "s3": "TT", "s4": "CA"},
            [("s1", "+", "s2", "+"), ("s1", "+", "s3", "+"), ("s2", "+", "s4", "+"), ("s3", "+", "s4", "+")],
        )

    def test_substrings_cross_junctions_and_both_alleles(self) -> None:
        self.assertEqual(self.graph.substrings(3), {"ACG", "CGT", "GTG", "GTT", "TGC", "TTT", "TTC", "TCA", "GCA"})
        self.assertTrue(self.graph.contains("ACGTGCA"))
        self.assertTrue(self.graph.contains("acgtttca"))
        self.assertFalse(self.graph.contains("ACGTTGCA"))  # mixes both alleles

    def test_inversion_link_reads_the_reverse_complement(self) -> None:
        graph = SequenceGraph({"a": "AAAC", "b": "GGT"}, [("a", "+", "b", "-")])
        self.assertTrue(graph.contains("AACACC"))  # AAAC + revcomp(GGT)=ACC
        self.assertEqual(reverse_complement("GGTN"), "NACC")

    def test_invalid_graphs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            SequenceGraph({}, [])
        with self.assertRaises(ValueError):
            SequenceGraph({"a": "ACGT"}, [("a", "+", "missing", "+")])
        with self.assertRaises(ValueError):
            self.graph.substrings(0)
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "x.gfa"
            gfa.write_text("S\ta\tACGT\nS\tb\tGG\nL\ta\t+\tb\t+\t2M\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsupported link overlap"):
                SequenceGraph.from_gfa(gfa)
            with self.assertRaises(FileNotFoundError):
                SequenceGraph.from_gfa(Path(tmp) / "missing.gfa")

    def test_rgfa_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            gfa = Path(tmp) / "x.gfa"
            gfa.write_text("S\ta\tacgt\tLN:i:4\tSN:Z:ref\nS\tb\tGG\tLN:i:2\nL\ta\t+\tb\t+\t0M\n", encoding="utf-8")
            graph = SequenceGraph.from_gfa(gfa)
        self.assertEqual(graph.total_bases(), 6)
        self.assertTrue(graph.contains("GTGG"))


class ScoringTests(unittest.TestCase):
    def test_walk_collapses_consecutive_repeats(self) -> None:
        walk = minimizer_walk("ACGTACGTACGT", 3, 4)
        self.assertEqual([step.minimizer for step in walk], ["TAC"])  # TAC at 3 and 7 collapse
        self.assertEqual(minimizer_walk("", 3, 4), [])

    def test_sequence_graph_scores(self) -> None:
        query = "ACGTGCA"
        walk = [WalkStep("ACG", 0), WalkStep("GTG", 2), WalkStep("GCA", 4)]
        graph = SequenceGraph({"s1": "ACGT", "s2": "G", "s3": "TT", "s4": "CA"},
                              [("s1", "+", "s2", "+"), ("s1", "+", "s3", "+"), ("s2", "+", "s4", "+"), ("s3", "+", "s4", "+")])
        scores = score_sequence_graph(walk, query, graph, 3)
        self.assertEqual((scores["node_hits"], scores["transition_hits"]), (3, 2))
        linear = score_sequence_graph(walk, query, SequenceGraph.linear("ACGTTCA"), 3)
        self.assertEqual((linear["node_hits"], linear["transition_hits"]), (1, 0))

    def test_proposed_scores_need_an_edge_not_just_nodes(self) -> None:
        walk = [WalkStep("AAA", 0), WalkStep("CCC", 3), WalkStep("GGG", 6)]
        nodes = {"AAA": "n1", "CCC": "n2", "GGG": "n3"}
        scores = score_proposed_graph(walk, nodes, {("n1", "n2")})
        self.assertEqual((scores["node_recall"], scores["transition_recall"]), (1.0, 0.5))
        empty = score_proposed_graph([], nodes, set())
        self.assertEqual((empty["node_recall"], empty["transition_recall"]), (0.0, 0.0))

    def test_aggregate_pools_hits_and_keeps_per_query_mean(self) -> None:
        rows = [
            {"gene": "G", "graph": "x", "query_minimizers": 10, "node_hits": 10, "node_recall": 1.0, "query_transitions": 9, "transition_hits": 9, "transition_recall": 1.0},
            {"gene": "G", "graph": "x", "query_minimizers": 30, "node_hits": 15, "node_recall": 0.5, "query_transitions": 29, "transition_hits": 0, "transition_recall": 0.0},
        ]
        by_gene = {(row["gene"], row["graph"]): row for row in aggregate(rows)}
        entry = by_gene[("G", "x")]
        self.assertEqual(entry["node_recall_pooled"], 25 / 40)
        self.assertEqual(entry["node_recall_query_mean"], 0.75)
        self.assertEqual(entry["transition_misses"], 29)
        self.assertIn(("ALL", "x"), by_gene)


if __name__ == "__main__":
    unittest.main()
