from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from theory.anchors_demo import Anchor
from theory.chaining_demo import co_linear_chain_dp, gap_cost, iterative_chaining
from theory.common import compute_reachability, extract_representative_subgraph, node_is_cyclic
from theory.dag_converter_demo import find_back_edges, remove_back_edges
from theory.path_cover_demo import greedy_path_cover
from theory.precedence_demo import anchor_precedes, build_precedence_graph
from theory.scc_demo import tarjan_scc

TRAIN_APP_GRAPH = ROOT / "graphs" / "train" / "app.combined.gfa"


def _graph(edges: list[tuple[str, str]], nodes: list[str] | None = None) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes or [])
    graph.add_edges_from(edges)
    return graph


def _random_graph(seed: int) -> nx.DiGraph:
    return nx.gnp_random_graph(30, 0.08, seed=seed, directed=True)


class TarjanSccTests(unittest.TestCase):
    def test_single_cycle_is_one_component(self) -> None:
        graph = _graph([("a", "b"), ("b", "c"), ("c", "a")])
        self.assertEqual(tarjan_scc(graph), [["a", "b", "c"]])

    def test_dag_has_only_singleton_components(self) -> None:
        graph = _graph([("a", "b"), ("b", "c")])
        self.assertEqual(sorted(tarjan_scc(graph)), [["a"], ["b"], ["c"]])

    def test_two_cycles_joined_by_one_edge_stay_separate(self) -> None:
        graph = _graph([("a", "b"), ("b", "a"), ("c", "d"), ("d", "c"), ("b", "c")])
        self.assertEqual(sorted(tarjan_scc(graph)), [["a", "b"], ["c", "d"]])

    def test_self_loop_is_singleton_but_cyclic(self) -> None:
        graph = _graph([("x", "x"), ("x", "y")])
        self.assertEqual(sorted(tarjan_scc(graph)), [["x"], ["y"]])
        self.assertTrue(node_is_cyclic(graph, "x"))
        self.assertFalse(node_is_cyclic(graph, "y"))

    def test_matches_networkx_on_random_graphs(self) -> None:
        for seed in range(5):
            graph = _random_graph(seed)
            expected = {frozenset(component) for component in nx.strongly_connected_components(graph)}
            self.assertEqual({frozenset(component) for component in tarjan_scc(graph)}, expected)


class BackEdgeRemovalTests(unittest.TestCase):
    def test_cycle_has_exactly_one_back_edge(self) -> None:
        graph = _graph([("a", "b"), ("b", "c"), ("c", "a")])
        self.assertEqual(find_back_edges(graph), [("c", "a")])

    def test_dag_has_no_back_edges(self) -> None:
        graph = _graph([("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])
        self.assertEqual(find_back_edges(graph), [])

    def test_removal_yields_dag_without_mutating_input(self) -> None:
        for seed in range(5):
            graph = _random_graph(seed)
            original_edges = set(graph.edges())
            back_edges = find_back_edges(graph)
            dag = remove_back_edges(graph, back_edges)
            self.assertTrue(nx.is_directed_acyclic_graph(dag))
            self.assertTrue(set(back_edges) <= original_edges)
            self.assertEqual(set(dag.edges()), original_edges - set(back_edges))
            self.assertEqual(set(graph.edges()), original_edges)


class GreedyPathCoverTests(unittest.TestCase):
    def _assert_valid_cover(self, dag: nx.DiGraph, paths: list[list]) -> None:
        covered = [node for path in paths for node in path]
        self.assertEqual(len(covered), len(set(covered)), "paths must be vertex-disjoint")
        self.assertEqual(set(covered), set(dag.nodes()))
        for path in paths:
            for left, right in zip(path, path[1:]):
                self.assertTrue(dag.has_edge(left, right), f"{left}->{right} is not a DAG edge")

    def test_empty_graph_has_empty_cover(self) -> None:
        self.assertEqual(greedy_path_cover(nx.DiGraph()), [])

    def test_chain_is_covered_by_one_path(self) -> None:
        self.assertEqual(greedy_path_cover(_graph([("a", "b"), ("b", "c")])), [["a", "b", "c"]])

    def test_diamond_needs_two_paths(self) -> None:
        dag = _graph([("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])
        paths = greedy_path_cover(dag)
        self._assert_valid_cover(dag, paths)
        self.assertEqual(len(paths), 2)

    def test_cover_is_valid_on_random_dags(self) -> None:
        for seed in range(5):
            graph = _random_graph(seed)
            dag = remove_back_edges(graph, find_back_edges(graph))
            self._assert_valid_cover(dag, greedy_path_cover(dag))


class PrecedenceAndChainingTests(unittest.TestCase):
    """Worked example: u1 -> u2 -> u3 and u1 -> u4, where u4 cannot reach u3."""

    def setUp(self) -> None:
        self.graph = _graph([("u1", "u2"), ("u2", "u3"), ("u1", "u4")])
        for node, length in (("u1", 5), ("u2", 5), ("u3", 6), ("u4", 4)):
            self.graph.nodes[node]["length"] = length
        self.reachability = compute_reachability(self.graph)
        self.a = Anchor("A", "u1", 1, 5, 1, 5, 5, "test")
        self.b = Anchor("B", "u2", 1, 5, 7, 11, 5, "test")
        self.c = Anchor("C", "u4", 1, 4, 8, 11, 4, "test")
        self.d = Anchor("D", "u3", 1, 6, 14, 19, 6, "test")
        self.anchors = [self.a, self.b, self.c, self.d]

    def test_precedence_requires_query_order_and_reachability(self) -> None:
        precedes = lambda left, right: anchor_precedes(left, right, self.reachability, set())
        self.assertTrue(precedes(self.a, self.b))
        self.assertTrue(precedes(self.a, self.c))
        self.assertTrue(precedes(self.b, self.d))
        self.assertFalse(precedes(self.b, self.c), "query intervals overlap")
        self.assertFalse(precedes(self.c, self.d), "u4 cannot reach u3")
        self.assertFalse(precedes(self.d, self.a), "query order is reversed")
        self.assertFalse(precedes(self.a, self.a), "an anchor never precedes itself")

    def test_same_vertex_precedence_needs_graph_order_or_a_cycle(self) -> None:
        forward_left = Anchor("P", "u1", 1, 2, 1, 2, 2, "test")
        forward_right = Anchor("R", "u1", 4, 5, 5, 6, 2, "test")
        backward_right = Anchor("S", "u1", 1, 2, 5, 6, 2, "test")
        backward_left = Anchor("T", "u1", 4, 5, 1, 2, 2, "test")
        self.assertTrue(anchor_precedes(forward_left, forward_right, self.reachability, set()))
        self.assertFalse(anchor_precedes(backward_left, backward_right, self.reachability, set()))
        self.assertTrue(anchor_precedes(backward_left, backward_right, self.reachability, {"u1"}))

    def test_gap_cost_components(self) -> None:
        self.assertEqual(gap_cost(self.graph, self.a, self.b), (1, 1.0, 2.0))
        self.assertEqual(gap_cost(self.graph, self.b, self.d), (2, 1.0, 3.0))
        self.assertEqual(gap_cost(self.graph, self.a, self.d), (8, 2.0, 10.0))
        same_left = Anchor("P", "u1", 1, 2, 1, 2, 2, "test")
        same_right = Anchor("R", "u1", 4, 5, 5, 6, 2, "test")
        self.assertEqual(gap_cost(self.graph, same_left, same_right), (2, 1, 3))
        self.assertEqual(gap_cost(self.graph, self.c, self.d)[1], float("inf"))

    def test_dp_selects_best_co_linear_chain(self) -> None:
        precedence = build_precedence_graph(self.anchors, self.reachability, set())
        result = co_linear_chain_dp(self.graph, self.anchors, precedence)
        self.assertEqual([anchor.anchor_id for anchor in result["best_chain"]], ["A", "B", "D"])
        self.assertEqual(result["best_score"], 11.0)
        self.assertEqual(result["scores"], [5.0, 8.0, 6.0, 11.0])

    def test_iterative_chaining_converges_to_dp_scores(self) -> None:
        precedence = build_precedence_graph(self.anchors, self.reachability, set())
        history, logs = iterative_chaining(self.graph, self.anchors, precedence)
        self.assertEqual(history[-1], co_linear_chain_dp(self.graph, self.anchors, precedence)["scores"])
        self.assertTrue(logs[-1].startswith("Converged after"))


class RepresentativeSubgraphTests(unittest.TestCase):
    def test_small_graph_is_returned_whole(self) -> None:
        graph = _graph([("a", "b"), ("b", "c")])
        self.assertEqual(set(extract_representative_subgraph(graph).edges()), set(graph.edges()))

    def test_result_is_induced_subgraph_with_attributes(self) -> None:
        graph = nx.DiGraph()
        for index in range(60):
            graph.add_node(f"s{index}", sequence="ACGT", length=4)
        graph.add_edges_from((f"s{index}", f"s{index + 1}", {"synthetic": False}) for index in range(59))
        graph.add_edges_from((f"s{index}", f"s{index + 7}", {"synthetic": False}) for index in range(0, 50, 5))

        subgraph = extract_representative_subgraph(graph, limit=18)
        kept = list(subgraph.nodes())
        self.assertEqual(len(kept), 18)
        self.assertEqual(set(subgraph.edges()), set(graph.subgraph(kept).edges()))
        for node in kept:
            self.assertEqual(subgraph.nodes[node], graph.nodes[node])
        for left, right in subgraph.edges():
            self.assertEqual(subgraph.edges[left, right], graph.edges[left, right])


DETERMINISM_SCRIPT = """
import json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from theory.anchors_demo import generate_demo_anchors
from theory.chaining_demo import co_linear_chain_dp
from theory.common import build_demo_context, node_is_cyclic
from theory.precedence_demo import build_precedence_graph

context = build_demo_context(Path(sys.argv[2]))
graph = context.cyclic_demo_graph
anchors = generate_demo_anchors(context)
cyclic = {node for node in graph.nodes() if node_is_cyclic(graph, node)}
chain = co_linear_chain_dp(graph, anchors, build_precedence_graph(anchors, context.reachability, cyclic))
print(json.dumps({
    "nodes": list(context.representative_graph.nodes()),
    "synthetic_back_edges": context.synthetic_back_edges,
    "walk": context.walk_nodes,
    "query": context.query_sequence,
    "best_chain": [anchor.anchor_id for anchor in chain["best_chain"]],
    "best_score": chain["best_score"],
}))
"""


@unittest.skipUnless(TRAIN_APP_GRAPH.exists(), "graphs/train/app.combined.gfa is not available")
class TheoryDeterminismTests(unittest.TestCase):
    def _run_with_hash_seed(self, seed: str) -> dict:
        env = {**os.environ, "PYTHONHASHSEED": seed, "MPLBACKEND": "Agg"}
        completed = subprocess.run(
            [sys.executable, "-c", DETERMINISM_SCRIPT, str(SCRIPTS_DIR), str(TRAIN_APP_GRAPH)],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        return json.loads(completed.stdout)

    def test_demo_context_and_chain_do_not_depend_on_hash_seed(self) -> None:
        baseline = self._run_with_hash_seed("1")
        for seed in ("2", "3"):
            self.assertEqual(self._run_with_hash_seed(seed), baseline, f"output changed with PYTHONHASHSEED={seed}")


if __name__ == "__main__":
    unittest.main()
