from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.anchors_demo import Anchor, run_anchors_demo
from theory.common import THEORY_OUTPUT_DIR, build_demo_context, draw_graph, node_is_cyclic, write_text


def anchor_precedes(anchor_a: Anchor, anchor_b: Anchor, reachability: dict[str, set[str]], cyclic_nodes: set[str]) -> bool:
    if anchor_a.anchor_id == anchor_b.anchor_id or anchor_a.d >= anchor_b.c:
        return False
    if anchor_a.vertex != anchor_b.vertex:
        return anchor_b.vertex in reachability.get(anchor_a.vertex, set())
    return anchor_a.y < anchor_b.x or anchor_a.vertex in cyclic_nodes


def build_precedence_graph(anchors: list[Anchor], reachability: dict[str, set[str]], cyclic_nodes: set[str]) -> nx.DiGraph:
    graph = nx.DiGraph()
    for anchor in anchors:
        graph.add_node(anchor.anchor_id)
    for left in anchors:
        for right in anchors:
            if anchor_precedes(left, right, reachability, cyclic_nodes):
                graph.add_edge(left.anchor_id, right.anchor_id)
    return graph


def run_precedence_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    anchors_result = run_anchors_demo(graph_path, theory_output_dir)
    anchors = anchors_result["anchors"]
    cyclic_nodes = {node for node in context.cyclic_demo_graph.nodes() if node_is_cyclic(context.cyclic_demo_graph, node)}
    precedence_graph = build_precedence_graph(anchors, context.reachability, cyclic_nodes)

    labels = {
        anchor.anchor_id: f"{anchor.anchor_id}\n{anchor.vertex}[{anchor.c}..{anchor.d}]"
        for anchor in anchors
    }
    draw_graph(
        precedence_graph,
        theory_output_dir / "precedence_graph.png",
        "Precedence relation among anchors",
        labels=labels,
    )

    lines = [
        "Precedence Analysis",
        "===================",
        "Simplified rule: anchor A precedes anchor B if query coordinates increase and graph reachability exists.",
        f"Cyclic nodes considered for same-vertex precedence: {sorted(cyclic_nodes)}",
        "",
    ]
    for source, target in precedence_graph.edges():
        lines.append(f"{source} ≺ {target}")
    if precedence_graph.number_of_edges() == 0:
        lines.append("No precedence edges were found in the current anchor set.")
    write_text(theory_output_dir / "precedence_analysis.txt", "\n".join(lines) + "\n")
    return {"context": context, "anchors": anchors, "precedence_graph": precedence_graph, "cyclic_nodes": cyclic_nodes}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute simplified precedence relations among educational anchors.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_precedence_demo(args.graph.resolve() if args.graph else None)
    print(f"Precedence edges: {result['precedence_graph'].number_of_edges()}")


if __name__ == "__main__":
    main()

