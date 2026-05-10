from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.anchors_demo import Anchor
from theory.common import THEORY_OUTPUT_DIR, build_demo_context, ensure_dir, write_text
from theory.precedence_demo import run_precedence_demo


def educational_node_length(graph: nx.DiGraph, node: str) -> int:
    # The real graph stores full sequence lengths, which can dwarf the demo.
    # For educational chaining we cap the effective span so gap costs remain interpretable.
    return min(12, max(1, int(graph.nodes[node].get("length", 1))))


def shortest_path_char_distance(graph: nx.DiGraph, source: str, target: str) -> float:
    if source == target:
        return 0.0
    weighted = nx.DiGraph()
    for node in graph.nodes():
        weighted.add_node(node)
    for left, right in graph.edges():
        weighted.add_edge(left, right, weight=1.0)
    try:
        return float(nx.shortest_path_length(weighted, source, target, weight="weight"))
    except nx.NetworkXNoPath:
        return float("inf")


def gap_cost(graph: nx.DiGraph, left: Anchor, right: Anchor) -> tuple[float, float, float]:
    query_gap = max(0, right.c - left.d - 1)
    if left.vertex == right.vertex:
        graph_gap = max(0, right.x - left.y - 1)
    else:
        distance = shortest_path_char_distance(graph, left.vertex, right.vertex)
        left_suffix = max(0, min(2, educational_node_length(graph, left.vertex) - left.y))
        right_prefix = max(0, min(2, right.x - 1))
        graph_gap = distance + left_suffix + right_prefix if distance != float("inf") else float("inf")
    total = query_gap + graph_gap
    return query_gap, graph_gap, total


def co_linear_chain_dp(graph: nx.DiGraph, anchors: list[Anchor], precedence_graph: nx.DiGraph) -> dict:
    ordered = sorted(anchors, key=lambda anchor: (anchor.c, anchor.d, anchor.anchor_id))
    anchor_by_id = {anchor.anchor_id: anchor for anchor in ordered}
    scores = [float(anchor.weight) for anchor in ordered]
    predecessors: list[int | None] = [None] * len(ordered)

    for right_index, right in enumerate(ordered):
        for left_index in range(right_index):
            left = ordered[left_index]
            if not precedence_graph.has_edge(left.anchor_id, right.anchor_id):
                continue
            _, _, total_gap = gap_cost(graph, left, right)
            if total_gap == float("inf"):
                continue
            candidate = scores[left_index] + right.weight - total_gap
            if candidate > scores[right_index]:
                scores[right_index] = candidate
                predecessors[right_index] = left_index

    best_index = max(range(len(ordered)), key=lambda index: scores[index])
    chain: list[Anchor] = []
    current_index: int | None = best_index
    while current_index is not None:
        chain.append(ordered[current_index])
        current_index = predecessors[current_index]
    chain.reverse()
    return {
        "ordered_anchors": ordered,
        "scores": scores,
        "predecessors": predecessors,
        "best_chain": chain,
        "best_score": scores[best_index],
        "anchor_by_id": anchor_by_id,
    }


def iterative_chaining(graph: nx.DiGraph, anchors: list[Anchor], precedence_graph: nx.DiGraph) -> tuple[list[list[float]], list[str]]:
    ordered = sorted(anchors, key=lambda anchor: (anchor.c, anchor.d, anchor.anchor_id))
    scores = [float(anchor.weight) for anchor in ordered]
    history = [scores.copy()]
    logs = ["Iteration 0: " + ", ".join(f"{anchor.anchor_id}={score:.2f}" for anchor, score in zip(ordered, scores))]

    for iteration in range(1, len(ordered) + 1):
        updated_scores = scores.copy()
        for right_index, right in enumerate(ordered):
            for left_index in range(right_index):
                left = ordered[left_index]
                if not precedence_graph.has_edge(left.anchor_id, right.anchor_id):
                    continue
                _, _, total_gap = gap_cost(graph, left, right)
                if total_gap == float("inf"):
                    continue
                updated_scores[right_index] = max(updated_scores[right_index], scores[left_index] + right.weight - total_gap)
        history.append(updated_scores.copy())
        logs.append(
            f"Iteration {iteration}: "
            + ", ".join(f"{anchor.anchor_id}={score:.2f}" for anchor, score in zip(ordered, updated_scores))
        )
        if updated_scores == scores:
            logs.append(f"Converged after {iteration} iterations.")
            break
        scores = updated_scores
    return history, logs


def draw_best_chain_plot(best_chain: list[Anchor], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    plt.figure(figsize=(9, 4))
    x_values = [anchor.c for anchor in best_chain]
    y_values = [anchor.weight for anchor in best_chain]
    labels = [anchor.anchor_id for anchor in best_chain]
    plt.plot(x_values, y_values, marker="o", linewidth=2.0, color="teal")
    for x_value, y_value, label in zip(x_values, y_values, labels):
        plt.text(x_value, y_value + 0.15, label, ha="center", fontsize=8)
    plt.xlabel("Query coordinate")
    plt.ylabel("Anchor weight")
    plt.title("Best co-linear chain (educational DP demo)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def draw_convergence_plot(history: list[list[float]], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    best_per_iteration = [max(scores) for scores in history]
    plt.figure(figsize=(8, 4))
    plt.plot(range(len(best_per_iteration)), best_per_iteration, marker="o", color="purple")
    plt.xlabel("Iteration")
    plt.ylabel("Best chain score")
    plt.title("Iterative chaining convergence")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def run_chaining_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    precedence_result = run_precedence_demo(graph_path, theory_output_dir)
    anchors = precedence_result["anchors"]
    precedence_graph = precedence_result["precedence_graph"]

    dp_result = co_linear_chain_dp(context.cyclic_demo_graph, anchors, precedence_graph)
    history, iteration_logs = iterative_chaining(context.cyclic_demo_graph, anchors, precedence_graph)

    gap_lines = [
        "Gap Cost Analysis",
        "=================",
    ]
    for left, right in zip(dp_result["best_chain"], dp_result["best_chain"][1:]):
        query_gap, graph_gap, total_gap = gap_cost(context.cyclic_demo_graph, left, right)
        gap_lines.append(
            f"{left.anchor_id} -> {right.anchor_id}: query_gap={query_gap:.2f}, graph_gap={graph_gap:.2f}, total={total_gap:.2f}"
        )
    write_text(theory_output_dir / "gap_cost_analysis.txt", "\n".join(gap_lines) + "\n")

    chaining_lines = [
        "Simplified Co-linear Chaining DP",
        "================================",
        "DP recurrence: dp[j] = max score of a chain ending at anchor j.",
        "",
    ]
    for anchor, score in zip(dp_result["ordered_anchors"], dp_result["scores"]):
        chaining_lines.append(f"{anchor.anchor_id} ({anchor.vertex}, [{anchor.c}..{anchor.d}]): score={score:.2f}")
    chaining_lines.append("")
    chaining_lines.append("Best chain: " + " -> ".join(anchor.anchor_id for anchor in dp_result["best_chain"]))
    chaining_lines.append(f"Best chain score: {dp_result['best_score']:.2f}")
    write_text(theory_output_dir / "chaining_results.txt", "\n".join(chaining_lines) + "\n")

    write_text(theory_output_dir / "iteration_log.txt", "\n".join(iteration_logs) + "\n")
    draw_best_chain_plot(dp_result["best_chain"], theory_output_dir / "chaining_graph.png")
    draw_convergence_plot(history, theory_output_dir / "convergence_plot.png")
    return {
        "context": context,
        "best_chain": dp_result["best_chain"],
        "best_score": dp_result["best_score"],
        "history": history,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the simplified co-linear chaining and iterative demos.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_chaining_demo(args.graph.resolve() if args.graph else None)
    print(f"Best chain score: {result['best_score']:.2f}")


if __name__ == "__main__":
    main()
