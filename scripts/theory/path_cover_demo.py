from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.common import THEORY_OUTPUT_DIR, build_demo_context, draw_graph, write_text
from theory.dag_converter_demo import run_dag_converter_demo


def greedy_path_cover(dag: nx.DiGraph) -> list[list[str]]:
    if dag.number_of_nodes() == 0:
        return []

    uncovered = set(dag.nodes())
    order = list(nx.topological_sort(dag))
    rank = {node: index for index, node in enumerate(order)}
    paths: list[list[str]] = []

    while uncovered:
        sources = [node for node in order if node in uncovered and not any(pred in uncovered for pred in dag.predecessors(node))]
        start = sources[0] if sources else min(uncovered, key=lambda node: rank[node])
        current = start
        path = [current]
        uncovered.remove(current)

        while True:
            candidates = sorted(
                [node for node in dag.successors(current) if node in uncovered],
                key=lambda node: rank[node],
            )
            if not candidates:
                break
            current = candidates[0]
            path.append(current)
            uncovered.remove(current)
        paths.append(path)

    return paths


def run_path_cover_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    dag_result = run_dag_converter_demo(graph_path, theory_output_dir)
    dag = dag_result["dag"]
    paths = greedy_path_cover(dag)

    palette = ["#1b9e77", "#d95f02", "#7570b3", "#66a61e", "#e7298a", "#e6ab02"]
    node_colors: dict[str, str] = {}
    for index, path in enumerate(paths):
        for node in path:
            node_colors[node] = palette[index % len(palette)]
    draw_graph(dag, theory_output_dir / "path_cover_graph.png", "Greedy path cover on DAG approximation", node_colors=node_colors)

    covered_nodes = {node for path in paths for node in path}
    lines = [
        "Path Cover Demo",
        "===============",
        f"DAG node count: {dag.number_of_nodes()}",
        f"Path count: {len(paths)}",
        f"All nodes covered: {covered_nodes == set(dag.nodes())}",
        "",
    ]
    for index, path in enumerate(paths, start=1):
        lines.append(f"Path {index}: {' -> '.join(path)}")
    write_text(theory_output_dir / "path_cover.txt", "\n".join(lines) + "\n")
    return {"context": context, "dag": dag, "paths": paths}


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a simplified path cover on the DAG approximation.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_path_cover_demo(args.graph.resolve() if args.graph else None)
    print(f"Path cover size: {len(result['paths'])}")


if __name__ == "__main__":
    main()

