from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.common import THEORY_OUTPUT_DIR, build_demo_context, draw_graph, write_text


def find_back_edges(graph: nx.DiGraph) -> list[tuple[str, str]]:
    color: dict[str, str] = {node: "white" for node in graph.nodes()}
    back_edges: list[tuple[str, str]] = []

    def dfs(node: str) -> None:
        color[node] = "gray"
        for neighbor in graph.successors(node):
            if color[neighbor] == "white":
                dfs(neighbor)
            elif color[neighbor] == "gray":
                back_edges.append((node, neighbor))
        color[node] = "black"

    for node in graph.nodes():
        if color[node] == "white":
            dfs(node)
    return back_edges


def remove_back_edges(graph: nx.DiGraph, back_edges: list[tuple[str, str]]) -> nx.DiGraph:
    dag = graph.copy()
    dag.remove_edges_from(back_edges)
    return dag


def run_dag_converter_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    back_edges = find_back_edges(context.cyclic_demo_graph)
    dag = remove_back_edges(context.cyclic_demo_graph, back_edges)

    edge_colors = {(u, v): "tomato" for u, v in dag.edges()}
    draw_graph(dag, theory_output_dir / "dag_graph.png", "DAG approximation after DFS back-edge removal", edge_colors=edge_colors)

    lines = [
        "Back-edge Removal Demo",
        "======================",
        f"Back edges identified: {len(back_edges)}",
        "",
    ]
    if back_edges:
        for source, target in back_edges:
            lines.append(f"{source} -> {target}")
    else:
        lines.append("No back edges were detected. The representative graph was already acyclic.")
    lines.append("")
    lines.append(f"DAG verification: {nx.is_directed_acyclic_graph(dag)}")
    write_text(theory_output_dir / "back_edges.txt", "\n".join(lines) + "\n")
    return {"context": context, "back_edges": back_edges, "dag": dag}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run DFS-based back-edge removal on the PanAligner theory graph.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_dag_converter_demo(args.graph.resolve() if args.graph else None)
    print(f"Back edges removed: {len(result['back_edges'])}")


if __name__ == "__main__":
    main()

