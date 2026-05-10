from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.common import (
    THEORY_OUTPUT_DIR,
    DemoGraphContext,
    build_demo_context,
    draw_graph,
    write_context_snapshot,
    write_text,
)


def summarize_graph(graph: nx.DiGraph) -> dict[str, int | bool]:
    components = list(nx.strongly_connected_components(graph))
    return {
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "weakly_connected_components": nx.number_weakly_connected_components(graph),
        "strongly_connected_components": len(components),
        "largest_scc_size": max((len(component) for component in components), default=0),
        "cycle_detected": not nx.is_directed_acyclic_graph(graph),
    }


def build_reachability_section(context: DemoGraphContext) -> str:
    lines = [
        "Reachability analysis",
        "---------------------",
    ]
    for node in sorted(context.cyclic_demo_graph.nodes()):
        reachable = ", ".join(sorted(context.reachability[node]))
        lines.append(f"{node} -> {reachable}")
    return "\n".join(lines)


def run_graph_analysis(graph_path: Path | None = None, output_dir: Path | None = None) -> DemoGraphContext:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    write_context_snapshot(context)

    original_summary = summarize_graph(context.original_graph)
    representative_summary = summarize_graph(context.representative_graph)
    cyclic_summary = summarize_graph(context.cyclic_demo_graph)

    edge_colors = {(u, v): "tomato" for u, v in context.synthetic_back_edges}
    draw_graph(
        context.cyclic_demo_graph,
        theory_output_dir / "cyclic_graph.png",
        "Educational cyclic graph used for PanAligner theory demos",
        edge_colors=edge_colors,
    )

    lines = [
        "Graph Analysis Demo",
        "===================",
        f"Source graph: {context.source_label}",
        f"Source path: {context.source_graph_path}",
        "",
        "Original graph summary",
        f"- nodes: {original_summary['node_count']}",
        f"- edges: {original_summary['edge_count']}",
        f"- SCC count: {original_summary['strongly_connected_components']}",
        f"- largest SCC size: {original_summary['largest_scc_size']}",
        f"- cycle detected: {original_summary['cycle_detected']}",
        "",
        "Representative subgraph summary",
        f"- nodes: {representative_summary['node_count']}",
        f"- edges: {representative_summary['edge_count']}",
        f"- cycle detected: {representative_summary['cycle_detected']}",
        "",
        "Educational cyclic demo graph summary",
        f"- nodes: {cyclic_summary['node_count']}",
        f"- edges: {cyclic_summary['edge_count']}",
        f"- cycle detected: {cyclic_summary['cycle_detected']}",
        f"- synthetic back edges added: {context.synthetic_back_edges or 'none'}",
        "",
        "Interpretation",
        "The real minigraph-derived graphs in this project are mostly DAG-like for the current APP/PSEN1/PSEN2 loci.",
        "To demonstrate the PanAligner paper's cyclic graph algorithms without altering the real graphs, the theory layer builds a small representative subgraph and optionally adds synthetic back edges only inside the educational copy.",
        "",
        build_reachability_section(context),
    ]
    write_text(theory_output_dir / "graph_analysis.txt", "\n".join(lines) + "\n")
    return context


def main() -> None:
    parser = argparse.ArgumentParser(description="Run educational graph analysis for the PanAligner theory module.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA to analyze.")
    args = parser.parse_args()
    context = run_graph_analysis(args.graph.resolve() if args.graph else None)
    print(f"Theory graph prepared from {context.source_graph_path}")
    print(f"Cyclic demo graph written to {(THEORY_OUTPUT_DIR / 'cyclic_graph.png').resolve()}")


if __name__ == "__main__":
    main()

