from __future__ import annotations

import argparse
from pathlib import Path

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.common import THEORY_OUTPUT_DIR, build_demo_context, draw_graph, write_text


def tarjan_scc(graph) -> list[list[str]]:
    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    components: list[list[str]] = []

    def strongconnect(node: str) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for neighbor in graph.successors(node):
            if neighbor not in indices:
                strongconnect(neighbor)
                lowlinks[node] = min(lowlinks[node], lowlinks[neighbor])
            elif neighbor in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[neighbor])

        if lowlinks[node] == indices[node]:
            component: list[str] = []
            while stack:
                popped = stack.pop()
                on_stack.remove(popped)
                component.append(popped)
                if popped == node:
                    break
            components.append(sorted(component))

    for node in graph.nodes():
        if node not in indices:
            strongconnect(node)
    return components


def run_scc_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    components = tarjan_scc(context.cyclic_demo_graph)

    palette = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69"]
    node_colors = {}
    for index, component in enumerate(components):
        for node in component:
            node_colors[node] = palette[index % len(palette)]
    draw_graph(
        context.cyclic_demo_graph,
        theory_output_dir / "scc_graph.png",
        "Strongly connected components (educational demo)",
        node_colors=node_colors,
    )

    lines = [
        "Strongly Connected Component Analysis",
        "====================================",
        f"Graph source: {context.source_label}",
        f"SCC count: {len(components)}",
        "",
    ]
    for index, component in enumerate(components, start=1):
        lines.append(f"SCC {index}: {', '.join(component)}")
    write_text(theory_output_dir / "scc_analysis.txt", "\n".join(lines) + "\n")
    return {"context": context, "components": components}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Tarjan SCC analysis on the educational PanAligner demo graph.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_scc_demo(args.graph.resolve() if args.graph else None)
    print(f"SCC count: {len(result['components'])}")


if __name__ == "__main__":
    main()
