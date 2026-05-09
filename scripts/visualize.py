from __future__ import annotations

import argparse
from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, read_json, write_json
from parse_gaf import parse_gaf_file


def parse_gfa(gfa_path: Path) -> nx.DiGraph:
    graph = nx.DiGraph()
    with gfa_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            record_type = fields[0]
            if record_type == "S":
                node_id = fields[1]
                sequence = fields[2]
                graph.add_node(node_id, length=len(sequence))
            elif record_type == "L":
                from_node = fields[1]
                from_orient = fields[2]
                to_node = fields[3]
                to_orient = fields[4]
                graph.add_edge(from_node, to_node, from_orient=from_orient, to_orient=to_orient)
    return graph


def graph_statistics(graph: nx.DiGraph) -> dict:
    weak_components = list(nx.weakly_connected_components(graph))
    sccs = list(nx.strongly_connected_components(graph))
    cycle_detected = not nx.is_directed_acyclic_graph(graph)
    largest_component = max((len(component) for component in weak_components), default=0)
    largest_scc = max((len(component) for component in sccs), default=0)
    return {
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "weakly_connected_components": len(weak_components),
        "strongly_connected_components": len(sccs),
        "largest_component_size": largest_component,
        "largest_scc_size": largest_scc,
        "cycle_detected": cycle_detected,
        "self_loops": nx.number_of_selfloops(graph),
    }


def reduced_plot_graph(graph: nx.DiGraph, max_nodes: int = 400) -> nx.DiGraph:
    if graph.number_of_nodes() <= max_nodes:
        return graph.copy()
    ranked_nodes = sorted(graph.degree, key=lambda item: item[1], reverse=True)[:max_nodes]
    return graph.subgraph([node for node, _ in ranked_nodes]).copy()


def draw_graph(graph: nx.DiGraph, output_path: Path, highlighted_nodes: set[str] | None = None, title: str = "") -> None:
    ensure_dir(output_path.parent)
    plot_graph = reduced_plot_graph(graph)
    if plot_graph.number_of_nodes() == 0:
        return

    highlighted_nodes = highlighted_nodes or set()
    size_factor = max(1.0, sqrt(plot_graph.number_of_nodes()))
    pos = nx.spring_layout(plot_graph, seed=42, k=1 / size_factor, iterations=80)
    node_colors = ["tomato" if node in highlighted_nodes else "skyblue" for node in plot_graph.nodes()]
    node_sizes = [180 if node in highlighted_nodes else 120 for node in plot_graph.nodes()]

    plt.figure(figsize=(12, 9))
    nx.draw_networkx(
        plot_graph,
        pos=pos,
        with_labels=plot_graph.number_of_nodes() <= 80,
        node_color=node_colors,
        node_size=node_sizes,
        font_size=7,
        arrows=True,
        width=0.8,
    )
    plt.title(title or output_path.stem)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def visualize_graph(gfa_path: Path, output_prefix: Path) -> dict:
    graph = parse_gfa(gfa_path)
    stats = graph_statistics(graph)
    write_json(output_prefix.parent / f"{output_prefix.name}.stats.json", stats)
    draw_graph(graph, output_prefix.parent / f"{output_prefix.name}.png", title=output_prefix.name)
    return stats


def visualize_alignment(graph_gfa: Path, gaf_path: Path, output_png: Path) -> None:
    graph = parse_gfa(graph_gfa)
    results = parse_gaf_file(gaf_path)
    highlighted_nodes = set()
    for result in results[:1]:
        highlighted_nodes.update(result.traversed_nodes)
    draw_graph(graph, output_png, highlighted_nodes=highlighted_nodes, title=output_png.stem)


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize pangenome graphs and PanAligner alignments.")
    parser.add_argument(
        "--graph-manifest",
        type=Path,
        default=METADATA_DIR / "graph_manifest.json",
        help="Graph manifest produced by build_graph.py.",
    )
    parser.add_argument(
        "--gaf",
        type=Path,
        default=None,
        help="Optional GAF alignment file to overlay on the selected graph.",
    )
    parser.add_argument(
        "--gene",
        type=str,
        default=None,
        help="Optional gene name for alignment overlay. Required when --gaf is set and multiple graphs exist.",
    )
    args = parser.parse_args()

    graph_manifest = read_json(args.graph_manifest.resolve())
    graph_output_dir = ensure_dir(OUTPUTS_DIR / "graphs")

    for gene, entry in graph_manifest["graphs"].items():
        for label, graph_key in (("healthy", "healthy_graph"), ("unhealthy", "unhealthy_graph"), ("combined", "combined_graph")):
            output_prefix = graph_output_dir / f"{gene.lower()}.{label}"
            stats = visualize_graph(Path(entry[graph_key]), output_prefix)
            print(f"{gene} {label}: nodes={stats['node_count']} edges={stats['edge_count']} cycles={stats['cycle_detected']}")

    if args.gaf:
        if args.gene is None:
            raise ValueError("Pass --gene when using --gaf so the script knows which graph to overlay.")
        selected = graph_manifest["graphs"][args.gene.upper()]["combined_graph"]
        overlay_path = ensure_dir(OUTPUTS_DIR / "alignments") / f"{args.gene.lower()}.alignment.png"
        visualize_alignment(Path(selected), args.gaf.resolve(), overlay_path)
        print(f"Alignment visualization written to {overlay_path.resolve()}")


if __name__ == "__main__":
    main()
