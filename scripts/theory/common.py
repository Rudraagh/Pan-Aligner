from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, read_json, write_json


THEORY_OUTPUT_DIR = OUTPUTS_DIR / "theory"
REPORTS_OUTPUT_DIR = OUTPUTS_DIR / "reports"


@dataclass
class DemoGraphContext:
    source_graph_path: Path
    source_label: str
    original_graph: nx.DiGraph
    representative_graph: nx.DiGraph
    cyclic_demo_graph: nx.DiGraph
    synthetic_back_edges: list[tuple[str, str]]
    query_sequence: str
    walk_nodes: list[str]
    reachability: dict[str, set[str]]


def parse_gfa_with_sequences(gfa_path: Path) -> nx.DiGraph:
    graph = nx.DiGraph()
    with gfa_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if fields[0] == "S":
                node_id = fields[1]
                sequence = fields[2]
                graph.add_node(node_id, sequence=sequence, length=len(sequence))
            elif fields[0] == "L":
                from_node = fields[1]
                to_node = fields[3]
                graph.add_edge(from_node, to_node, synthetic=False)
    return graph


def choose_source_graph(graph_path: Path | None = None) -> tuple[Path, str]:
    if graph_path is not None:
        resolved = graph_path.resolve()
        return resolved, resolved.stem

    for manifest_name in ("train_graph_manifest.json", "graph_manifest.json"):
        manifest_path = METADATA_DIR / manifest_name
        if manifest_path.exists():
            manifest = read_json(manifest_path)
            if "APP" in manifest["graphs"]:
                selected = Path(manifest["graphs"]["APP"]["combined_graph"]).resolve()
                if selected.exists():
                    return selected, "APP combined"
            first_gene = sorted(manifest["graphs"])[0]
            selected = Path(manifest["graphs"][first_gene]["combined_graph"]).resolve()
            if selected.exists():
                return selected, f"{first_gene} combined"

    fallback = ROOT / "graphs" / "app.combined.gfa"
    if fallback.exists():
        return fallback.resolve(), "APP combined"

    available = sorted((ROOT / "graphs").glob("*.combined.gfa"))
    if not available:
        raise FileNotFoundError("No combined graph found for theory demos. Build graphs first or pass --theory-graph.")
    return available[0].resolve(), available[0].stem


def _topological_or_node_order(graph: nx.DiGraph) -> list[str]:
    if nx.is_directed_acyclic_graph(graph):
        return list(nx.topological_sort(graph))
    return list(graph.nodes())


def extract_representative_subgraph(graph: nx.DiGraph, limit: int = 18) -> nx.DiGraph:
    if graph.number_of_nodes() <= limit:
        return graph.copy()

    ordered_nodes = sorted(
        graph.nodes(),
        key=lambda node: (
            -(graph.in_degree(node) + graph.out_degree(node)),
            node,
        ),
    )
    selected: list[str] = []
    for node in ordered_nodes:
        if node not in selected:
            selected.append(node)
        if len(selected) >= min(6, limit):
            break

    # Include shortest-path connectors so the educational graph remains readable.
    for left, right in zip(selected, selected[1:]):
        try:
            connector = nx.shortest_path(graph, left, right)
        except nx.NetworkXNoPath:
            continue
        for node in connector:
            if node not in selected:
                selected.append(node)
            if len(selected) >= limit:
                break
        if len(selected) >= limit:
            break

    for node in _topological_or_node_order(graph):
        if node not in selected:
            selected.append(node)
        if len(selected) >= limit:
            break

    return graph.subgraph(selected[:limit]).copy()


def augment_with_synthetic_cycles(graph: nx.DiGraph) -> tuple[nx.DiGraph, list[tuple[str, str]]]:
    augmented = graph.copy()
    synthetic_edges: list[tuple[str, str]] = []
    if not nx.is_directed_acyclic_graph(augmented) or augmented.number_of_nodes() < 4:
        return augmented, synthetic_edges

    order = list(nx.topological_sort(augmented))
    candidate_pairs = [
        (order[-1], order[1]),
        (order[-2], order[0]),
    ]
    for source, target in candidate_pairs:
        if source != target and not augmented.has_edge(source, target):
            augmented.add_edge(source, target, synthetic=True)
            synthetic_edges.append((source, target))
    return augmented, synthetic_edges


def build_demo_walk(graph: nx.DiGraph, max_steps: int = 8) -> list[str]:
    if graph.number_of_nodes() == 0:
        return []

    try:
        cycle_edges = nx.find_cycle(graph, orientation="original")
    except nx.NetworkXNoCycle:
        cycle_edges = []

    if cycle_edges:
        walk = [cycle_edges[0][0]]
        for source, target, *_ in cycle_edges:
            if walk[-1] != source:
                walk.append(source)
            walk.append(target)
        return walk[: max_steps + 1]

    order = _topological_or_node_order(graph)
    return order[:max_steps]


def build_demo_query(graph: nx.DiGraph, walk_nodes: list[str], segment_length: int = 10) -> str:
    segments: list[str] = []
    for index, node in enumerate(walk_nodes):
        sequence = graph.nodes[node].get("sequence", "A")
        segments.append(sequence[: min(segment_length, len(sequence))] or "A")
        if index < len(walk_nodes) - 1 and index % 2 == 0:
            segments.append("N")
    return "".join(segments)


def compute_reachability(graph: nx.DiGraph) -> dict[str, set[str]]:
    reachability: dict[str, set[str]] = {}
    for node in graph.nodes():
        reachability[node] = nx.descendants(graph, node) | {node}
    return reachability


def build_demo_context(graph_path: Path | None = None) -> DemoGraphContext:
    source_graph_path, source_label = choose_source_graph(graph_path)
    original_graph = parse_gfa_with_sequences(source_graph_path)
    representative_graph = extract_representative_subgraph(original_graph)
    cyclic_demo_graph, synthetic_back_edges = augment_with_synthetic_cycles(representative_graph)
    walk_nodes = build_demo_walk(cyclic_demo_graph)
    query_sequence = build_demo_query(cyclic_demo_graph, walk_nodes)
    reachability = compute_reachability(cyclic_demo_graph)
    return DemoGraphContext(
        source_graph_path=source_graph_path,
        source_label=source_label,
        original_graph=original_graph,
        representative_graph=representative_graph,
        cyclic_demo_graph=cyclic_demo_graph,
        synthetic_back_edges=synthetic_back_edges,
        query_sequence=query_sequence,
        walk_nodes=walk_nodes,
        reachability=reachability,
    )


def draw_graph(
    graph: nx.DiGraph,
    output_path: Path,
    title: str,
    node_colors: dict[str, str] | None = None,
    edge_colors: dict[tuple[str, str], str] | None = None,
    labels: dict[str, str] | None = None,
) -> None:
    ensure_dir(output_path.parent)
    plt.figure(figsize=(11, 8))
    pos = nx.spring_layout(graph, seed=42)
    node_colors = node_colors or {}
    edge_colors = edge_colors or {}
    labels = labels or {node: node for node in graph.nodes()}

    nx.draw_networkx_nodes(
        graph,
        pos,
        node_color=[node_colors.get(node, "skyblue") for node in graph.nodes()],
        node_size=800,
        edgecolors="black",
        linewidths=0.8,
    )
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=8)
    nx.draw_networkx_edges(
        graph,
        pos,
        edge_color=[edge_colors.get((u, v), "gray") for u, v in graph.edges()],
        arrows=True,
        width=1.8,
        arrowsize=16,
        connectionstyle="arc3,rad=0.08",
    )
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def node_is_cyclic(graph: nx.DiGraph, node: str) -> bool:
    for component in nx.strongly_connected_components(graph):
        if node in component:
            return len(component) > 1 or graph.has_edge(node, node)
    return False


def write_context_snapshot(context: DemoGraphContext) -> None:
    snapshot = {
        "source_graph_path": str(context.source_graph_path),
        "source_label": context.source_label,
        "original_graph_nodes": context.original_graph.number_of_nodes(),
        "original_graph_edges": context.original_graph.number_of_edges(),
        "representative_graph_nodes": context.representative_graph.number_of_nodes(),
        "representative_graph_edges": context.representative_graph.number_of_edges(),
        "cyclic_demo_nodes": context.cyclic_demo_graph.number_of_nodes(),
        "cyclic_demo_edges": context.cyclic_demo_graph.number_of_edges(),
        "synthetic_back_edges": context.synthetic_back_edges,
        "walk_nodes": context.walk_nodes,
        "query_sequence": context.query_sequence,
    }
    write_json(THEORY_OUTPUT_DIR / "theory_context.json", snapshot)
