from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from .models import ProposedGraph


LOW_CONSISTENCY_THRESHOLD = 0.50


def _mean(values: list[float], default: float = 0.0) -> float:
    return sum(values) / len(values) if values else default


def _dominant_fraction(values: list[str]) -> float:
    if not values:
        return 1.0
    return max(Counter(values).values()) / len(values)


def _node_status(node) -> str:
    if node.is_rescued and node.support == 1:
        return "rescued_private"
    if node.is_rescued:
        return "rescued_low_support"
    if node.support == 1:
        return "haplotype_specific"
    return "shared"


def calculate_graph_path_consistency(graph: ProposedGraph) -> dict:
    """Annotate a graph using only its observed HaplotypePath traversals.

    Node score: 0.25 * path coverage + 0.25 * predecessor agreement +
    0.25 * successor agreement + 0.25 * local transition validity.
    Edge score: 0.50 * path support fraction + 0.50 * source-transition
    agreement. Path score: 0.40 * valid-transition fraction + 0.30 * mean
    node score + 0.30 * mean edge score.
    """

    path_count = len(graph.paths)
    node_by_id = {node.node_id: node for node in graph.nodes}
    edge_by_pair = {(edge.source, edge.target): edge for edge in graph.edges}
    node_path_ids: dict[str, set[str]] = defaultdict(set)
    predecessor_contexts: dict[str, list[str]] = defaultdict(list)
    successor_contexts: dict[str, list[str]] = defaultdict(list)
    local_validity: dict[str, list[float]] = defaultdict(list)
    observed_edge_paths: dict[tuple[str, str], set[str]] = defaultdict(set)

    for path in graph.paths:
        for index, node_id in enumerate(path.node_ids):
            if node_id not in node_by_id:
                continue
            node_path_ids[node_id].add(path.assembly_id)
            predecessor_contexts[node_id].append(path.node_ids[index - 1] if index else "START")
            successor_contexts[node_id].append(
                path.node_ids[index + 1] if index + 1 < len(path.node_ids) else "END"
            )
            expected_transitions = 0
            valid_transitions = 0
            if index:
                expected_transitions += 1
                valid_transitions += int((path.node_ids[index - 1], node_id) in edge_by_pair)
            if index + 1 < len(path.node_ids):
                expected_transitions += 1
                valid_transitions += int((node_id, path.node_ids[index + 1]) in edge_by_pair)
            local_validity[node_id].append(
                1.0 if expected_transitions == 0 else valid_transitions / expected_transitions
            )
        for source, target in zip(path.node_ids, path.node_ids[1:]):
            observed_edge_paths[(source, target)].add(path.assembly_id)

    for node in graph.nodes:
        node.path_coverage = len(node_path_ids[node.node_id]) / max(1, path_count)
        node.predecessor_consistency = _dominant_fraction(predecessor_contexts[node.node_id])
        node.successor_consistency = _dominant_fraction(successor_contexts[node.node_id])
        node.local_path_consistency = _mean(local_validity[node.node_id], default=0.0)
        node.consistency_score = (
            0.25 * node.path_coverage
            + 0.25 * node.predecessor_consistency
            + 0.25 * node.successor_consistency
            + 0.25 * node.local_path_consistency
        )

    source_path_counts = {node_id: len(path_ids) for node_id, path_ids in node_path_ids.items()}
    for edge in graph.edges:
        observed_paths = observed_edge_paths[(edge.source, edge.target)]
        edge.haplotypes = sorted(observed_paths)
        edge.support = len(observed_paths)
        edge.path_coverage = edge.support / max(1, path_count)
        edge.transition_consistency = edge.support / max(1, source_path_counts.get(edge.source, 0))
        edge.consistency_score = 0.50 * edge.path_coverage + 0.50 * edge.transition_consistency

    for path in graph.paths:
        transitions = list(zip(path.node_ids, path.node_ids[1:]))
        valid_edges = [edge_by_pair[pair] for pair in transitions if pair in edge_by_pair]
        path.valid_transition_count = len(valid_edges)
        path.transition_fraction = len(valid_edges) / len(transitions) if transitions else 1.0
        path.mean_node_consistency = _mean(
            [node_by_id[node_id].consistency_score for node_id in path.node_ids if node_id in node_by_id],
            default=0.0,
        )
        path.mean_edge_consistency = _mean(
            [edge.consistency_score for edge in valid_edges],
            default=1.0,
        )
        path.rescued_node_count = sum(
            node_by_id[node_id].is_rescued for node_id in path.node_ids if node_id in node_by_id
        )
        path.questionable_node_count = sum(
            node_by_id[node_id].consistency_score < LOW_CONSISTENCY_THRESHOLD
            for node_id in path.node_ids
            if node_id in node_by_id
        )
        path.consistency_score = (
            0.40 * path.transition_fraction
            + 0.30 * path.mean_node_consistency
            + 0.30 * path.mean_edge_consistency
        )

    summary = {
        "formula": {
            "node": "0.25*path_coverage + 0.25*predecessor_consistency + 0.25*successor_consistency + 0.25*local_path_consistency",
            "edge": "0.50*path_support_fraction + 0.50*source_transition_agreement",
            "path": "0.40*valid_transition_fraction + 0.30*mean_node_consistency + 0.30*mean_edge_consistency",
            "graph": "0.40*mean_node_consistency + 0.30*mean_edge_consistency + 0.30*mean_path_consistency",
        },
        "mean_node_consistency": _mean([node.consistency_score for node in graph.nodes]),
        "mean_edge_consistency": _mean([edge.consistency_score for edge in graph.edges], default=1.0),
        "mean_path_consistency": _mean([path.consistency_score for path in graph.paths]),
        "low_consistency_node_fraction": (
            sum(node.consistency_score < LOW_CONSISTENCY_THRESHOLD for node in graph.nodes) / len(graph.nodes)
            if graph.nodes
            else 0.0
        ),
        "low_consistency_edge_fraction": (
            sum(edge.consistency_score < LOW_CONSISTENCY_THRESHOLD for edge in graph.edges) / len(graph.edges)
            if graph.edges
            else 0.0
        ),
        "node_count": len(graph.nodes),
        "edge_count": len(graph.edges),
        "path_count": path_count,
    }
    summary["graph_path_consistency"] = (
        0.40 * summary["mean_node_consistency"]
        + 0.30 * summary["mean_edge_consistency"]
        + 0.30 * summary["mean_path_consistency"]
    )
    graph.metadata["graph_consistency"] = summary
    return summary


def write_graph_consistency_reports(graph: ProposedGraph, output_dir: Path) -> dict:
    """Write explainable path-derived consistency reports for one graph bundle."""

    output_dir.mkdir(parents=True, exist_ok=True)
    summary = calculate_graph_path_consistency(graph)
    node_path = output_dir / "node_consistency.csv"
    edge_path = output_dir / "edge_consistency.csv"
    path_path = output_dir / "path_consistency.csv"
    graph_path = output_dir / "graph_consistency.json"

    with node_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "node_id", "haplotype_count", "path_fraction", "predecessor_consistency",
            "successor_consistency", "local_consistency", "consistency_score", "status",
        ])
        writer.writeheader()
        for node in graph.nodes:
            writer.writerow({
                "node_id": node.node_id,
                "haplotype_count": len(node.haplotypes),
                "path_fraction": node.path_coverage,
                "predecessor_consistency": node.predecessor_consistency,
                "successor_consistency": node.successor_consistency,
                "local_consistency": node.local_path_consistency,
                "consistency_score": node.consistency_score,
                "status": _node_status(node),
            })
    with edge_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "source", "target", "haplotype_count", "path_support_fraction",
            "source_transition_agreement", "consistency_score",
        ])
        writer.writeheader()
        for edge in graph.edges:
            writer.writerow({
                "source": edge.source,
                "target": edge.target,
                "haplotype_count": len(edge.haplotypes),
                "path_support_fraction": edge.path_coverage,
                "source_transition_agreement": edge.transition_consistency,
                "consistency_score": edge.consistency_score,
            })
    with path_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "assembly_id", "haplotype_label", "path_length", "valid_transition_count",
            "transition_fraction", "mean_node_consistency", "mean_edge_consistency",
            "rescued_node_count", "questionable_node_count", "consistency_score",
        ])
        writer.writeheader()
        for path in graph.paths:
            writer.writerow({
                "assembly_id": path.assembly_id,
                "haplotype_label": path.haplotype_label,
                "path_length": len(path.node_ids),
                "valid_transition_count": path.valid_transition_count,
                "transition_fraction": path.transition_fraction,
                "mean_node_consistency": path.mean_node_consistency,
                "mean_edge_consistency": path.mean_edge_consistency,
                "rescued_node_count": path.rescued_node_count,
                "questionable_node_count": path.questionable_node_count,
                "consistency_score": path.consistency_score,
            })
    graph_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return {
        "node_consistency_csv": str(node_path.resolve()),
        "edge_consistency_csv": str(edge_path.resolve()),
        "path_consistency_csv": str(path_path.resolve()),
        "graph_consistency_json": str(graph_path.resolve()),
    }
