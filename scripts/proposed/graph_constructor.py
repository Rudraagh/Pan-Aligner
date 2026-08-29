from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

from common import ensure_dir, write_json

from .anchor_selection import (
    DEFAULT_ALPHA,
    DEFAULT_BETA,
    DEFAULT_GAMMA,
    anchor_summary,
    build_graph_anchor_candidates,
    select_informative_anchors,
    validate_anchor_weights,
    write_anchor_reports,
)
from .graph_consistency import calculate_graph_path_consistency, write_graph_consistency_reports
from .minimizer_engine import compute_minimizers
from .models import (
    AnchorCandidate,
    AssemblyInput,
    HaplotypePath,
    MinimizerRecord,
    ProposedGraph,
    ProposedGraphEdge,
    ProposedGraphNode,
)


def _records_from_payload(payload: dict) -> list[MinimizerRecord]:
    return [MinimizerRecord(**item) for item in payload["records"]]


def _mean(values: list[float]) -> float:
    return sum(values) / max(1, len(values))


def validate_proposed_graph(graph: ProposedGraph) -> list[str]:
    """Return path-integrity errors without mutating the graph."""

    node_ids = {node.node_id for node in graph.nodes}
    edge_pairs = {(edge.source, edge.target) for edge in graph.edges}
    errors: list[str] = []
    observed_paths: set[str] = set()

    for path in graph.paths:
        if path.assembly_id in observed_paths:
            errors.append(f"Duplicate haplotype path: {path.assembly_id}")
        observed_paths.add(path.assembly_id)
        if len(path.node_ids) != len(path.positions):
            errors.append(f"Path {path.assembly_id} has mismatched node and position counts.")
        if path.positions != sorted(path.positions):
            errors.append(f"Path {path.assembly_id} positions are not ordered.")
        for node_id in path.node_ids:
            if node_id not in node_ids:
                errors.append(f"Path {path.assembly_id} references missing node {node_id}.")
        for source, target in zip(path.node_ids, path.node_ids[1:]):
            if (source, target) not in edge_pairs:
                errors.append(f"Path {path.assembly_id} transition {source}->{target} is missing an edge.")
    return errors


def build_proposed_pangenome_graph(
    assemblies: list[AssemblyInput],
    *,
    k: int = 9,
    window_size: int = 7,
    min_node_support: int = 2,
    max_anchors: int = 12,
    min_anchor_spacing: float = 4.0,
    anchor_alpha: float = DEFAULT_ALPHA,
    anchor_beta: float = DEFAULT_BETA,
    anchor_gamma: float = DEFAULT_GAMMA,
    backend: str = "auto",
) -> ProposedGraph:
    if not assemblies:
        raise ValueError("At least one assembly is required.")

    gene_names = {item.gene for item in assemblies}
    if len(gene_names) != 1:
        raise ValueError("All assemblies passed to the proposed graph builder must belong to the same gene.")
    validate_anchor_weights(anchor_alpha, anchor_beta, anchor_gamma)

    minimizers_by_assembly: dict[str, list[MinimizerRecord]] = {}
    backend_used = "cpu"
    for assembly in assemblies:
        payload = compute_minimizers(assembly.sequence, k, window_size, assembly.assembly_id, backend=backend)
        backend_used = payload["backend_used"]
        minimizers_by_assembly[assembly.assembly_id] = _records_from_payload(payload)

    grouped_records: dict[str, list[MinimizerRecord]] = defaultdict(list)
    for records in minimizers_by_assembly.values():
        seen_in_assembly: set[str] = set()
        for record in records:
            if record.minimizer in seen_in_assembly:
                continue
            grouped_records[record.minimizer].append(record)
            seen_in_assembly.add(record.minimizer)

    nodes: list[ProposedGraphNode] = []
    node_lookup: dict[str, ProposedGraphNode] = {}
    for minimizer, records in sorted(grouped_records.items()):
        support = len({item.source_assembly for item in records})
        # Preserve every observed minimizer walk. Low-support nodes are retained as
        # path-specific rescue nodes instead of silently removing divergent sequence.
        is_rescued = support < min_node_support

        mean_position = _mean([float(item.position) for item in records])
        node = ProposedGraphNode(
            node_id=f"n{len(nodes) + 1}",
            sequence=minimizer,
            support=support,
            mean_position=mean_position,
            consistency_score=0.0,
            haplotypes=sorted({item.source_assembly for item in records}),
            is_rescued=is_rescued,
        )
        nodes.append(node)
        node_lookup[minimizer] = node

    edge_supports: dict[tuple[str, str], set[str]] = defaultdict(set)
    paths: list[HaplotypePath] = []
    for assembly in assemblies:
        records = minimizers_by_assembly[assembly.assembly_id]
        ordered_nodes: list[str] = []
        ordered_positions: list[int] = []
        for record in records:
            node_id = node_lookup[record.minimizer].node_id
            if not ordered_nodes or ordered_nodes[-1] != node_id:
                ordered_nodes.append(node_id)
                ordered_positions.append(record.position)
        paths.append(
            HaplotypePath(
                assembly_id=assembly.assembly_id,
                haplotype_label=assembly.haplotype_label,
                node_ids=ordered_nodes,
                positions=ordered_positions,
            )
        )
        for source, target in zip(ordered_nodes, ordered_nodes[1:]):
            edge_supports[(source, target)].add(assembly.assembly_id)

    edges: list[ProposedGraphEdge] = []
    for (source, target), haplotypes in sorted(edge_supports.items()):
        support = len(haplotypes)
        edges.append(
            ProposedGraphEdge(
                source=source,
                target=target,
                support=support,
                consistency_score=0.0,
                haplotypes=sorted(haplotypes),
            )
        )

    shared_nodes = [node for node in nodes if node.support > 1]
    haplotype_specific_nodes = [node for node in nodes if node.support == 1]
    rescued_nodes = [node for node in nodes if node.is_rescued]
    metadata = {
        "path_retention_policy": (
            "Every observed assembly minimizer walk is retained. Nodes below "
            "min_node_support are kept as rescued path-specific nodes."
        ),
        "haplotype_count": len(assemblies),
        "path_count": len(paths),
        "shared_node_count": len(shared_nodes),
        "haplotype_specific_node_count": len(haplotype_specific_nodes),
        "rescued_node_count": len(rescued_nodes),
        "rescued_private_node_count": len([node for node in rescued_nodes if node.support == 1]),
        "discarded": {
            "minimizer_records": 0,
            "nodes": 0,
            "reason": "Observed minimizers are retained to preserve haplotype paths.",
        },
    }

    graph = ProposedGraph(
        gene=next(iter(gene_names)),
        assembly_count=len(assemblies),
        parameters={
            "k": k,
            "window_size": window_size,
            "min_node_support": min_node_support,
            "max_anchors": max_anchors,
            "min_anchor_spacing": min_anchor_spacing,
            "anchor_alpha": anchor_alpha,
            "anchor_beta": anchor_beta,
            "anchor_gamma": anchor_gamma,
        },
        backend_used=backend_used,
        nodes=nodes,
        edges=edges,
        anchors=[],
        quality={},
        paths=paths,
        metadata=metadata,
    )
    consistency_summary = calculate_graph_path_consistency(graph)
    graph.anchor_candidates = build_graph_anchor_candidates(
        graph,
        alpha=anchor_alpha,
        beta=anchor_beta,
        gamma=anchor_gamma,
    )
    graph.anchors = select_informative_anchors(
        graph.anchor_candidates,
        max_anchors=max_anchors,
        min_spacing=min_anchor_spacing,
    )
    graph.metadata["anchor_scoring"] = anchor_summary(graph.anchor_candidates)
    anchor_haplotype_coverage = sorted({haplotype for anchor in graph.anchors for haplotype in anchor.haplotypes})
    outgoing_counts: dict[str, int] = defaultdict(int)
    for edge in edges:
        outgoing_counts[edge.source] += 1
    graph.quality = {
        "node_count": len(nodes),
        "edge_count": len(edges),
        "mean_node_consistency": consistency_summary["mean_node_consistency"],
        "mean_edge_consistency": consistency_summary["mean_edge_consistency"],
        "mean_path_consistency": consistency_summary["mean_path_consistency"],
        "graph_path_consistency": consistency_summary["graph_path_consistency"],
        "low_consistency_node_fraction": consistency_summary["low_consistency_node_fraction"],
        "low_consistency_edge_fraction": consistency_summary["low_consistency_edge_fraction"],
        "mean_node_support_ratio": _mean([node.support / max(1, len(assemblies)) for node in nodes]) if nodes else 0.0,
        "anchor_count": len(graph.anchors),
        "anchor_coverage_ratio": len(anchor_haplotype_coverage) / max(1, len(assemblies)),
        "branching_node_ratio": (
            len([node_id for node_id, count in outgoing_counts.items() if count > 1]) / max(1, len(nodes))
            if nodes
            else 0.0
        ),
    }
    graph.quality["overall_quality_score"] = (
        0.25 * graph.quality["mean_node_consistency"]
        + 0.20 * graph.quality["mean_edge_consistency"]
        + 0.15 * graph.quality["mean_path_consistency"]
        + 0.15 * graph.quality["mean_node_support_ratio"]
        + 0.15 * graph.quality["anchor_coverage_ratio"]
        + 0.10 * (1.0 - min(1.0, graph.quality["branching_node_ratio"]))
    )
    errors = validate_proposed_graph(graph)
    if errors:
        raise ValueError("Invalid proposed graph: " + "; ".join(errors))
    return graph


def write_proposed_gfa(graph: ProposedGraph, output_path: Path) -> Path:
    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("H\tVN:Z:1.0\n")
        for node in graph.nodes:
            handle.write(
                f"S\t{node.node_id}\t{node.sequence}\tRC:i:{node.support}\tCS:f:{node.consistency_score:.6f}\n"
            )
        for edge in graph.edges:
            handle.write(
                f"L\t{edge.source}\t+\t{edge.target}\t+\t0M\tRC:i:{edge.support}\tCS:f:{edge.consistency_score:.6f}\n"
            )
        for path in graph.paths:
            segments = ",".join(f"{node_id}+" for node_id in path.node_ids) or "*"
            handle.write(
                f"P\t{path.assembly_id}\t{segments}\t*\tHP:Z:{path.haplotype_label}\n"
            )
    return output_path


def write_proposed_graph_bundle(graph: ProposedGraph, output_dir: Path) -> dict:
    ensure_dir(output_dir)
    graph_json = output_dir / f"{graph.gene.lower()}.proposed_graph.json"
    graph_gfa = output_dir / f"{graph.gene.lower()}.proposed_graph.gfa"
    anchors_json = output_dir / f"{graph.gene.lower()}.selected_anchors.json"
    quality_json = output_dir / f"{graph.gene.lower()}.quality.json"

    write_json(
        graph_json,
        {
            "gene": graph.gene,
            "assembly_count": graph.assembly_count,
            "parameters": graph.parameters,
            "backend_used": graph.backend_used,
            "nodes": [asdict(node) for node in graph.nodes],
            "edges": [asdict(edge) for edge in graph.edges],
            "paths": [asdict(path) for path in graph.paths],
            "metadata": graph.metadata,
            "anchors": [asdict(anchor) for anchor in graph.anchors],
            "anchor_candidates": [asdict(anchor) for anchor in graph.anchor_candidates],
            "quality": graph.quality,
        },
    )
    write_json(anchors_json, [asdict(anchor) for anchor in graph.anchors])
    write_json(quality_json, graph.quality)
    write_proposed_gfa(graph, graph_gfa)
    consistency_reports = write_graph_consistency_reports(graph, output_dir)
    anchor_reports = write_anchor_reports(graph.anchor_candidates, output_dir)
    return {
        "graph_json": str(graph_json.resolve()),
        "graph_gfa": str(graph_gfa.resolve()),
        "anchors_json": str(anchors_json.resolve()),
        "quality_json": str(quality_json.resolve()),
        **consistency_reports,
        **anchor_reports,
    }
