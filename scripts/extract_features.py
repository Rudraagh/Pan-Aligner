from __future__ import annotations

import argparse
import csv
from pathlib import Path

from align import run_panaligner
from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, read_json, resolve_binary, write_json
from parse_gaf import AlignmentResult, best_alignment
from visualize import parse_gfa


FEATURE_COLUMNS = [
    "sequence_id",
    "split",
    "gene",
    "label_name",
    "label",
    "alignment_score",
    "mapq",
    "alignment_length",
    "query_coverage",
    "matched_graph_nodes",
    "traversed_graph_paths",
    "healthy_overlap_score",
    "unhealthy_overlap_score",
    "healthy_unhealthy_score_ratio",
    "mismatches",
    "graph_branches_traversed",
    "graph_path_complexity",
    "identity_percentage",
    "path_traversal_frequency",
    "node_revisit_count",
]


def label_to_int(label_name: str) -> int:
    return 1 if label_name.upper() == "UNHEALTHY" else 0


def safe_alignment_value(result: AlignmentResult | None, attribute: str, default: float = 0.0) -> float:
    if result is None:
        return default
    return float(getattr(result, attribute))


def get_mismatch_count(result: AlignmentResult | None) -> int:
    if result is None:
        return 0
    try:
        return int(result.tags.get("NM", 0))
    except ValueError:
        return 0


def compute_graph_path_features(graph, result: AlignmentResult | None) -> dict[str, float]:
    if result is None:
        return {
            "matched_graph_nodes": 0,
            "traversed_graph_paths": 0,
            "graph_branches_traversed": 0,
            "graph_path_complexity": 0.0,
            "path_traversal_frequency": 0.0,
            "node_revisit_count": 0,
        }

    traversed_nodes = result.traversed_nodes
    unique_nodes = list(dict.fromkeys(traversed_nodes))
    revisit_count = max(0, len(traversed_nodes) - len(unique_nodes))
    branching_nodes = [
        node
        for node in unique_nodes
        if graph.has_node(node) and (graph.out_degree(node) > 1 or graph.in_degree(node) > 1)
    ]
    if unique_nodes:
        average_degree = sum((graph.out_degree(node) + graph.in_degree(node)) for node in unique_nodes if graph.has_node(node)) / len(unique_nodes)
    else:
        average_degree = 0.0

    return {
        "matched_graph_nodes": len(unique_nodes),
        "traversed_graph_paths": len(traversed_nodes),
        "graph_branches_traversed": len(branching_nodes),
        "graph_path_complexity": average_degree,
        "path_traversal_frequency": len(traversed_nodes) / max(1, len(unique_nodes)),
        "node_revisit_count": revisit_count,
    }


def extract_sequence_features(
    sequence_entry: dict,
    split_name: str,
    gene: str,
    graph_manifest: dict,
    panaligner_bin: Path,
    threads: int,
    output_dir: Path,
    graph_cache: dict[str, object],
) -> dict:
    gene_graphs = graph_manifest["graphs"][gene]
    combined_graph_path = Path(gene_graphs["combined_graph"])
    healthy_graph_path = Path(gene_graphs["healthy_graph"])
    unhealthy_graph_path = Path(gene_graphs["unhealthy_graph"])
    sequence_path = Path(sequence_entry["path"])

    sequence_output_dir = ensure_dir(output_dir / split_name / gene.lower() / sequence_entry["id"])
    alignments = {}

    for label_name, graph_path in (
        ("healthy", healthy_graph_path),
        ("unhealthy", unhealthy_graph_path),
        ("combined", combined_graph_path),
    ):
        gaf_path = sequence_output_dir / f"{label_name}.gaf"
        alignments[label_name] = best_alignment(run_panaligner(graph_path, sequence_path, gaf_path, panaligner_bin, threads))

    combined_graph_key = str(combined_graph_path.resolve())
    if combined_graph_key not in graph_cache:
        graph_cache[combined_graph_key] = parse_gfa(combined_graph_path)
    combined_graph = graph_cache[combined_graph_key]

    combined_alignment = alignments["combined"]
    graph_features = compute_graph_path_features(combined_graph, combined_alignment)
    healthy_overlap_score = safe_alignment_value(alignments["healthy"], "normalized_score")
    unhealthy_overlap_score = safe_alignment_value(alignments["unhealthy"], "normalized_score")

    return {
        "sequence_id": sequence_entry["id"],
        "split": split_name,
        "gene": gene,
        "label_name": sequence_entry["label"],
        "label": label_to_int(sequence_entry["label"]),
        "alignment_score": safe_alignment_value(combined_alignment, "alignment_score"),
        "mapq": safe_alignment_value(combined_alignment, "mapping_quality"),
        "alignment_length": safe_alignment_value(combined_alignment, "alignment_block_length"),
        "query_coverage": safe_alignment_value(combined_alignment, "coverage"),
        "matched_graph_nodes": graph_features["matched_graph_nodes"],
        "traversed_graph_paths": graph_features["traversed_graph_paths"],
        "healthy_overlap_score": healthy_overlap_score,
        "unhealthy_overlap_score": unhealthy_overlap_score,
        "healthy_unhealthy_score_ratio": healthy_overlap_score / max(unhealthy_overlap_score, 1e-9),
        "mismatches": get_mismatch_count(combined_alignment),
        "graph_branches_traversed": graph_features["graph_branches_traversed"],
        "graph_path_complexity": graph_features["graph_path_complexity"],
        "identity_percentage": safe_alignment_value(combined_alignment, "identity") * 100.0,
        "path_traversal_frequency": graph_features["path_traversal_frequency"],
        "node_revisit_count": graph_features["node_revisit_count"],
    }


def generate_feature_dataset(
    train_manifest_path: Path,
    test_manifest_path: Path,
    train_graph_manifest_path: Path,
    panaligner_bin: Path,
    threads: int,
    output_dir: Path,
) -> dict:
    train_manifest = read_json(train_manifest_path)
    test_manifest = read_json(test_manifest_path)
    graph_manifest = read_json(train_graph_manifest_path)
    ml_output_dir = ensure_dir(output_dir)
    alignments_dir = ensure_dir(ml_output_dir / "alignments")

    rows: list[dict] = []
    graph_cache: dict[str, object] = {}

    for split_name, manifest in (("train", train_manifest), ("test", test_manifest)):
        for gene, payload in manifest["genes"].items():
            for label_key in ("healthy", "unhealthy"):
                for sequence_entry in payload[label_key]:
                    rows.append(
                        extract_sequence_features(
                            sequence_entry,
                            split_name,
                            gene,
                            graph_manifest,
                            panaligner_bin,
                            threads,
                            alignments_dir,
                            graph_cache,
                        )
                    )

    features_path = ml_output_dir / "features.csv"
    with features_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FEATURE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    for split_name in ("train", "test"):
        split_rows = [row for row in rows if row["split"] == split_name]
        split_path = ml_output_dir / f"{split_name}_features.csv"
        with split_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FEATURE_COLUMNS)
            writer.writeheader()
            writer.writerows(split_rows)

    summary = {
        "feature_columns": FEATURE_COLUMNS,
        "total_rows": len(rows),
        "train_rows": sum(1 for row in rows if row["split"] == "train"),
        "test_rows": sum(1 for row in rows if row["split"] == "test"),
    }
    write_json(ml_output_dir / "feature_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract ML-ready features from PanAligner graph-alignment outputs.")
    parser.add_argument("--train-manifest", type=Path, default=METADATA_DIR / "train_manifest.json", help="Path to train manifest JSON.")
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json", help="Path to test manifest JSON.")
    parser.add_argument(
        "--train-graph-manifest",
        type=Path,
        default=METADATA_DIR / "train_graph_manifest.json",
        help="Path to train graph manifest JSON.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument(
        "--panaligner-bin",
        type=str,
        default=None,
        help="Path to the PanAligner binary. Defaults to ./PanAligner/PanAligner if present.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUTS_DIR / "ml",
        help="Directory where extracted features and alignment artifacts will be written.",
    )
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    summary = generate_feature_dataset(
        args.train_manifest.resolve(),
        args.test_manifest.resolve(),
        args.train_graph_manifest.resolve(),
        panaligner_bin,
        args.threads,
        args.output_dir.resolve(),
    )
    print(f"Feature extraction complete. Rows={summary['total_rows']}")
    print(f"Features written to {(args.output_dir.resolve() / 'features.csv')}")


if __name__ == "__main__":
    main()

