from __future__ import annotations

import argparse
from pathlib import Path

from align import run_panaligner
from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, read_json, resolve_binary, write_json
from parse_gaf import best_alignment
from visualize import visualize_alignment


def evaluate_panaligner_workflow(
    test_manifest_path: Path,
    train_graph_manifest_path: Path,
    panaligner_bin: Path,
    threads: int,
    output_dir: Path,
) -> dict:
    test_manifest = read_json(test_manifest_path)
    graph_manifest = read_json(train_graph_manifest_path)
    evaluation_root = ensure_dir(output_dir)
    alignments_root = ensure_dir(evaluation_root / "alignments")

    records: list[dict] = []
    per_gene_summary: dict[str, dict[str, float | int]] = {}

    for gene, payload in test_manifest["genes"].items():
        combined_graph = Path(graph_manifest["graphs"][gene]["combined_graph"])
        per_gene_summary[gene] = {
            "sequence_count": 0,
            "aligned_sequences": 0,
            "mean_identity": 0.0,
            "mean_coverage": 0.0,
            "mean_mapq": 0.0,
            "mean_alignment_score": 0.0,
        }
        first_visualized = False

        for bucket_name in ("healthy", "unhealthy"):
            for entry in payload[bucket_name]:
                sequence_path = Path(entry["path"])
                sequence_output_dir = ensure_dir(alignments_root / gene.lower() / entry["id"])
                gaf_path = sequence_output_dir / "combined.gaf"
                results = run_panaligner(combined_graph, sequence_path, gaf_path, panaligner_bin, threads)
                best = best_alignment(results)

                record = {
                    "gene": gene,
                    "sequence_id": entry["id"],
                    "bucket": bucket_name,
                    "query_fasta": str(sequence_path.resolve()),
                    "gaf_path": str(gaf_path.resolve()),
                    "alignment_found": best is not None,
                    "alignment_score": best.alignment_score if best else 0.0,
                    "normalized_score": best.normalized_score if best else 0.0,
                    "identity": best.identity if best else 0.0,
                    "coverage": best.coverage if best else 0.0,
                    "mapping_quality": best.mapping_quality if best else 0,
                    "path_span": best.path_span if best else 0,
                    "matched_nodes": len(best.traversed_nodes) if best else 0,
                }
                records.append(record)

                summary = per_gene_summary[gene]
                summary["sequence_count"] += 1
                if best is not None:
                    summary["aligned_sequences"] += 1
                    summary["mean_identity"] += best.identity
                    summary["mean_coverage"] += best.coverage
                    summary["mean_mapq"] += best.mapping_quality
                    summary["mean_alignment_score"] += best.alignment_score
                    if not first_visualized:
                        visualize_alignment(
                            combined_graph,
                            gaf_path,
                            evaluation_root / f"{gene.lower()}.example_alignment.png",
                        )
                        first_visualized = True

    for gene, summary in per_gene_summary.items():
        aligned = max(1, int(summary["aligned_sequences"]))
        if summary["aligned_sequences"]:
            summary["mean_identity"] = float(summary["mean_identity"]) / aligned
            summary["mean_coverage"] = float(summary["mean_coverage"]) / aligned
            summary["mean_mapq"] = float(summary["mean_mapq"]) / aligned
            summary["mean_alignment_score"] = float(summary["mean_alignment_score"]) / aligned

    aligned_records = [record for record in records if record["alignment_found"]]
    overall = {
        "query_count": len(records),
        "aligned_query_count": len(aligned_records),
        "alignment_rate": len(aligned_records) / max(1, len(records)),
        "mean_identity": sum(record["identity"] for record in aligned_records) / max(1, len(aligned_records)),
        "mean_coverage": sum(record["coverage"] for record in aligned_records) / max(1, len(aligned_records)),
        "mean_mapq": sum(record["mapping_quality"] for record in aligned_records) / max(1, len(aligned_records)),
        "mean_alignment_score": sum(record["alignment_score"] for record in aligned_records) / max(1, len(aligned_records)),
    }

    metrics = {
        "overall": overall,
        "per_gene": per_gene_summary,
    }
    write_json(evaluation_root / "alignment_metrics.json", metrics)
    write_json(evaluation_root / "alignment_records.json", records)

    lines = [
        "PanAligner Reproduction Evaluation",
        "=================================",
        "This evaluation summarizes sequence-to-graph alignment behavior on held-out test sequences.",
        "",
        "Overall metrics",
        f"- query count: {overall['query_count']}",
        f"- aligned query count: {overall['aligned_query_count']}",
        f"- alignment rate: {overall['alignment_rate']:.4f}",
        f"- mean identity: {overall['mean_identity']:.4f}",
        f"- mean coverage: {overall['mean_coverage']:.4f}",
        f"- mean MAPQ: {overall['mean_mapq']:.4f}",
        f"- mean alignment score: {overall['mean_alignment_score']:.4f}",
        "",
        "Per-gene metrics",
    ]
    for gene, summary in per_gene_summary.items():
        lines.extend(
            [
                f"{gene}",
                f"- query count: {summary['sequence_count']}",
                f"- aligned query count: {summary['aligned_sequences']}",
                f"- mean identity: {summary['mean_identity']:.4f}",
                f"- mean coverage: {summary['mean_coverage']:.4f}",
                f"- mean MAPQ: {summary['mean_mapq']:.4f}",
                f"- mean alignment score: {summary['mean_alignment_score']:.4f}",
                "",
            ]
        )
    (evaluation_root / "evaluation_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate the PanAligner reproduction workflow on held-out graph alignment queries.")
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json", help="Path to the test manifest JSON.")
    parser.add_argument("--train-graph-manifest", type=Path, default=METADATA_DIR / "train_graph_manifest.json", help="Path to the train graph manifest JSON.")
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument("--panaligner-bin", type=str, default=None, help="Path to PanAligner binary.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "evaluation", help="Output directory for evaluation artifacts.")
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    metrics = evaluate_panaligner_workflow(
        args.test_manifest.resolve(),
        args.train_graph_manifest.resolve(),
        panaligner_bin,
        args.threads,
        args.output_dir.resolve(),
    )
    print(f"Evaluation complete. Alignment rate={metrics['overall']['alignment_rate']:.4f}")


if __name__ == "__main__":
    main()
