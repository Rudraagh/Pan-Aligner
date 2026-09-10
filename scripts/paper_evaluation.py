from __future__ import annotations

import argparse
from pathlib import Path
import csv
import shutil

import matplotlib.pyplot as plt

from align import run_minigraph, run_panaligner
from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, read_json, resolve_binary, write_json
from parse_gaf import best_alignment
from visualize import visualize_alignment
from visualize import graph_statistics, parse_gfa


def _safe_mean(values: list[float]) -> float:
    return sum(values) / max(1, len(values))


def summarize_records(records: list[dict]) -> dict[str, float | int]:
    aligned = [record for record in records if record["alignment_found"]]
    return {
        "query_count": len(records),
        "aligned_query_count": len(aligned),
        "unaligned_query_count": len(records) - len(aligned),
        "alignment_rate": len(aligned) / max(1, len(records)),
        "unaligned_rate": (len(records) - len(aligned)) / max(1, len(records)),
        "mean_identity": _safe_mean([float(record["identity"]) for record in aligned]),
        "mean_coverage": _safe_mean([float(record["coverage"]) for record in aligned]),
        "mean_mapq": _safe_mean([float(record["mapping_quality"]) for record in aligned]),
        "mean_alignment_score": _safe_mean([float(record["alignment_score"]) for record in aligned]),
        "mean_normalized_score": _safe_mean([float(record["normalized_score"]) for record in aligned]),
        "mean_path_span": _safe_mean([float(record["path_span"]) for record in aligned]),
        "mean_matched_nodes": _safe_mean([float(record["matched_nodes"]) for record in aligned]),
    }


def build_mq_cutoff_sweep(records: list[dict], cutoffs: list[int] | None = None) -> list[dict[str, float | int]]:
    thresholds = cutoffs or list(range(0, 61))
    total_queries = len(records)
    sweep: list[dict[str, float | int]] = []
    for cutoff in thresholds:
        kept = [
            record
            for record in records
            if record["alignment_found"] and int(record["mapping_quality"]) >= cutoff
        ]
        sweep.append(
            {
                "mapq_cutoff": cutoff,
                "retained_alignment_count": len(kept),
                "retained_alignment_fraction": len(kept) / max(1, total_queries),
                "mean_identity": _safe_mean([float(record["identity"]) for record in kept]),
                "mean_coverage": _safe_mean([float(record["coverage"]) for record in kept]),
                "mean_alignment_score": _safe_mean([float(record["alignment_score"]) for record in kept]),
                "mean_normalized_score": _safe_mean([float(record["normalized_score"]) for record in kept]),
            }
        )
    return sweep


def draw_mq_sweep_plot(sweep: list[dict[str, float | int]], output_path: Path, title: str) -> None:
    ensure_dir(output_path.parent)
    cutoffs = [int(item["mapq_cutoff"]) for item in sweep]
    retained_fraction = [float(item["retained_alignment_fraction"]) for item in sweep]
    mean_identity = [float(item["mean_identity"]) for item in sweep]
    mean_coverage = [float(item["mean_coverage"]) for item in sweep]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(cutoffs, retained_fraction, color="#2A9D8F", linewidth=2)
    axes[0].set_ylabel("Aligned Fraction")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].grid(alpha=0.25)
    axes[0].set_title(title)

    axes[1].plot(cutoffs, mean_identity, color="#264653", linewidth=2, label="Mean Identity")
    axes[1].plot(cutoffs, mean_coverage, color="#E76F51", linewidth=2, label="Mean Coverage")
    axes[1].set_xlabel("MAPQ Cutoff")
    axes[1].set_ylabel("Metric Value")
    axes[1].set_ylim(0.0, 1.05)
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_csv_rows(output_path: Path, rows: list[dict]) -> None:
    ensure_dir(output_path.parent)
    if not rows:
        output_path.write_text("", encoding="utf-8")
        return
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_local_graph_properties(graph_manifest: dict) -> list[dict[str, int | str | bool]]:
    rows: list[dict[str, int | str | bool]] = []
    for gene, entry in sorted(graph_manifest["graphs"].items()):
        combined_graph = Path(entry["combined_graph"])
        stats = graph_statistics(parse_gfa(combined_graph))
        rows.append(
            {
                "gene": gene,
                "graph_path": str(combined_graph.resolve()),
                "node_count": stats["node_count"],
                "edge_count": stats["edge_count"],
                "weakly_connected_components": stats["weakly_connected_components"],
                "strongly_connected_components": stats["strongly_connected_components"],
                "largest_component_size": stats["largest_component_size"],
                "largest_scc_size": stats["largest_scc_size"],
                "cycle_detected": stats["cycle_detected"],
                "self_loops": stats["self_loops"],
            }
        )
    return rows


def write_dataset_specific_report(
    evaluation_root: Path,
    overall: dict,
    per_gene_summary: dict,
    per_bucket_summary: dict,
    graph_rows: list[dict],
) -> Path:
    lines = [
        "Dataset-Specific Reproduction Report",
        "===================================",
        "This report captures the paper-style outputs that are feasible with the local APP/PSEN1/PSEN2 dataset.",
        "",
        "What is included here",
        "- Local combined-graph structural properties per gene.",
        "- Held-out alignment summary across all test sequences.",
        "- Per-gene and per-bucket alignment behavior.",
        "- Mapping-quality cutoff sweep artifacts and plots.",
        "",
        "Overall local evaluation",
        f"- query count: {overall['query_count']}",
        f"- aligned query count: {overall['aligned_query_count']}",
        f"- unaligned query count: {overall['unaligned_query_count']}",
        f"- alignment rate: {overall['alignment_rate']:.4f}",
        f"- mean identity: {overall['mean_identity']:.4f}",
        f"- mean coverage: {overall['mean_coverage']:.4f}",
        f"- mean MAPQ: {overall['mean_mapq']:.4f}",
        f"- mean alignment score: {overall['mean_alignment_score']:.4f}",
        "",
        "Local graph properties",
    ]
    for row in graph_rows:
        lines.extend(
            [
                f"{row['gene']}",
                f"- nodes: {row['node_count']}",
                f"- edges: {row['edge_count']}",
                f"- weak components: {row['weakly_connected_components']}",
                f"- SCCs: {row['strongly_connected_components']}",
                f"- cycle detected: {row['cycle_detected']}",
                "",
            ]
        )
    lines.append("Per-bucket alignment summary")
    for bucket, summary in sorted(per_bucket_summary.items()):
        lines.extend(
            [
                f"{bucket}",
                f"- query count: {summary['query_count']}",
                f"- alignment rate: {summary['alignment_rate']:.4f}",
                f"- mean identity: {summary['mean_identity']:.4f}",
                f"- mean coverage: {summary['mean_coverage']:.4f}",
                f"- mean MAPQ: {summary['mean_mapq']:.4f}",
                "",
            ]
        )
    lines.extend(
        [
            "Important limitation",
            "These outputs are dataset-specific approximations, not the paper's 10H/40H/80H/95H multi-aligner benchmark tables.",
        ]
    )
    report_path = evaluation_root / "dataset_specific_reproduction_report.txt"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def is_easy_minigraph_alignment(results: list) -> bool:
    if len(results) != 1:
        return False
    best = best_alignment(results)
    if best is None:
        return False
    return best.mapping_quality == 60 and len(best.traversed_nodes) == 1 and best.coverage >= 0.9


def run_alignment_with_local_hybrid(
    graph_gfa: Path,
    query_fasta: Path,
    output_gaf: Path,
    panaligner_bin: Path,
    minigraph_bin: Path | None,
    threads: int,
) -> tuple[object | None, dict]:
    if minigraph_bin is None:
        panaligner_results = run_panaligner(graph_gfa, query_fasta, output_gaf, panaligner_bin, threads)
        return best_alignment(panaligner_results), {
            "hybrid_used": False,
            "selected_aligner": "PanAligner",
            "minigraph_easy_case": False,
            "minigraph_alignment_count": 0,
        }

    minigraph_output = output_gaf.with_name(output_gaf.stem + ".minigraph.gaf")
    minigraph_results = run_minigraph(graph_gfa, query_fasta, minigraph_output, minigraph_bin, threads)
    if is_easy_minigraph_alignment(minigraph_results):
        shutil.copyfile(minigraph_output, output_gaf)
        return best_alignment(minigraph_results), {
            "hybrid_used": True,
            "selected_aligner": "Minigraph",
            "minigraph_easy_case": True,
            "minigraph_alignment_count": len(minigraph_results),
        }

    panaligner_output = output_gaf.with_name(output_gaf.stem + ".panaligner.gaf")
    panaligner_results = run_panaligner(graph_gfa, query_fasta, panaligner_output, panaligner_bin, threads)
    shutil.copyfile(panaligner_output, output_gaf)
    return best_alignment(panaligner_results), {
        "hybrid_used": True,
        "selected_aligner": "PanAligner",
        "minigraph_easy_case": False,
        "minigraph_alignment_count": len(minigraph_results),
    }


def evaluate_panaligner_workflow(
    test_manifest_path: Path,
    train_graph_manifest_path: Path,
    panaligner_bin: Path,
    minigraph_bin: Path | None,
    threads: int,
    output_dir: Path,
) -> dict:
    test_manifest = read_json(test_manifest_path)
    graph_manifest = read_json(train_graph_manifest_path)
    evaluation_root = ensure_dir(output_dir)
    alignments_root = ensure_dir(evaluation_root / "alignments")

    records: list[dict] = []
    per_gene_summary: dict[str, dict[str, float | int]] = {}
    hybrid_summary = {
        "enabled": minigraph_bin is not None,
        "minigraph_selected_count": 0,
        "panaligner_selected_count": 0,
        "easy_case_count": 0,
    }

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
                best, hybrid_info = run_alignment_with_local_hybrid(
                    combined_graph,
                    sequence_path,
                    gaf_path,
                    panaligner_bin,
                    minigraph_bin,
                    threads,
                )

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
                    "selected_aligner": hybrid_info["selected_aligner"],
                    "hybrid_used": hybrid_info["hybrid_used"],
                    "minigraph_easy_case": hybrid_info["minigraph_easy_case"],
                    "minigraph_alignment_count": hybrid_info["minigraph_alignment_count"],
                }
                records.append(record)
                if hybrid_info["selected_aligner"] == "Minigraph":
                    hybrid_summary["minigraph_selected_count"] += 1
                else:
                    hybrid_summary["panaligner_selected_count"] += 1
                if hybrid_info["minigraph_easy_case"]:
                    hybrid_summary["easy_case_count"] += 1

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

    overall = summarize_records(records)
    per_bucket_summary = {
        bucket_name: summarize_records([record for record in records if record["bucket"] == bucket_name])
        for bucket_name in ("healthy", "unhealthy")
    }
    mq_cutoff_sweep = build_mq_cutoff_sweep(records)
    per_gene_mq_cutoff_sweep = {
        gene: build_mq_cutoff_sweep([record for record in records if record["gene"] == gene])
        for gene in sorted(per_gene_summary)
    }
    graph_rows = build_local_graph_properties(graph_manifest)

    metrics = {
        "overall": overall,
        "per_gene": per_gene_summary,
        "per_bucket": per_bucket_summary,
        "mq_cutoff_sweep": mq_cutoff_sweep,
        "per_gene_mq_cutoff_sweep": per_gene_mq_cutoff_sweep,
        "local_graph_properties": graph_rows,
        "hybrid_summary": hybrid_summary,
    }
    write_json(evaluation_root / "alignment_metrics.json", metrics)
    write_json(evaluation_root / "alignment_records.json", records)
    write_json(evaluation_root / "paper_style_results.json", metrics)
    write_csv_rows(evaluation_root / "mq_cutoff_sweep.csv", mq_cutoff_sweep)
    write_csv_rows(evaluation_root / "local_graph_properties.csv", graph_rows)
    draw_mq_sweep_plot(
        mq_cutoff_sweep,
        evaluation_root / "mq_cutoff_sweep.png",
        "Local Dataset MQ Cutoff Sweep",
    )
    for gene, sweep in per_gene_mq_cutoff_sweep.items():
        draw_mq_sweep_plot(
            sweep,
            evaluation_root / f"{gene.lower()}.mq_cutoff_sweep.png",
            f"{gene} MQ Cutoff Sweep",
        )
    dataset_report_path = write_dataset_specific_report(
        evaluation_root,
        overall,
        per_gene_summary,
        per_bucket_summary,
        graph_rows,
    )

    lines = [
        "PanAligner Reproduction Evaluation",
        "=================================",
        "This evaluation summarizes sequence-to-graph alignment behavior on held-out test sequences.",
        "",
        "Overall metrics",
        f"- query count: {overall['query_count']}",
        f"- aligned query count: {overall['aligned_query_count']}",
        f"- unaligned query count: {overall['unaligned_query_count']}",
        f"- alignment rate: {overall['alignment_rate']:.4f}",
        f"- mean identity: {overall['mean_identity']:.4f}",
        f"- mean coverage: {overall['mean_coverage']:.4f}",
        f"- mean MAPQ: {overall['mean_mapq']:.4f}",
        f"- mean alignment score: {overall['mean_alignment_score']:.4f}",
        f"- mean normalized score: {overall['mean_normalized_score']:.4f}",
        "",
        "Hybrid approximation summary",
        f"- enabled: {hybrid_summary['enabled']}",
        f"- minigraph-selected queries: {hybrid_summary['minigraph_selected_count']}",
        f"- PanAligner-selected queries: {hybrid_summary['panaligner_selected_count']}",
        f"- easy-case queries: {hybrid_summary['easy_case_count']}",
        "",
        "Per-bucket metrics",
    ]
    for bucket_name, summary in per_bucket_summary.items():
        lines.extend(
            [
                f"{bucket_name}",
                f"- query count: {summary['query_count']}",
                f"- aligned query count: {summary['aligned_query_count']}",
                f"- unaligned query count: {summary['unaligned_query_count']}",
                f"- alignment rate: {summary['alignment_rate']:.4f}",
                f"- mean identity: {summary['mean_identity']:.4f}",
                f"- mean coverage: {summary['mean_coverage']:.4f}",
                f"- mean MAPQ: {summary['mean_mapq']:.4f}",
                "",
            ]
        )
    lines.extend(
        [
        "Per-gene metrics",
        ]
    )
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
    lines.extend(
        [
            "Additional paper-style local artifacts",
            f"- mq cutoff sweep csv: {(evaluation_root / 'mq_cutoff_sweep.csv').resolve()}",
            f"- mq cutoff sweep plot: {(evaluation_root / 'mq_cutoff_sweep.png').resolve()}",
            f"- local graph properties csv: {(evaluation_root / 'local_graph_properties.csv').resolve()}",
            f"- local reproduction report: {dataset_report_path.resolve()}",
            "",
            "Note",
            "These artifacts approximate the paper's evaluation style with the local dataset, but they do not replace the paper's full external benchmark suite.",
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
    parser.add_argument("--minigraph-bin", type=str, default=None, help="Optional minigraph binary for local hybrid approximation.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "evaluation", help="Output directory for evaluation artifacts.")
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph")]) if args.minigraph_bin else None
    metrics = evaluate_panaligner_workflow(
        args.test_manifest.resolve(),
        args.train_graph_manifest.resolve(),
        panaligner_bin,
        minigraph_bin,
        args.threads,
        args.output_dir.resolve(),
    )
    print(f"Evaluation complete. Alignment rate={metrics['overall']['alignment_rate']:.4f}")


if __name__ == "__main__":
    main()
