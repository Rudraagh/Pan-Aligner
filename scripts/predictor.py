from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

from align import run_panaligner
from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, read_json, resolve_binary, run_command, write_json, write_wrapped_fasta
from parse_gaf import AlignmentResult, best_alignment, parse_gaf_file


def prepare_query_fasta(query_sequence: str | None, query_fasta: Path | None, output_path: Path) -> Path:
    if query_fasta:
        return query_fasta.resolve()
    if not query_sequence:
        raise ValueError("Provide either --query-fasta or --query-sequence.")
    write_wrapped_fasta([("query", query_sequence.upper())], output_path)
    return output_path.resolve()


def summarize_alignment(result: AlignmentResult | None) -> dict:
    if result is None:
        return {
            "available": False,
            "normalized_score": 0.0,
            "identity": 0.0,
            "coverage": 0.0,
            "mapping_quality": 0,
            "path": "",
            "traversed_nodes": [],
        }
    return {
        "available": True,
        "query_name": result.query_name,
        "normalized_score": result.normalized_score,
        "identity": result.identity,
        "coverage": result.coverage,
        "mapping_quality": result.mapping_quality,
        "path": result.path,
        "traversed_nodes": result.traversed_nodes,
        "path_span": result.path_span,
        "alignment_score": result.alignment_score,
    }


def composite_score(result: AlignmentResult | None) -> float:
    if result is None:
        return 0.0
    return (
        0.45 * result.identity
        + 0.35 * result.coverage
        + 0.15 * min(result.mapping_quality / 60.0, 1.0)
        + 0.05 * min(result.normalized_score, 1.0)
    )


def build_prediction_payload(
    query_fasta: Path,
    selected_gene: str,
    class_results: dict[str, AlignmentResult | None],
    gene_detection: dict,
) -> dict:
    healthy_score = composite_score(class_results["HEALTHY"])
    unhealthy_score = composite_score(class_results["UNHEALTHY"])
    predicted_label = "HEALTHY" if healthy_score >= unhealthy_score else "UNHEALTHY"
    score_gap = abs(healthy_score - unhealthy_score)
    total_score = healthy_score + unhealthy_score + 1e-9
    confidence = min(1.0, score_gap / total_score)
    tie_like = score_gap < 1e-6

    if tie_like:
        explanation = (
            f"Gene inference selected {selected_gene}, but the healthy and unhealthy graph alignments are effectively indistinguishable "
            f"(healthy={healthy_score:.4f}, unhealthy={unhealthy_score:.4f}). This suggests the current graph pair does not separate the query clearly."
        )
    else:
        explanation = (
            f"Gene inference selected {selected_gene}. "
            f"The best {predicted_label.lower()}-graph alignment had higher composite evidence "
            f"(healthy={healthy_score:.4f}, unhealthy={unhealthy_score:.4f}) using identity, coverage, MAPQ, and normalized alignment score."
        )

    prediction = {
        "query_fasta": str(query_fasta.resolve()),
        "selected_gene": selected_gene,
        "prediction": predicted_label,
        "confidence": confidence,
        "healthy_score": healthy_score,
        "unhealthy_score": unhealthy_score,
        "tie_like": tie_like,
        "gene_detection": gene_detection,
        "combined_alignment": summarize_alignment(class_results["COMBINED"]),
        "class_alignments": {
            "HEALTHY": summarize_alignment(class_results["HEALTHY"]),
            "UNHEALTHY": summarize_alignment(class_results["UNHEALTHY"]),
        },
        "matched_regions": {
            "predicted_class_nodes": summarize_alignment(class_results[predicted_label])["traversed_nodes"],
            "combined_graph_nodes": summarize_alignment(class_results["COMBINED"])["traversed_nodes"],
        },
        "explanation": explanation,
    }
    return prediction


def infer_gene(query_fasta: Path, graph_manifest: dict, panaligner_bin: Path, threads: int, output_dir: Path) -> tuple[str, dict]:
    per_gene_scores: dict[str, dict] = {}
    for gene, entry in graph_manifest["graphs"].items():
        gaf_path = output_dir / f"{gene.lower()}.combined.gaf"
        results = run_panaligner(Path(entry["combined_graph"]), query_fasta, gaf_path, panaligner_bin, threads)
        best = best_alignment(results)
        per_gene_scores[gene] = summarize_alignment(best)

    best_gene = max(
        per_gene_scores,
        key=lambda gene_name: (
            per_gene_scores[gene_name]["normalized_score"],
            per_gene_scores[gene_name]["identity"],
            per_gene_scores[gene_name]["coverage"],
            per_gene_scores[gene_name]["mapping_quality"],
        ),
    )
    return best_gene, per_gene_scores


def predict_health_state(
    query_fasta: Path,
    graph_manifest_path: Path,
    panaligner_bin: Path,
    threads: int,
    gene: str | None = None,
) -> dict:
    graph_manifest = read_json(graph_manifest_path)
    output_dir = ensure_dir(OUTPUTS_DIR / "alignments")

    selected_gene, gene_detection = infer_gene(query_fasta, graph_manifest, panaligner_bin, threads, output_dir) if gene is None else (gene.upper(), {})
    if selected_gene not in graph_manifest["graphs"]:
        raise ValueError(f"Unknown gene '{selected_gene}'. Available genes: {', '.join(sorted(graph_manifest['graphs']))}")

    gene_entry = graph_manifest["graphs"][selected_gene]
    class_results: dict[str, AlignmentResult | None] = {}

    for label, graph_key in (("HEALTHY", "healthy_graph"), ("UNHEALTHY", "unhealthy_graph"), ("COMBINED", "combined_graph")):
        gaf_path = output_dir / f"{selected_gene.lower()}.{label.lower()}.gaf"
        results = run_panaligner(Path(gene_entry[graph_key]), query_fasta, gaf_path, panaligner_bin, threads)
        class_results[label] = best_alignment(results)

    prediction = build_prediction_payload(query_fasta, selected_gene, class_results, gene_detection)

    prediction_path = OUTPUTS_DIR / "alignments" / "prediction.json"
    write_json(prediction_path, prediction)
    return prediction


def predict_from_existing_gafs(
    query_fasta: Path,
    selected_gene: str,
    healthy_gaf: Path,
    unhealthy_gaf: Path,
    combined_gaf: Path,
) -> dict:
    class_results = {
        "HEALTHY": best_alignment(parse_gaf_file(healthy_gaf.resolve())),
        "UNHEALTHY": best_alignment(parse_gaf_file(unhealthy_gaf.resolve())),
        "COMBINED": best_alignment(parse_gaf_file(combined_gaf.resolve())),
    }
    prediction = build_prediction_payload(query_fasta, selected_gene.upper(), class_results, {})
    prediction_path = OUTPUTS_DIR / "alignments" / "prediction.json"
    write_json(prediction_path, prediction)
    return prediction


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict healthy vs unhealthy status for a query using PanAligner class-specific graphs.")
    parser.add_argument("--query-fasta", type=Path, default=None, help="Existing query FASTA.")
    parser.add_argument("--query-sequence", type=str, default=None, help="Inline query sequence.")
    parser.add_argument("--gene", type=str, default=None, help="Optional gene override (APP, PSEN1, PSEN2).")
    parser.add_argument("--healthy-gaf", type=Path, default=None, help="Existing healthy-graph GAF for score-only mode.")
    parser.add_argument("--unhealthy-gaf", type=Path, default=None, help="Existing unhealthy-graph GAF for score-only mode.")
    parser.add_argument("--combined-gaf", type=Path, default=None, help="Existing combined-graph GAF for score-only mode.")
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument(
        "--graph-manifest",
        type=Path,
        default=METADATA_DIR / "graph_manifest.json",
        help="Graph manifest produced by build_graph.py.",
    )
    parser.add_argument(
        "--panaligner-bin",
        type=str,
        default=None,
        help="Path to the PanAligner binary. Defaults to ./PanAligner/PanAligner if present.",
    )
    args = parser.parse_args()

    query_fasta = prepare_query_fasta(
        args.query_sequence,
        args.query_fasta,
        ROOT / "data" / "queries" / "query.fa",
    )
    if args.healthy_gaf and args.unhealthy_gaf and args.combined_gaf and args.gene:
        prediction = predict_from_existing_gafs(
            query_fasta,
            args.gene,
            args.healthy_gaf,
            args.unhealthy_gaf,
            args.combined_gaf,
        )
    else:
        panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
        prediction = predict_health_state(query_fasta, args.graph_manifest.resolve(), panaligner_bin, args.threads, args.gene)
    print(f"Prediction: {prediction['prediction']}")
    print(f"Selected gene: {prediction['selected_gene']}")
    print(f"Confidence: {prediction['confidence']:.4f}")
    print(prediction["explanation"])


if __name__ == "__main__":
    main()
