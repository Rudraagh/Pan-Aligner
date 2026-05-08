from __future__ import annotations

import argparse
from pathlib import Path

from build_graph import build_graphs
from common import METADATA_DIR, OUTPUTS_DIR, ROOT, resolve_binary
from predictor import predict_health_state, prepare_query_fasta
from preprocess import preprocess_fasta_files
from visualize import visualize_alignment, visualize_graph


def default_fasta_inputs() -> list[Path]:
    return sorted(ROOT.glob("*_combined.fasta"))


def run_pipeline(
    input_fastas: list[Path],
    minigraph_bin: Path,
    panaligner_bin: Path,
    threads: int,
    query_fasta: Path | None,
    query_sequence: str | None,
    gene: str | None,
) -> dict | None:
    preprocess_fasta_files(input_fastas)
    graph_manifest = build_graphs(METADATA_DIR / "preprocess_manifest.json", minigraph_bin, threads)

    for gene_name, graph_entry in graph_manifest["graphs"].items():
        visualize_graph(Path(graph_entry["combined_graph"]), OUTPUTS_DIR / "graphs" / f"{gene_name.lower()}.combined")
        visualize_graph(Path(graph_entry["healthy_graph"]), OUTPUTS_DIR / "graphs" / f"{gene_name.lower()}.healthy")
        visualize_graph(Path(graph_entry["unhealthy_graph"]), OUTPUTS_DIR / "graphs" / f"{gene_name.lower()}.unhealthy")

    if query_fasta or query_sequence:
        prepared_query = prepare_query_fasta(query_sequence, query_fasta, ROOT / "data" / "queries" / "query.fa")
        prediction = predict_health_state(prepared_query, METADATA_DIR / "graph_manifest.json", panaligner_bin, threads, gene)
        selected_gene = prediction["selected_gene"]
        combined_gaf = OUTPUTS_DIR / "alignments" / f"{selected_gene.lower()}.combined.gaf"
        if combined_gaf.exists():
            visualize_alignment(Path(graph_manifest["graphs"][selected_gene]["combined_graph"]), combined_gaf, OUTPUTS_DIR / "alignments" / f"{selected_gene.lower()}.alignment.png")
        return prediction
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the complete Alzheimer’s pangenome graph prediction pipeline.")
    parser.add_argument(
        "--input-fastas",
        nargs="+",
        type=Path,
        default=default_fasta_inputs(),
        help="Input gene FASTA files. Defaults to all '*_combined.fasta' files in the project root.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Thread count for minigraph and PanAligner.")
    parser.add_argument("--query-fasta", type=Path, default=None, help="Optional query FASTA.")
    parser.add_argument("--query-sequence", type=str, default=None, help="Optional inline query sequence.")
    parser.add_argument("--gene", type=str, default=None, help="Optional gene override.")
    parser.add_argument("--minigraph-bin", type=str, default=None, help="Path to minigraph. Defaults to ./minigraph/minigraph.")
    parser.add_argument("--panaligner-bin", type=str, default=None, help="Path to PanAligner. Defaults to ./PanAligner/PanAligner.")
    args = parser.parse_args()

    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph")])
    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    prediction = run_pipeline(
        [path.resolve() for path in args.input_fastas],
        minigraph_bin,
        panaligner_bin,
        args.threads,
        args.query_fasta.resolve() if args.query_fasta else None,
        args.query_sequence,
        args.gene,
    )

    print("Pipeline complete.")
    if prediction:
        print(f"Prediction: {prediction['prediction']}")
        print(f"Confidence: {prediction['confidence']:.4f}")
        print(f"Selected gene: {prediction['selected_gene']}")


if __name__ == "__main__":
    main()

