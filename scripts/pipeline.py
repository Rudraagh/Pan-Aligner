from __future__ import annotations

import argparse
from pathlib import Path

from build_graph import build_graphs
from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, resolve_binary
from evaluate import evaluate_test_set
from extract_features import generate_feature_dataset
from predict_ml import predict_with_ml
from predictor import predict_health_state, prepare_query_fasta
from preprocess import preprocess_fasta_files
from split_dataset import split_dataset
from train_ml import train_and_evaluate_models
from visualize import visualize_alignment, visualize_graph


def default_fasta_inputs() -> list[Path]:
    return sorted(ROOT.glob("*_combined.fasta"))


def run_pipeline(
    input_fastas: list[Path],
    minigraph_bin: Path,
    panaligner_bin: Path,
    threads: int,
    split_seed: int,
    test_fraction: float,
    evaluate_split: bool,
    query_fasta: Path | None,
    query_sequence: str | None,
    gene: str | None,
    predict_ml_query: bool,
    ml_model_path: Path,
) -> dict | None:
    preprocess_fasta_files(input_fastas)
    split_dataset(METADATA_DIR / "preprocess_manifest.json", test_fraction, split_seed)
    train_graph_dir = ensure_dir(ROOT / "graphs" / "train")
    train_graph_manifest_path = METADATA_DIR / "train_graph_manifest.json"
    graph_manifest = build_graphs(
        METADATA_DIR / "train_manifest.json",
        minigraph_bin,
        threads,
        train_graph_dir,
        train_graph_manifest_path,
    )

    for gene_name, graph_entry in graph_manifest["graphs"].items():
        train_graph_output_dir = ensure_dir(OUTPUTS_DIR / "graphs" / "train")
        visualize_graph(Path(graph_entry["combined_graph"]), train_graph_output_dir / f"{gene_name.lower()}.combined")
        visualize_graph(Path(graph_entry["healthy_graph"]), train_graph_output_dir / f"{gene_name.lower()}.healthy")
        visualize_graph(Path(graph_entry["unhealthy_graph"]), train_graph_output_dir / f"{gene_name.lower()}.unhealthy")

    if evaluate_split:
        baseline_metrics = evaluate_test_set(
            METADATA_DIR / "test_manifest.json",
            train_graph_manifest_path,
            panaligner_bin,
            minigraph_bin,
            threads,
            OUTPUTS_DIR / "evaluation",
        )
        generate_feature_dataset(
            METADATA_DIR / "train_manifest.json",
            METADATA_DIR / "test_manifest.json",
            train_graph_manifest_path,
            panaligner_bin,
            threads,
            OUTPUTS_DIR / "ml",
        )
        ml_metrics = train_and_evaluate_models(
            OUTPUTS_DIR / "ml" / "features.csv",
            OUTPUTS_DIR / "ml",
            ROOT / "models",
            split_seed,
            OUTPUTS_DIR / "evaluation" / "metrics.json",
            OUTPUTS_DIR / "evaluation" / "records.json",
        )
        return {
            "baseline_accuracy": baseline_metrics["accuracy"],
            "baseline_f1_score": baseline_metrics["f1_score"],
            "primary_model": ml_metrics["primary_model"],
            "primary_model_accuracy": ml_metrics["models"][ml_metrics["primary_model"]]["accuracy"],
            "primary_model_f1_score": ml_metrics["models"][ml_metrics["primary_model"]]["f1_score"],
        }

    if query_fasta or query_sequence:
        prepared_query = prepare_query_fasta(query_sequence, query_fasta, ROOT / "data" / "queries" / "query.fa")
        if predict_ml_query:
            return predict_with_ml(
                prepared_query,
                gene,
                ml_model_path,
                train_graph_manifest_path,
                panaligner_bin,
                threads,
                OUTPUTS_DIR / "ml" / "predict",
            )
        prediction = predict_health_state(
            prepared_query,
            train_graph_manifest_path,
            panaligner_bin,
            threads,
            gene,
            OUTPUTS_DIR / "alignments",
        )
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
    parser.add_argument("--split-seed", type=int, default=42, help="Random seed used for reproducible train/test splitting.")
    parser.add_argument("--test-fraction", type=float, default=0.20, help="Fraction of each class assigned to the test set.")
    parser.add_argument("--evaluate", action="store_true", help="Run full train/test evaluation using test sequences as unseen queries.")
    parser.add_argument("--query-fasta", type=Path, default=None, help="Optional query FASTA.")
    parser.add_argument("--query-sequence", type=str, default=None, help="Optional inline query sequence.")
    parser.add_argument("--gene", type=str, default=None, help="Optional gene override.")
    parser.add_argument("--predict-ml", action="store_true", help="Use the trained ML classifier for query prediction instead of the rule-based scorer.")
    parser.add_argument("--ml-model-path", type=Path, default=ROOT / "models" / "random_forest.pkl", help="Serialized ML model used when --predict-ml is set.")
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
        args.split_seed,
        args.test_fraction,
        args.evaluate,
        args.query_fasta.resolve() if args.query_fasta else None,
        args.query_sequence,
        args.gene,
        args.predict_ml,
        args.ml_model_path.resolve(),
    )

    print("Pipeline complete.")
    if prediction:
        if args.evaluate:
            print(f"Rule-based accuracy: {prediction['baseline_accuracy']:.4f}")
            print(f"Rule-based F1-score: {prediction['baseline_f1_score']:.4f}")
            print(f"Primary ML model: {prediction['primary_model']}")
            print(f"Primary ML accuracy: {prediction['primary_model_accuracy']:.4f}")
            print(f"Primary ML F1-score: {prediction['primary_model_f1_score']:.4f}")
        else:
            print(f"Prediction: {prediction['prediction']}")
            print(f"Confidence: {prediction['confidence']:.4f}")
            if "selected_gene" in prediction:
                print(f"Selected gene: {prediction['selected_gene']}")


if __name__ == "__main__":
    main()
