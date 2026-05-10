from __future__ import annotations

import argparse
from pathlib import Path

import joblib

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, read_json, resolve_binary, write_json
from extract_features import extract_sequence_features
from predictor import infer_gene
from train_ml import prepare_model_input


def rank_feature_contributions(model_pipeline, feature_row: dict) -> list[tuple[str, float]]:
    vectorizer = model_pipeline.named_steps["vectorizer"]
    transformed = vectorizer.transform([prepare_model_input(feature_row)])[0]
    feature_names = vectorizer.get_feature_names_out().tolist()
    classifier = model_pipeline.named_steps["classifier"]

    if hasattr(classifier, "feature_importances_"):
        scores = [abs(value) * importance for value, importance in zip(transformed, classifier.feature_importances_)]
    elif hasattr(classifier, "coef_"):
        coefficients = classifier.coef_[0]
        scores = [abs(value * coefficient) for value, coefficient in zip(transformed, coefficients)]
    else:
        scores = [abs(value) for value in transformed]

    ranked = sorted(zip(feature_names, scores), key=lambda item: item[1], reverse=True)
    return ranked[:5]


def predict_with_ml(
    query_fasta: Path,
    gene: str | None,
    model_path: Path,
    train_graph_manifest: Path,
    panaligner_bin: Path,
    threads: int,
    output_dir: Path,
) -> dict:
    model_pipeline = joblib.load(model_path)
    graph_manifest = read_json(train_graph_manifest)
    resolved_output_dir = ensure_dir(output_dir)
    selected_gene = gene.upper() if gene else infer_gene(query_fasta, graph_manifest, panaligner_bin, threads, resolved_output_dir)[0]
    feature_row = extract_sequence_features(
        {
            "id": query_fasta.stem,
            "path": str(query_fasta.resolve()),
            "label": "UNKNOWN",
        },
        "query",
        selected_gene,
        graph_manifest,
        panaligner_bin,
        threads,
        resolved_output_dir,
        {},
    )
    probabilities = model_pipeline.predict_proba([prepare_model_input(feature_row)])[0]
    predicted_label = int(model_pipeline.predict([prepare_model_input(feature_row)])[0])
    confidence = float(probabilities[predicted_label])
    ranked_features = rank_feature_contributions(model_pipeline, feature_row)
    result = {
        "query_fasta": str(query_fasta.resolve()),
        "selected_gene": selected_gene,
        "prediction": "UNHEALTHY" if predicted_label == 1 else "HEALTHY",
        "confidence": confidence,
        "top_features": ranked_features,
        "feature_row": feature_row,
    }
    write_json(resolved_output_dir / "prediction.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict healthy vs unhealthy status using the trained ML classifier on PanAligner-derived features.")
    parser.add_argument("--query-fasta", type=Path, required=True, help="Query FASTA file to classify.")
    parser.add_argument("--gene", type=str, default=None, help="Optional gene name for the query (APP, PSEN1, PSEN2). If omitted, the best-matching gene is inferred from combined train graphs.")
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument("--model-path", type=Path, default=ROOT / "models" / "random_forest.pkl", help="Serialized model to load.")
    parser.add_argument("--train-graph-manifest", type=Path, default=METADATA_DIR / "train_graph_manifest.json", help="Train graph manifest JSON.")
    parser.add_argument("--panaligner-bin", type=str, default=None, help="Path to PanAligner binary.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "ml" / "predict", help="Directory where query alignment artifacts will be written.")
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    result = predict_with_ml(
        args.query_fasta.resolve(),
        args.gene,
        args.model_path.resolve(),
        args.train_graph_manifest.resolve(),
        panaligner_bin,
        args.threads,
        args.output_dir.resolve(),
    )
    print(f"Prediction: {result['prediction']}")
    print(f"Selected gene: {result['selected_gene']}")
    print(f"Confidence: {result['confidence'] * 100:.1f}%")
    print("Top contributing features:")
    for feature_name, score in result["top_features"][:3]:
        print(f"- {feature_name}: {score:.6f}")


if __name__ == "__main__":
    main()
