from __future__ import annotations

import argparse
import csv
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from common import OUTPUTS_DIR, ROOT, ensure_dir, read_json, write_json


MODEL_INPUT_COLUMNS = [
    "gene",
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


def load_feature_rows(features_csv: Path) -> list[dict]:
    with features_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def prepare_model_input(row: dict) -> dict:
    prepared = {"gene": row["gene"]}
    for column in MODEL_INPUT_COLUMNS:
        if column == "gene":
            continue
        prepared[column] = float(row[column])
    return prepared


def split_feature_rows(rows: list[dict]) -> tuple[list[dict], list[int], list[dict], list[int]]:
    train_rows = [row for row in rows if row["split"] == "train"]
    test_rows = [row for row in rows if row["split"] == "test"]
    x_train = [prepare_model_input(row) for row in train_rows]
    y_train = [int(row["label"]) for row in train_rows]
    x_test = [prepare_model_input(row) for row in test_rows]
    y_test = [int(row["label"]) for row in test_rows]
    return x_train, y_train, x_test, y_test


def build_model_pipelines(random_seed: int) -> dict[str, Pipeline]:
    return {
        "random_forest": Pipeline(
            [
                ("vectorizer", DictVectorizer(sparse=False)),
                ("classifier", RandomForestClassifier(n_estimators=300, random_state=random_seed, class_weight="balanced")),
            ]
        ),
        "svm": Pipeline(
            [
                ("vectorizer", DictVectorizer(sparse=False)),
                ("scaler", StandardScaler()),
                ("classifier", SVC(probability=True, kernel="rbf", class_weight="balanced", random_state=random_seed)),
            ]
        ),
        "logistic_regression": Pipeline(
            [
                ("vectorizer", DictVectorizer(sparse=False)),
                ("scaler", StandardScaler()),
                ("classifier", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=random_seed)),
            ]
        ),
    }


def evaluate_model(pipeline: Pipeline, x_test: list[dict], y_test: list[int]) -> tuple[dict, list[int], list[float]]:
    predictions = pipeline.predict(x_test)
    probabilities = pipeline.predict_proba(x_test)[:, 1]
    metrics = {
        "accuracy": accuracy_score(y_test, predictions),
        "precision": precision_score(y_test, predictions, zero_division=0),
        "recall": recall_score(y_test, predictions, zero_division=0),
        "f1_score": f1_score(y_test, predictions, zero_division=0),
        "roc_auc": roc_auc_score(y_test, probabilities) if len(set(y_test)) > 1 else 0.0,
        "confusion_matrix": confusion_matrix(y_test, predictions).tolist(),
        "classification_report": classification_report(y_test, predictions, target_names=["healthy", "unhealthy"], zero_division=0),
    }
    return metrics, predictions.tolist(), probabilities.tolist()


def compute_safe_roc_curve(y_true: list[int], probabilities: list[float]) -> tuple[list[float], list[float]] | None:
    if len(set(y_true)) < 2:
        return None
    fpr, tpr, _ = roc_curve(y_true, probabilities)
    return fpr.tolist(), tpr.tolist()


def draw_confusion_matrix(matrix: list[list[int]], output_path: Path, title: str) -> None:
    ensure_dir(output_path.parent)
    plt.figure(figsize=(6, 5))
    sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", xticklabels=["healthy", "unhealthy"], yticklabels=["healthy", "unhealthy"])
    plt.xlabel("Predicted label")
    plt.ylabel("True label")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def draw_roc_curves(curves: dict[str, tuple[list[float], list[float], float]], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    plt.figure(figsize=(7, 6))
    for model_name, (fpr, tpr, auc_value) in curves.items():
        plt.plot(fpr, tpr, label=f"{model_name} (AUC={auc_value:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="chance")
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("ROC Curve Comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def draw_class_distribution(rows: list[dict], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    counts = {}
    for split_name in ("train", "test"):
        for label_name in ("HEALTHY", "UNHEALTHY"):
            counts[(split_name, label_name)] = sum(1 for row in rows if row["split"] == split_name and row["label_name"] == label_name)

    categories = ["train healthy", "train unhealthy", "test healthy", "test unhealthy"]
    values = [
        counts[("train", "HEALTHY")],
        counts[("train", "UNHEALTHY")],
        counts[("test", "HEALTHY")],
        counts[("test", "UNHEALTHY")],
    ]
    plt.figure(figsize=(8, 5))
    sns.barplot(x=categories, y=values, hue=categories, palette="Set2", legend=False)
    plt.ylabel("Sequence count")
    plt.title("Class Distribution Across Train/Test Splits")
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def draw_prediction_probability_plot(probabilities: list[float], y_test: list[int], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    healthy_scores = [score for score, label in zip(probabilities, y_test) if label == 0]
    unhealthy_scores = [score for score, label in zip(probabilities, y_test) if label == 1]
    plt.figure(figsize=(8, 5))
    sns.histplot(healthy_scores, color="steelblue", label="healthy", kde=True, stat="density", bins=12, alpha=0.45)
    sns.histplot(unhealthy_scores, color="tomato", label="unhealthy", kde=True, stat="density", bins=12, alpha=0.45)
    plt.xlabel("Predicted probability of unhealthy")
    plt.title("Prediction Probability Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def draw_feature_importance_plot(feature_names: list[str], importances: list[float], output_path: Path, top_n: int = 12) -> list[tuple[str, float]]:
    ensure_dir(output_path.parent)
    ranked = sorted(zip(feature_names, importances), key=lambda item: item[1], reverse=True)[:top_n]
    labels = [item[0] for item in ranked][::-1]
    values = [item[1] for item in ranked][::-1]
    plt.figure(figsize=(9, 6))
    sns.barplot(x=values, y=labels, hue=labels, palette="viridis", legend=False)
    plt.xlabel("Importance")
    plt.title("Random Forest Feature Importance")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()
    return ranked


def infer_feature_importance_text(top_features: list[tuple[str, float]]) -> str:
    interpretation = {
        "healthy_overlap_score": "Strong healthy-graph alignment evidence tends to support the healthy class.",
        "unhealthy_overlap_score": "Strong unhealthy-graph alignment evidence tends to support the disease-associated class.",
        "alignment_score": "Overall graph-alignment quality reflects how well the query fits the training-derived pangenome.",
        "graph_path_complexity": "More complex traversal patterns can reflect variant-rich or branch-heavy alignments.",
        "identity_percentage": "Higher identity indicates stronger base-level agreement with graph-supported paths.",
        "mismatches": "Mismatch burden highlights divergence from graph segments traversed by the alignment.",
    }
    lines = ["Random Forest Feature Importance", "==============================="]
    for name, value in top_features:
        lines.append(f"{name}: {value:.6f}")
        if name in interpretation:
            lines.append(f"  Interpretation: {interpretation[name]}")
    return "\n".join(lines) + "\n"


def build_ml_analysis_report(model_metrics: dict, baseline_metrics: dict | None, top_features: list[tuple[str, float]]) -> str:
    lines = [
        "ML Analysis Report",
        "==================",
        "PanAligner was used as a graph-based biological feature extraction engine, while machine learning models were trained on graph-alignment-derived features for downstream Alzheimer’s Disease prediction.",
        "",
        "Model comparison",
    ]
    if baseline_metrics:
        lines.append(
            f"Rule-based baseline: accuracy={baseline_metrics['accuracy']:.4f}, precision={baseline_metrics['precision']:.4f}, "
            f"recall={baseline_metrics['recall']:.4f}, f1={baseline_metrics['f1_score']:.4f}, roc_auc={baseline_metrics.get('roc_auc', 0.0):.4f}"
        )
    for model_name, metrics in model_metrics.items():
        lines.append(
            f"{model_name}: accuracy={metrics['accuracy']:.4f}, precision={metrics['precision']:.4f}, "
            f"recall={metrics['recall']:.4f}, f1={metrics['f1_score']:.4f}, roc_auc={metrics['roc_auc']:.4f}"
        )

    lines.extend(
        [
            "",
            "Feature importance analysis",
            "The random forest rankings identify which PanAligner-derived alignment signals contribute most to healthy vs unhealthy separation.",
        ]
    )
    for feature_name, importance in top_features[:5]:
        lines.append(f"- {feature_name}: {importance:.6f}")

    lines.extend(
        [
            "",
            "Strengths",
            "- The models operate on biologically grounded graph-alignment features rather than raw sequence embeddings.",
            "- Random forest, SVM, and logistic regression remain suitable for small datasets and are easier to explain during project review.",
            "",
            "Weaknesses",
            "- Training samples are also used in graph construction, so the feature distribution can remain optimistic for in-graph samples.",
            "- Small class counts, especially for PSEN2, limit statistical stability.",
            "",
            "Overfitting discussion",
            "- With only ~220 labeled sequences overall, complex models are intentionally avoided. Random forest depth and linear margins should still be interpreted carefully.",
            "",
            "Dataset limitations",
            "- If healthy and unhealthy graph alignments are nearly identical for a gene, both rule-based and ML stages will struggle to separate classes.",
            "",
            "Biological interpretation",
            "- The ML stage does not learn directly from DNA strings. It learns from PanAligner-derived alignment quality, path traversal, and graph-overlap signals.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_prediction_records(
    rows: list[dict],
    y_test: list[int],
    results: dict[str, dict],
    output_path: Path,
) -> None:
    test_rows = [row for row in rows if row["split"] == "test"]
    fieldnames = [
        "sequence_id",
        "gene",
        "true_label",
        "random_forest_prediction",
        "random_forest_probability_unhealthy",
        "svm_prediction",
        "svm_probability_unhealthy",
        "logistic_regression_prediction",
        "logistic_regression_probability_unhealthy",
    ]
    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, row in enumerate(test_rows):
            writer.writerow(
                {
                    "sequence_id": row["sequence_id"],
                    "gene": row["gene"],
                    "true_label": "UNHEALTHY" if y_test[index] == 1 else "HEALTHY",
                    "random_forest_prediction": "UNHEALTHY" if results["random_forest"]["predictions"][index] == 1 else "HEALTHY",
                    "random_forest_probability_unhealthy": results["random_forest"]["probabilities"][index],
                    "svm_prediction": "UNHEALTHY" if results["svm"]["predictions"][index] == 1 else "HEALTHY",
                    "svm_probability_unhealthy": results["svm"]["probabilities"][index],
                    "logistic_regression_prediction": "UNHEALTHY" if results["logistic_regression"]["predictions"][index] == 1 else "HEALTHY",
                    "logistic_regression_probability_unhealthy": results["logistic_regression"]["probabilities"][index],
                }
            )


def train_and_evaluate_models(
    features_csv: Path,
    output_dir: Path,
    models_dir: Path,
    random_seed: int,
    baseline_metrics_path: Path | None = None,
    baseline_records_path: Path | None = None,
) -> dict:
    rows = load_feature_rows(features_csv)
    x_train, y_train, x_test, y_test = split_feature_rows(rows)
    model_pipelines = build_model_pipelines(random_seed)
    ml_output_dir = ensure_dir(output_dir)
    models_root = ensure_dir(models_dir)

    results: dict[str, dict] = {}
    roc_curves: dict[str, tuple[list[float], list[float], float]] = {}
    classification_report_sections: list[str] = []
    comparison_rows: list[dict] = []

    for model_name, pipeline in model_pipelines.items():
        pipeline.fit(x_train, y_train)
        metrics, predictions, probabilities = evaluate_model(pipeline, x_test, y_test)
        results[model_name] = metrics | {"predictions": predictions, "probabilities": probabilities}
        joblib.dump(pipeline, models_root / f"{model_name}.pkl")
        curve_points = compute_safe_roc_curve(y_test, probabilities)
        if curve_points is not None:
            roc_curves[model_name] = (curve_points[0], curve_points[1], metrics["roc_auc"])
        classification_report_sections.append(f"[{model_name}]\n{metrics['classification_report']}\n")
        comparison_rows.append(
            {
                "model": model_name,
                "accuracy": metrics["accuracy"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1_score": metrics["f1_score"],
                "roc_auc": metrics["roc_auc"],
            }
        )

    baseline_metrics = read_json(baseline_metrics_path) if baseline_metrics_path and baseline_metrics_path.exists() else None
    baseline_records = read_json(baseline_records_path) if baseline_records_path and baseline_records_path.exists() else None
    if baseline_metrics:
        comparison_rows.insert(
            0,
            {
                "model": "rule_based_panaligner",
                "accuracy": baseline_metrics["accuracy"],
                "precision": baseline_metrics["precision"],
                "recall": baseline_metrics["recall"],
                "f1_score": baseline_metrics["f1_score"],
                "roc_auc": baseline_metrics.get("roc_auc", 0.0),
            },
        )
        if baseline_records:
            baseline_labels = [1 if item["true_label"] == "UNHEALTHY" else 0 for item in baseline_records]
            baseline_scores = [item.get("unhealthy_probability_score", 0.0) for item in baseline_records]
            curve_points = compute_safe_roc_curve(baseline_labels, baseline_scores)
            if curve_points is not None:
                roc_curves["rule_based_panaligner"] = (curve_points[0], curve_points[1], baseline_metrics.get("roc_auc", 0.0))

    primary_model = "random_forest"
    primary_metrics = results[primary_model]
    draw_confusion_matrix(primary_metrics["confusion_matrix"], ml_output_dir / "confusion_matrix.png", "Random Forest Confusion Matrix")
    draw_roc_curves({name: (values[0], values[1], values[2]) for name, values in roc_curves.items()}, ml_output_dir / "roc_curve.png")
    draw_class_distribution(rows, ml_output_dir / "class_distribution.png")
    draw_prediction_probability_plot(primary_metrics["probabilities"], y_test, ml_output_dir / "prediction_probability.png")

    rf_pipeline = model_pipelines["random_forest"]
    vectorizer = rf_pipeline.named_steps["vectorizer"]
    classifier = rf_pipeline.named_steps["classifier"]
    feature_names = vectorizer.get_feature_names_out().tolist()
    importances = classifier.feature_importances_.tolist()
    top_features = draw_feature_importance_plot(feature_names, importances, ml_output_dir / "feature_importance.png")
    (ml_output_dir / "feature_importance.txt").write_text(infer_feature_importance_text(top_features), encoding="utf-8")

    metrics_payload = {
        "primary_model": primary_model,
        "models": {name: {key: value for key, value in metrics.items() if key not in {"predictions", "probabilities"}} for name, metrics in results.items()},
        "baseline": baseline_metrics,
    }
    write_json(ml_output_dir / "metrics.json", metrics_payload)
    (ml_output_dir / "classification_report.txt").write_text("\n".join(classification_report_sections), encoding="utf-8")
    (ml_output_dir / "ml_analysis_report.txt").write_text(build_ml_analysis_report(metrics_payload["models"], baseline_metrics, top_features), encoding="utf-8")
    write_prediction_records(rows, y_test, results, ml_output_dir / "prediction_records.csv")

    comparison_path = ml_output_dir / "comparison_table.csv"
    with comparison_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["model", "accuracy", "precision", "recall", "f1_score", "roc_auc"])
        writer.writeheader()
        writer.writerows(comparison_rows)

    return metrics_payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Train lightweight ML models on PanAligner-derived graph alignment features.")
    parser.add_argument("--features-csv", type=Path, default=OUTPUTS_DIR / "ml" / "features.csv", help="ML feature CSV produced by extract_features.py.")
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for reproducible model training.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "ml", help="Directory where ML reports and plots will be written.")
    parser.add_argument("--models-dir", type=Path, default=ROOT / "models", help="Directory where trained models will be serialized.")
    parser.add_argument("--baseline-metrics", type=Path, default=OUTPUTS_DIR / "evaluation" / "metrics.json", help="Optional rule-based baseline metrics JSON.")
    parser.add_argument("--baseline-records", type=Path, default=OUTPUTS_DIR / "evaluation" / "records.json", help="Optional rule-based baseline records JSON.")
    args = parser.parse_args()

    metrics = train_and_evaluate_models(
        args.features_csv.resolve(),
        args.output_dir.resolve(),
        args.models_dir.resolve(),
        args.random_seed,
        args.baseline_metrics.resolve() if args.baseline_metrics else None,
        args.baseline_records.resolve() if args.baseline_records else None,
    )
    print(f"ML training complete. Primary model={metrics['primary_model']}")
    print(f"Metrics written to {(args.output_dir.resolve() / 'metrics.json')}")


if __name__ == "__main__":
    main()
