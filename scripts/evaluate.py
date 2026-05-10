from __future__ import annotations

import argparse
import os
import platform
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, read_json, resolve_binary, write_json
from predictor import predict_health_state


def compute_binary_metrics(records: list[dict], positive_label: str = "UNHEALTHY") -> dict:
    negative_label = "HEALTHY" if positive_label == "UNHEALTHY" else "UNHEALTHY"
    tp = sum(1 for item in records if item["true_label"] == positive_label and item["predicted_label"] == positive_label)
    tn = sum(1 for item in records if item["true_label"] == negative_label and item["predicted_label"] == negative_label)
    fp = sum(1 for item in records if item["true_label"] == negative_label and item["predicted_label"] == positive_label)
    fn = sum(1 for item in records if item["true_label"] == positive_label and item["predicted_label"] == negative_label)

    total = max(1, len(records))
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    accuracy = (tp + tn) / total

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1_score": f1_score,
        "confusion_matrix": {
            "labels": [negative_label, positive_label],
            "matrix": [[tn, fp], [fn, tp]],
        },
        "counts": {
            "tp": tp,
            "tn": tn,
            "fp": fp,
            "fn": fn,
        },
    }


def draw_confusion_matrix(matrix: list[list[int]], labels: list[str], output_path: Path) -> None:
    ensure_dir(output_path.parent)
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(matrix, cmap="Blues")
    plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xticks(range(len(labels)), labels=labels)
    ax.set_yticks(range(len(labels)), labels=labels)
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    ax.set_title("Confusion Matrix")

    for row_index, row in enumerate(matrix):
        for column_index, value in enumerate(row):
            ax.text(column_index, row_index, str(value), ha="center", va="center", color="black")

    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close(fig)


def get_panaligner_version(panaligner_bin: Path) -> str:
    try:
        completed = subprocess.run(
            [
                str(panaligner_bin),
                "-cx",
                "lr",
                str((ROOT / "PanAligner" / "test" / "MT.gfa").resolve()),
                str((ROOT / "PanAligner" / "test" / "MT-orangA.fa").resolve()),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        combined = f"{completed.stdout}\n{completed.stderr}"
        match = re.search(r"PanAligner v[0-9.]+", combined)
        if match:
            return match.group(0)
    except Exception:
        pass
    return f"PanAligner binary: {panaligner_bin}"


def get_minigraph_version(minigraph_bin: Path) -> str:
    try:
        completed = subprocess.run([str(minigraph_bin), "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        output = (completed.stdout or completed.stderr).strip()
        return output or f"minigraph binary: {minigraph_bin}"
    except Exception:
        return f"minigraph binary: {minigraph_bin}"


def get_total_memory_gb() -> str:
    meminfo_path = Path("/proc/meminfo")
    if meminfo_path.exists():
        for line in meminfo_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                parts = line.split()
                if len(parts) >= 2:
                    kb = int(parts[1])
                    return f"{kb / (1024 * 1024):.2f} GB"
    return "Unavailable"


def build_experimental_setup(
    split_manifest: dict,
    metrics: dict,
    panaligner_bin: Path,
    minigraph_bin: Path,
    threads: int,
) -> str:
    lines = [
        "Experimental Setup",
        "==================",
        "1. Dataset used",
        "The evaluation used the APP, PSEN1, and PSEN2 Alzheimer’s-related gene sequence collections provided in the project workspace.",
        "",
        "2. Training/testing split",
        f"A reproducible stratified split was applied independently within each gene and label using seed {split_manifest['random_seed']}.",
        f"Train/test ratio: {(1.0 - split_manifest['test_fraction']):.0%}/{split_manifest['test_fraction']:.0%}.",
        "",
        "3. Hardware specifications",
        f"Operating system: {platform.platform()}",
        f"Processor: {platform.processor() or 'Unavailable'}",
        f"CPU cores visible to Python: {os.cpu_count()}",
        f"Total memory: {get_total_memory_gb()}",
        "",
        "4. Software tools used",
        "Python orchestration scripts, Biopython, NetworkX, Matplotlib, PanAligner, and minigraph.",
        "",
        "5. PanAligner version",
        get_panaligner_version(panaligner_bin),
        "",
        "6. Minigraph version",
        get_minigraph_version(minigraph_bin),
        "",
        "7. Parameter settings",
        "PanAligner: -cx lr",
        "minigraph: -cxggs -l1k -L1",
        f"Threads: {threads}",
        "Training phase meaning: graph construction on train sequences only.",
        "Testing phase meaning: unseen test sequences aligned and classified against train-derived graphs only.",
        "",
        "8. Evaluation metrics",
        "Accuracy, precision, recall, F1-score, confusion matrix, average alignment score, and average mapping quality.",
        "",
        "9. Number of experiments conducted",
        f"Total evaluated test sequences: {metrics['test_sequence_count']}",
    ]
    return "\n".join(lines) + "\n"


def build_result_analysis_template() -> str:
    return (
        "Result Analysis Template\n"
        "========================\n"
        "1. Quantitative analysis\n"
        "- Report accuracy, precision, recall, F1-score, and confusion matrix trends.\n\n"
        "2. Comparative analysis\n"
        "- Compare gene-wise performance and healthy vs unhealthy prediction behavior.\n\n"
        "3. Strengths\n"
        "- Highlight where graph-based alignment produced stable and interpretable classifications.\n\n"
        "4. Weaknesses\n"
        "- Note cases where healthy and unhealthy graphs were not well separated.\n\n"
        "5. Robustness\n"
        "- Discuss reproducibility under the fixed random split and consistency across genes.\n\n"
        "6. Computational efficiency\n"
        "- Summarize alignment runtime, graph size, and resource usage observations.\n\n"
        "7. Generalization capability\n"
        "- Explain how well the train-derived graph references handled unseen test sequences.\n\n"
        "8. Limitations\n"
        "- Describe class ambiguity, limited sample count for some genes, and graph collapse risks.\n\n"
        "9. Failure cases\n"
        "- Document misclassifications, tie-like predictions, and low-confidence examples.\n"
    )


def build_evaluation_report(metrics: dict, records: list[dict]) -> str:
    lines = [
        "Evaluation Report",
        "=================",
        f"Evaluated test sequences: {metrics['test_sequence_count']}",
        f"Accuracy: {metrics['accuracy']:.4f}",
        f"Precision: {metrics['precision']:.4f}",
        f"Recall: {metrics['recall']:.4f}",
        f"F1-score: {metrics['f1_score']:.4f}",
        f"Average alignment score: {metrics['average_alignment_score']:.4f}",
        f"Average mapping quality: {metrics['average_mapping_quality']:.4f}",
        "",
        "Confusion Matrix",
        str(metrics["confusion_matrix"]["matrix"]),
        "",
        "Per-gene counts",
    ]
    for gene, counts in metrics["per_gene_counts"].items():
        lines.append(f"{gene}: total={counts['total']} correct={counts['correct']}")

    lines.extend(["", "Sequence-level results"])
    for record in records:
        lines.append(
            f"{record['gene']} {record['sequence_id']}: true={record['true_label']} "
            f"predicted={record['predicted_label']} confidence={record['confidence']:.4f} "
            f"alignment_score={record['alignment_score']:.4f} mapq={record['mapping_quality']}"
        )
    return "\n".join(lines) + "\n"


def evaluate_test_set(
    test_manifest_path: Path,
    train_graph_manifest_path: Path,
    panaligner_bin: Path,
    minigraph_bin: Path,
    threads: int,
    output_dir: Path,
) -> dict:
    test_manifest = read_json(test_manifest_path)
    split_manifest = read_json(METADATA_DIR / "split_manifest.json")
    evaluation_root = ensure_dir(output_dir)
    alignments_root = ensure_dir(evaluation_root / "alignments")

    records: list[dict] = []
    per_gene_counts: dict[str, dict[str, int]] = {}

    for gene, payload in test_manifest["genes"].items():
        per_gene_counts[gene] = {"total": 0, "correct": 0}
        for label_key, label_name in (("healthy", "HEALTHY"), ("unhealthy", "UNHEALTHY")):
            for entry in payload[label_key]:
                sequence_output_dir = ensure_dir(alignments_root / gene.lower() / entry["id"])
                prediction = predict_health_state(
                    Path(entry["path"]),
                    train_graph_manifest_path,
                    panaligner_bin,
                    threads,
                    gene=gene,
                    output_dir=sequence_output_dir,
                )
                combined_alignment = prediction["combined_alignment"]
                record = {
                    "gene": gene,
                    "sequence_id": entry["id"],
                    "true_label": label_name,
                    "predicted_label": prediction["prediction"],
                    "confidence": prediction["confidence"],
                    "healthy_score": prediction["healthy_score"],
                    "unhealthy_score": prediction["unhealthy_score"],
                    "unhealthy_probability_score": prediction["unhealthy_score"] / max(
                        prediction["healthy_score"] + prediction["unhealthy_score"],
                        1e-9,
                    ),
                    "alignment_score": combined_alignment.get("alignment_score", 0.0) if combined_alignment.get("available") else 0.0,
                    "mapping_quality": combined_alignment.get("mapping_quality", 0) if combined_alignment.get("available") else 0,
                    "tie_like": prediction["tie_like"],
                }
                records.append(record)
                per_gene_counts[gene]["total"] += 1
                if record["true_label"] == record["predicted_label"]:
                    per_gene_counts[gene]["correct"] += 1

    metrics = compute_binary_metrics(records)
    labels = [1 if item["true_label"] == "UNHEALTHY" else 0 for item in records]
    scores = [item["unhealthy_probability_score"] for item in records]
    metrics["roc_auc"] = roc_auc_score(labels, scores) if len(set(labels)) > 1 else 0.0
    metrics["test_sequence_count"] = len(records)
    metrics["average_alignment_score"] = sum(item["alignment_score"] for item in records) / max(1, len(records))
    metrics["average_mapping_quality"] = sum(item["mapping_quality"] for item in records) / max(1, len(records))
    metrics["per_gene_counts"] = per_gene_counts
    metrics["tie_like_predictions"] = sum(1 for item in records if item["tie_like"])

    draw_confusion_matrix(
        metrics["confusion_matrix"]["matrix"],
        metrics["confusion_matrix"]["labels"],
        evaluation_root / "confusion_matrix.png",
    )
    write_json(evaluation_root / "metrics.json", metrics)
    (evaluation_root / "evaluation_report.txt").write_text(build_evaluation_report(metrics, records), encoding="utf-8")
    (evaluation_root / "experimental_setup.txt").write_text(
        build_experimental_setup(split_manifest, metrics, panaligner_bin, minigraph_bin, threads),
        encoding="utf-8",
    )
    (evaluation_root / "result_analysis_template.txt").write_text(build_result_analysis_template(), encoding="utf-8")
    write_json(evaluation_root / "records.json", records)

    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate train-graph / test-query performance for the AD PanAligner pipeline.")
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json", help="Path to the test manifest JSON.")
    parser.add_argument(
        "--train-graph-manifest",
        type=Path,
        default=METADATA_DIR / "train_graph_manifest.json",
        help="Path to the train graph manifest JSON.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument(
        "--panaligner-bin",
        type=str,
        default=None,
        help="Path to the PanAligner binary. Defaults to ./PanAligner/PanAligner if present.",
    )
    parser.add_argument(
        "--minigraph-bin",
        type=str,
        default=None,
        help="Path to the minigraph binary. Defaults to ./minigraph/minigraph if present.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUTS_DIR / "evaluation",
        help="Directory where evaluation metrics and reports will be written.",
    )
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph")])
    metrics = evaluate_test_set(
        args.test_manifest.resolve(),
        args.train_graph_manifest.resolve(),
        panaligner_bin,
        minigraph_bin,
        args.threads,
        args.output_dir.resolve(),
    )
    print(f"Evaluation complete. Accuracy={metrics['accuracy']:.4f} F1={metrics['f1_score']:.4f}")
    print(f"Metrics written to {(args.output_dir.resolve() / 'metrics.json')}")


if __name__ == "__main__":
    main()
