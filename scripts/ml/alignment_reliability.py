"""Model A: alignment reliability (correct / incorrect / insufficient evidence).

Trains on `outputs/ml_alignment/alignment_features.csv.gz` (build_alignment_dataset.py),
whose labels come from simulation truth. Only quantities visible after running an
aligner are features: the GAF columns and tags, plus which aligner produced them.
Simulation settings, the decoy flag and the genes are never features.

Evaluation is group-disjoint: outer GroupKFold over source haplotypes. Inside each
training fold a haplotype-disjoint quarter is held back to calibrate probabilities
(isotonic, `CalibratedClassifierCV(FrozenEstimator(model))`) and to choose the two
decision thresholds:

  correct                p(incorrect) <= low, where low keeps the calibration-set error
                         among accepted alignments at or below `max_accept_error`
  incorrect              p(incorrect) >= high, where high keeps the calibration-set
                         precision of rejections at or above `min_reject_precision`
  insufficient_evidence  everything in between (an abstention band, not a trained class)

A MAPQ-only model, calibrated the same way, is the baseline every model is compared to.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Runnable as a script: put scripts/ on the import path for common, truth_eval, ml.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import OUTPUTS_DIR, ensure_dir, project_relative, write_json
from ml.build_alignment_dataset import FEATURE_COLUMNS


LABELLED = ("correct", "incorrect")
DECISIONS = ("correct", "incorrect", "insufficient_evidence")
MODEL_FEATURES = {
    "mapq_only": ["mapq"],
    "logistic": [*FEATURE_COLUMNS, "aligner_panaligner"],
    "gradient_boosting": [*FEATURE_COLUMNS, "aligner_panaligner"],
}
DEFAULT_DIR = OUTPUTS_DIR / "ml_alignment"


def model_feature_names() -> list[str]:
    return list(MODEL_FEATURES["logistic"])


def build_matrix(rows: list[dict], names: list[str]) -> np.ndarray:
    def value(row: dict, name: str) -> float:
        if name == "aligner_panaligner":
            return 1.0 if row["aligner"] == "panaligner" else 0.0
        return float(row[name])
    return np.array([[value(row, name) for name in names] for row in rows], dtype=float)


def make_estimator(model_name: str, seed: int):
    if model_name == "gradient_boosting":
        return HistGradientBoostingClassifier(random_state=seed, class_weight="balanced")
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))


def group_folds(groups: np.ndarray, n_splits: int):
    return GroupKFold(n_splits=n_splits).split(np.zeros(len(groups)), groups=groups)


def expected_calibration_error(y_true: np.ndarray, prob: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    index = np.clip(np.digitize(prob, edges[1:-1], right=True), 0, bins - 1)
    total = 0.0
    for b in range(bins):
        members = index == b
        if members.any():
            total += members.sum() / len(prob) * abs(prob[members].mean() - y_true[members].mean())
    return float(total)


def choose_thresholds(prob: np.ndarray, y: np.ndarray, *, max_accept_error: float, min_reject_precision: float) -> tuple[float, float]:
    """(low, high) on calibration data; -inf / +inf when a target cannot be met (so nothing is accepted / rejected)."""
    candidates = np.unique(prob)
    low = -math.inf
    for threshold in candidates:
        accepted = prob <= threshold
        if y[accepted].mean() <= max_accept_error:
            low = float(threshold)
    high = math.inf
    for threshold in candidates[::-1]:
        rejected = prob >= threshold
        if y[rejected].mean() >= min_reject_precision:
            high = float(threshold)
    if low >= high:  # overlapping bands: no abstention, split at the midpoint
        low = high = (low + high) / 2.0
    return low, high


def three_way(prob: np.ndarray, low: float, high: float) -> np.ndarray:
    decision = np.full(len(prob), "insufficient_evidence", dtype=object)
    decision[prob <= low] = "correct"
    decision[prob >= high] = "incorrect"
    return decision


def cross_validate(rows: list[dict], model_name: str, *, n_splits: int = 5, seed: int = 0,
                   max_accept_error: float = 0.01, min_reject_precision: float = 0.9) -> dict:
    labelled = [row for row in rows if row["label"] in LABELLED]
    names = MODEL_FEATURES[model_name]
    X = build_matrix(labelled, names)
    y = np.array([int(row["label"] == "incorrect") for row in labelled])
    groups = np.array([row["group"] for row in labelled])
    prob = np.zeros(len(labelled))
    decision = np.empty(len(labelled), dtype=object)
    thresholds = []
    for fold, (train, test) in enumerate(group_folds(groups, n_splits)):
        inner = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed + fold)
        fit_part, cal_part = next(inner.split(X[train], y[train], groups[train]))
        fit_idx, cal_idx = train[fit_part], train[cal_part]
        for part, idx in (("fit", fit_idx), ("calibration", cal_idx)):
            if len(set(y[idx])) < 2:
                raise ValueError(f"fold {fold}: the {part} part has only one class; more incorrect examples are needed")
        estimator = make_estimator(model_name, seed).fit(X[fit_idx], y[fit_idx])
        calibrated = CalibratedClassifierCV(FrozenEstimator(estimator), method="isotonic").fit(X[cal_idx], y[cal_idx])
        low, high = choose_thresholds(calibrated.predict_proba(X[cal_idx])[:, 1], y[cal_idx],
                                      max_accept_error=max_accept_error, min_reject_precision=min_reject_precision)
        prob[test] = calibrated.predict_proba(X[test])[:, 1]
        decision[test] = three_way(prob[test], low, high)
        thresholds.append({"fold": fold, "low": low, "high": high, "test_groups": sorted(set(groups[test]))})
    return {"model": model_name, "features": names, "rows": labelled, "y": y, "prob_incorrect": prob,
            "decision": decision, "thresholds": thresholds}


def _decision_metrics(y: np.ndarray, decision: np.ndarray) -> dict:
    accepted, rejected = decision == "correct", decision == "incorrect"
    return {
        "rows": int(len(y)),
        "incorrect_rows": int(y.sum()),
        "decisions": {name: int((decision == name).sum()) for name in DECISIONS},
        "coverage": float((accepted | rejected).mean()) if len(y) else 0.0,
        "accepted_error_rate": float(y[accepted].mean()) if accepted.any() else None,
        "reject_precision": float(y[rejected].mean()) if rejected.any() else None,
        "incorrect_caught": float(rejected[y == 1].mean()) if y.any() else None,
        "incorrect_accepted": int((accepted & (y == 1)).sum()),
    }


def summarize(result: dict) -> dict:
    y, prob, decision, rows = result["y"], result["prob_incorrect"], result["decision"], result["rows"]
    overall = {
        **_decision_metrics(y, decision),
        "pr_auc_incorrect": float(average_precision_score(y, prob)),
        "roc_auc": float(roc_auc_score(y, prob)),
        "brier": float(brier_score_loss(y, prob)),
        "ece": expected_calibration_error(y, prob),
        "prevalence_incorrect": float(y.mean()),
    }
    breakdown = {}
    for key in ("setting", "aligner", "decoy"):
        values = np.array([str(row[key]) for row in rows])
        breakdown[key] = {value: _decision_metrics(y[values == value], decision[values == value]) for value in sorted(set(values))}
    return {"model": result["model"], "features": result["features"], "overall": overall, "breakdown": breakdown,
            "thresholds": result["thresholds"]}


def risk_coverage(y: np.ndarray, prob: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Accept the alignments with the lowest p(incorrect) first: (coverage, error among accepted)."""
    order = np.argsort(prob, kind="stable")
    errors = np.cumsum(y[order])
    accepted = np.arange(1, len(y) + 1)
    return accepted / len(y), errors / accepted


def load_rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def plot(results: list[dict], output_dir: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []
    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    for result in results:
        coverage, risk = risk_coverage(result["y"], result["prob_incorrect"])
        axis.plot(coverage, risk, label=result["model"])
    axis.set(xlabel="Fraction of alignments accepted (lowest p(incorrect) first)", ylabel="Error rate among accepted",
             title="Model A risk-coverage (held-out haplotype folds)")
    axis.legend(fontsize=8); axis.grid(alpha=0.3); figure.tight_layout()
    path = output_dir / "risk_coverage.png"; figure.savefig(path, dpi=150); plt.close(figure); written.append(str(path))

    figure, axis = plt.subplots(figsize=(5.2, 5.2))
    edges = np.linspace(0, 1, 11)
    for result in results:
        y, prob = result["y"], result["prob_incorrect"]
        index = np.clip(np.digitize(prob, edges[1:-1], right=True), 0, 9)
        points = [(prob[index == b].mean(), y[index == b].mean()) for b in range(10) if (index == b).any()]
        axis.plot(*zip(*points), marker="o", label=result["model"])
    axis.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1)
    axis.set(xlabel="Predicted p(incorrect)", ylabel="Observed fraction incorrect", title="Model A reliability")
    axis.legend(fontsize=8); axis.grid(alpha=0.3); figure.tight_layout()
    path = output_dir / "reliability.png"; figure.savefig(path, dpi=150); plt.close(figure); written.append(str(path))
    return written


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Train and evaluate Model A (alignment reliability) with group-disjoint CV.")
    parser.add_argument("--features", type=Path, default=DEFAULT_DIR / "alignment_features.csv.gz")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DIR)
    parser.add_argument("--models", nargs="+", choices=tuple(MODEL_FEATURES), default=list(MODEL_FEATURES))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-accept-error", type=float, default=0.01)
    parser.add_argument("--min-reject-precision", type=float, default=0.9)
    args = parser.parse_args(argv)

    rows = load_rows(args.features)
    output_dir = ensure_dir(args.output_dir.resolve())
    results = [cross_validate(rows, model, n_splits=args.folds, seed=args.seed, max_accept_error=args.max_accept_error,
                              min_reject_precision=args.min_reject_precision) for model in args.models]
    summaries = [summarize(result) for result in results]
    plots = plot(results, output_dir)
    write_json(output_dir / "model_results.json", {
        "command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "features_file": project_relative(args.features),
        "label_counts": dict(Counter(row["label"] for row in rows)),
        "targets": {"max_accept_error": args.max_accept_error, "min_reject_precision": args.min_reject_precision},
        "evaluation": f"GroupKFold({args.folds}) over source haplotypes; calibration and thresholds on a haplotype-disjoint quarter of each training fold",
        "models": summaries,
        "plots": [project_relative(Path(path)) for path in plots],
    })
    print(json.dumps({summary["model"]: {key: summary["overall"][key] for key in
                      ("pr_auc_incorrect", "brier", "ece", "coverage", "accepted_error_rate", "incorrect_caught")}
                      for summary in summaries}, indent=2))


if __name__ == "__main__":
    main()
