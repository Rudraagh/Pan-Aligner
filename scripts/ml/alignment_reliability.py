"""Model A: alignment reliability (correct / incorrect / insufficient evidence).

Trains on `outputs/ml_alignment/alignment_features.csv.gz` (build_alignment_dataset.py),
whose labels come from simulation truth. Only quantities visible after running an
aligner are features: the GAF columns and tags, plus which aligner produced them.
Simulation settings, the decoy flag and the genes are never features.

Evaluation is group-disjoint, with three groupings:

  window              (default) GroupKFold over 20 kb windows of the read's true position,
                      with reads from adjacent windows removed from training
  leave_one_gene_out  test on one gene's reads; train without that gene's reads or graph
  haplotype           GroupKFold over source haplotypes. Weak: every haplotype is the same
                      reference plus one variant, so the same locus is on both sides. Kept
                      only to compare with the first results in outputs/ml_alignment/.

Inside each training fold a group-disjoint quarter is held back to calibrate probabilities
(isotonic, `CalibratedClassifierCV(FrozenEstimator(model))`) and to choose the two
decision thresholds:

  correct                p(incorrect) <= low, where low keeps the calibration-set error
                         among accepted alignments at or below `max_accept_error`
  incorrect              p(incorrect) >= high, where high keeps the calibration-set
                         precision of rejections at or above `min_reject_precision`
  insufficient_evidence  everything in between (an abstention band, not a trained class)

Baselines: a MAPQ-only model calibrated the same way, and a plain MAPQ threshold sweep
(accept when MAPQ >= Q), which involves no training at all.
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
WINDOW_BP = 20000
GROUPINGS = {
    "window": f"GroupKFold over {WINDOW_BP // 1000} kb windows of the read's true position, with adjacent windows removed from training",
    "leave_one_gene_out": "each gene's reads tested in turn, with that gene's reads and graph removed from training",
    "haplotype": "GroupKFold over source haplotypes (weak: every haplotype is the same reference plus one variant)",
}
MAPQ_THRESHOLDS = (0, 1, 5, 10, 20, 30, 40, 50, 60)


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


def window_index(row: dict) -> int:
    """20 kb window of the read's true start on its own gene's reference."""
    return int(row["read"].split("!")[2]) // WINDOW_BP


def evaluation_folds(rows: list[dict], grouping: str, n_splits: int) -> tuple[list[tuple[np.ndarray, np.ndarray]], np.ndarray]:
    """(train, test) index folds, and the groups used to keep calibration data disjoint inside a fold."""
    if grouping not in GROUPINGS:
        raise ValueError(f"grouping must be one of {sorted(GROUPINGS)}")
    if grouping == "haplotype":
        groups = np.array([row["group"] for row in rows])
        return list(group_folds(groups, n_splits)), groups
    own_genes = np.array([row["own_gene"] for row in rows])
    indices = np.array([window_index(row) for row in rows])
    groups = np.array([f"{gene}:{index}" for gene, index in zip(own_genes, indices)])
    if grouping == "leave_one_gene_out":
        graph_genes = np.array([row["graph_gene"] for row in rows])
        return [(np.flatnonzero((own_genes != gene) & (graph_genes != gene)), np.flatnonzero(own_genes == gene))
                for gene in sorted(set(own_genes))], groups
    folds = []
    for train, test in group_folds(groups, n_splits):
        neighbours = {f"{own_genes[i]}:{indices[i] + offset}" for i in test for offset in (-1, 1)}
        folds.append((train[~np.isin(groups[train], sorted(neighbours))], test))
    return folds, groups


def output_dir_for(output_dir: Path | None, grouping: str) -> Path:
    """Results go to their own folder, so the first results in outputs/ml_alignment/ stay as committed."""
    return output_dir if output_dir is not None else DEFAULT_DIR / "regrouped" / grouping


def mapq_threshold_sweep(y: np.ndarray, mapq: np.ndarray, thresholds=MAPQ_THRESHOLDS) -> list[dict]:
    """Accept alignments with MAPQ >= Q: no training, no calibration."""
    rows = []
    for threshold in thresholds:
        accepted = mapq >= threshold
        rows.append({
            "threshold": int(threshold),
            "coverage": float(accepted.mean()),
            "accepted_error_rate": float(y[accepted].mean()) if accepted.any() else None,
            "incorrect_caught": float((~accepted)[y == 1].mean()) if y.any() else None,
        })
    return rows


def model_error_at_coverage(y: np.ndarray, prob: np.ndarray, coverage: float) -> float | None:
    """Error among the `coverage` fraction of alignments with the lowest p(incorrect)."""
    count = round(coverage * len(y))
    if count == 0:
        return None
    return float(y[np.argsort(prob, kind="stable")[:count]].mean())


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


def fold_parts(train: np.ndarray, y: np.ndarray, groups: np.ndarray, *, seed: int, separate_threshold_set: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(fit, calibration, thresholds) indices inside a training fold, disjoint by group.

    By default a quarter is held back and used for both calibration and thresholds. With
    `separate_threshold_set`, 60% fits, 20% calibrates and 20% only sets thresholds, so
    thresholds never see the rows the calibrator was fitted on.
    """
    if not separate_threshold_set:
        inner = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed)
        fit_part, cal_part = next(inner.split(train, y[train], groups[train]))
        return train[fit_part], train[cal_part], train[cal_part]
    inner = GroupShuffleSplit(n_splits=1, test_size=0.4, random_state=seed)
    fit_part, held = next(inner.split(train, y[train], groups[train]))
    held_idx = train[held]
    second = GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=seed + 1000)
    cal_part, thr_part = next(second.split(held_idx, y[held_idx], groups[held_idx]))
    return train[fit_part], held_idx[cal_part], held_idx[thr_part]


def coverage_threshold(prob: np.ndarray, coverage: float) -> float:
    """Highest p(incorrect) to accept so that `coverage` of these rows are accepted."""
    return float(np.quantile(prob, coverage, method="inverted_cdf"))


def apply_targets(result: dict, *, max_accept_error: float, min_reject_precision: float, target_coverage: float | None = None) -> tuple[np.ndarray, list[dict]]:
    """Decisions for other targets without refitting; thresholds come only from each fold's validation part."""
    decision = np.empty(len(result["y"]), dtype=object)
    thresholds = []
    for fold in result["validation"]:
        low, high = choose_thresholds(fold["prob"], fold["y"], max_accept_error=max_accept_error, min_reject_precision=min_reject_precision)
        if target_coverage is not None:
            low = coverage_threshold(fold["prob"], target_coverage)
        decision[fold["test"]] = three_way(result["prob_incorrect"][fold["test"]], low, high)
        thresholds.append({"fold": fold["fold"], "low": low, "high": high})
    return decision, thresholds


def cross_validate(rows: list[dict], model_name: str, *, n_splits: int = 5, seed: int = 0,
                   max_accept_error: float = 0.01, min_reject_precision: float = 0.9, grouping: str = "haplotype",
                   separate_threshold_set: bool = False) -> dict:
    labelled = [row for row in rows if row["label"] in LABELLED]
    names = MODEL_FEATURES[model_name]
    X = build_matrix(labelled, names)
    y = np.array([int(row["label"] == "incorrect") for row in labelled])
    folds, groups = evaluation_folds(labelled, grouping, n_splits)
    prob = np.zeros(len(labelled))
    decision = np.empty(len(labelled), dtype=object)
    thresholds = []
    validation = []
    for fold, (train, test) in enumerate(folds):
        fit_idx, cal_idx, thr_idx = fold_parts(train, y, groups, seed=seed + fold, separate_threshold_set=separate_threshold_set)
        for part, idx in (("fit", fit_idx), ("calibration", cal_idx), ("threshold", thr_idx)):
            if len(set(y[idx])) < 2:
                raise ValueError(f"fold {fold}: the {part} part has only one class; more incorrect examples are needed")
        estimator = make_estimator(model_name, seed).fit(X[fit_idx], y[fit_idx])
        calibrated = CalibratedClassifierCV(FrozenEstimator(estimator), method="isotonic").fit(X[cal_idx], y[cal_idx])
        validation_prob = calibrated.predict_proba(X[thr_idx])[:, 1]
        low, high = choose_thresholds(validation_prob, y[thr_idx],
                                      max_accept_error=max_accept_error, min_reject_precision=min_reject_precision)
        prob[test] = calibrated.predict_proba(X[test])[:, 1]
        decision[test] = three_way(prob[test], low, high)
        validation.append({"fold": fold, "prob": validation_prob, "y": y[thr_idx], "test": test})
        thresholds.append({"fold": fold, "low": low, "high": high, "test_groups": sorted(set(groups[test])),
                           "calibration_groups": sorted(set(groups[cal_idx])), "threshold_groups": sorted(set(groups[thr_idx]))})
    return {"model": model_name, "grouping": grouping, "features": names, "rows": labelled, "y": y, "groups": groups, "prob_incorrect": prob,
            "decision": decision, "thresholds": thresholds, "validation": validation}


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
    # Incorrect alignments split into two different problems; correct rows are shared by both.
    decoy = np.array([int(row["decoy"]) for row in rows])
    for name, keep in (("within_locus", (decoy == 0)), ("off_target", (y == 0) | (decoy == 1))):
        subset_y, subset_prob = y[keep], prob[keep]
        overall[name] = {
            "rows": int(keep.sum()), "incorrect_rows": int(subset_y.sum()),
            "prevalence_incorrect": float(subset_y.mean()) if keep.any() else None,
            "pr_auc_incorrect": float(average_precision_score(subset_y, subset_prob)) if 0 < subset_y.sum() < len(subset_y) else None,
            "roc_auc": float(roc_auc_score(subset_y, subset_prob)) if 0 < subset_y.sum() < len(subset_y) else None,
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


def plot(results: list[dict], output_dir: Path, grouping: str = "haplotype") -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []
    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    for result in results:
        coverage, risk = risk_coverage(result["y"], result["prob_incorrect"])
        axis.plot(coverage, risk, label=result["model"])
    axis.set(xlabel="Fraction of alignments accepted (lowest p(incorrect) first)", ylabel="Error rate among accepted",
             title=f"Model A risk-coverage ({grouping} folds)")
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


def write_report(path: Path, rows: list[dict], summaries: list[dict], targets: dict, grouping: str = "haplotype", mapq: dict | None = None) -> None:
    """Markdown report whose every number comes from the rows and summaries passed in."""
    counts = Counter(row["label"] for row in rows)
    incorrect = [row for row in rows if row["label"] == "incorrect"]
    decoy_incorrect = sum(int(row["decoy"]) == 1 for row in incorrect)

    def number(value, digits=3):
        return "–" if value is None else f"{value:.{digits}f}"

    table = []
    split_table = []
    notes = []
    for summary in summaries:
        o = summary["overall"]
        table.append(f"| {summary['model']} | {number(o['pr_auc_incorrect'])} | {number(o['roc_auc'])} | {number(o['brier'], 4)} | "
                     f"{number(o['ece'], 4)} | {number(o['coverage'])} | {number(o['accepted_error_rate'], 4)} | {number(o['incorrect_caught'])} |")
        split_table.append(f"| {summary['model']} | {number(o['within_locus']['pr_auc_incorrect'])} | {number(o['within_locus']['roc_auc'])} | "
                           f"{number(o['off_target']['pr_auc_incorrect'])} | {number(o['off_target']['roc_auc'])} |")
        finite = [t for t in summary["thresholds"] if math.isfinite(t["low"]) and math.isfinite(t["high"])]
        collapsed = sum(t["low"] == t["high"] for t in finite)
        if collapsed:
            notes.append(f"- `{summary['model']}`: abstention band collapsed in {collapsed} of {len(summary['thresholds'])} folds "
                         "(the accept and reject targets overlapped), so `insufficient_evidence` was not used there.")
        own = summary["breakdown"].get("decoy", {}).get("0")
        if own and own["incorrect_rows"]:
            notes.append(f"- `{summary['model']}` on within-locus alignments: {own['incorrect_rows']} incorrect of {own['rows']}; "
                         f"caught {number(own['incorrect_caught'])}, {own['incorrect_accepted']} accepted as correct, "
                         f"{own['decisions']['incorrect']} flagged incorrect in total.")
    lines = [
        "# Model A: alignment reliability", "",
        "Labels come from simulation truth (PBSIM2 reads from the held-out haplotypes), never from the aligners being judged.", "",
        f"Rows: {counts.get('correct', 0)} correct, {counts.get('incorrect', 0)} incorrect, {counts.get('unaligned', 0)} unaligned "
        f"(not used), {counts.get('unprojectable', 0)} unprojectable (not used).", "",
        f"{decoy_incorrect} of {len(incorrect)} incorrect alignments are off-target: a read aligned to another gene's graph, "
        "where its true locus was not available to the aligner. The rest are within-locus misplacements.", "",
        f"Evaluation: {GROUPINGS[grouping]}; isotonic calibration and thresholds on a group-disjoint part of each "
        f"training fold. Targets: error among accepted <= {targets['max_accept_error']}, rejection precision >= {targets['min_reject_precision']}.", "",
        "| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Coverage | Error among accepted | Incorrect caught |",
        "|---|---|---|---|---|---|---|---|", *table, "",
        "## Within-locus and off-target errors separately", "",
        "Both columns share the same correct rows, but within-locus incorrect rows are far rarer, so compare each column with its own prevalence.", "",
        "| Model | Within-locus PR-AUC | Within-locus ROC-AUC | Off-target PR-AUC | Off-target ROC-AUC |", "|---|---|---|---|---|", *split_table, "",
        *_mapq_section(mapq, number),
        "## Notes", "", *(notes or ["- none"]), "",
        "## Limits", "",
        *[f"- {limit}" for limit in report_limits(rows)],
        "- Simulated reads (PBSIM2, locus scale); real sequencing data and whole-genome graphs are not covered.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def report_limits(rows: list[dict]) -> list[str]:
    """Limits that depend on what the incorrect class actually contains."""
    incorrect = [row for row in rows if row["label"] == "incorrect"]
    if not incorrect:
        return []
    off_target = sum(int(row["decoy"]) == 1 for row in incorrect) / len(incorrect)
    if off_target >= 0.9:
        limits = ["The incorrect class is dominated by off-target reads on single-locus graphs, so the models mainly learn that low identity "
                  "and coverage mean wrong. This is evidence about detecting off-target or contaminating reads, not about whole-genome alignment."]
    else:
        limits = [f"{off_target:.0%} of incorrect alignments are off-target reads on single-locus graphs, where the true locus was not available "
                  "to the aligner. This says nothing direct about whole-genome alignment."]
    within = [row for row in incorrect if int(row["decoy"]) == 0]
    if within:
        aligner, count = Counter(row["aligner"] for row in within).most_common(1)[0]
        limits.append(f"Of {len(within)} within-locus errors, {count / len(within):.0%} from {aligner}, so those results mostly describe that aligner. "
                      "Errors from the same regions repeat, so the effective sample is smaller than the count.")
    return limits


def _mapq_section(mapq: dict | None, number) -> list[str]:
    if not mapq:
        return []
    lines = ["## Plain MAPQ thresholds (no training)", "",
             f"Accept when MAPQ >= Q. `{mapq['compared_model']} at same coverage` is that model's error among the same share of alignments, chosen by lowest p(incorrect).", "",
             f"| Q | Coverage | Error among accepted | Incorrect caught | {mapq['compared_model']} at same coverage |", "|---|---|---|---|---|"]
    lines += [f"| {row['threshold']} | {number(row['coverage'])} | {number(row['accepted_error_rate'], 4)} | {number(row['incorrect_caught'])} | {number(row['model_error_at_same_coverage'], 4)} |"
              for row in mapq["all"]]
    for aligner, rows in mapq["by_aligner"].items():
        lines += ["", f"`{aligner}` only:", "", "| Q | Coverage | Error among accepted | Incorrect caught |", "|---|---|---|---|"]
        lines += [f"| {row['threshold']} | {number(row['coverage'])} | {number(row['accepted_error_rate'], 4)} | {number(row['incorrect_caught'])} |" for row in rows]
    return [*lines, ""]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Train and evaluate Model A (alignment reliability) with group-disjoint CV.")
    parser.add_argument("--features", type=Path, default=DEFAULT_DIR / "alignment_features.csv.gz")
    parser.add_argument("--output-dir", type=Path, default=None, help="Default: outputs/ml_alignment/regrouped/<grouping>.")
    parser.add_argument("--grouping", choices=tuple(GROUPINGS), default="window")
    parser.add_argument("--models", nargs="+", choices=tuple(MODEL_FEATURES), default=list(MODEL_FEATURES))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-accept-error", type=float, default=0.01)
    parser.add_argument("--min-reject-precision", type=float, default=0.9)
    args = parser.parse_args(argv)

    rows = load_rows(args.features)
    output_dir = ensure_dir(output_dir_for(args.output_dir, args.grouping).resolve())
    results = [cross_validate(rows, model, n_splits=args.folds, seed=args.seed, max_accept_error=args.max_accept_error,
                              min_reject_precision=args.min_reject_precision, grouping=args.grouping) for model in args.models]
    summaries = [summarize(result) for result in results]
    plots = plot(results, output_dir, args.grouping)
    compared = results[-1]
    y, labelled = compared["y"], compared["rows"]
    mapq_values = np.array([float(row["mapq"]) for row in labelled])
    aligners = np.array([row["aligner"] for row in labelled])
    mapq = {
        "compared_model": compared["model"],
        "all": [{**row, "model_error_at_same_coverage": model_error_at_coverage(y, compared["prob_incorrect"], row["coverage"])}
                for row in mapq_threshold_sweep(y, mapq_values)],
        "by_aligner": {aligner: mapq_threshold_sweep(y[aligners == aligner], mapq_values[aligners == aligner]) for aligner in sorted(set(aligners))},
    }
    write_json(output_dir / "model_results.json", {
        "command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "features_file": project_relative(args.features),
        "label_counts": dict(Counter(row["label"] for row in rows)),
        "targets": {"max_accept_error": args.max_accept_error, "min_reject_precision": args.min_reject_precision},
        "grouping": args.grouping,
        "evaluation": f"{GROUPINGS[args.grouping]}; calibration and thresholds on a group-disjoint quarter of each training fold",
        "mapq_threshold_sweep": mapq,
        "models": summaries,
        "plots": [project_relative(Path(path)) for path in plots],
    })
    write_report(output_dir / "REPORT.md", rows, summaries,
                 {"max_accept_error": args.max_accept_error, "min_reject_precision": args.min_reject_precision}, args.grouping, mapq)
    print(json.dumps({summary["model"]: {key: summary["overall"][key] for key in
                      ("pr_auc_incorrect", "brier", "ece", "coverage", "accepted_error_rate", "incorrect_caught")}
                      for summary in summaries}, indent=2))


if __name__ == "__main__":
    main()
