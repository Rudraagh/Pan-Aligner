"""Model A strict-target evaluation on the v2 dataset.

Asks three questions of alignment_reliability's models, fitted once per grouping:

  abstention   Do stricter risk or coverage targets open a usable insufficient_evidence
               band? A band is `meaningful` only if it opens in every fold, abstains on at
               least MIN_ABSTENTION_SHARE of alignments, and the abstained alignments' error
               rate lies between those of accepted and rejected alignments.
  thresholds   Each training fold is split by group into fit (60%), calibration (20%) and
               threshold (20%) parts, so thresholds see neither the test fold nor the rows
               the calibrator was fitted on.
  MAPQ         Gradient boosting against a plain MAPQ threshold at the same accepted share.
               Isotonic calibration produces tied probabilities, so a tied tier is accepted
               fractionally rather than in row order. Intervals resample whole 20 kb windows.

Within-locus and off-target results are always reported separately. The report states
what the intervals show and does not treat the model as a substitute for MAPQ.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

# Runnable as a script: put scripts/ on the import path for common, truth_eval, ml.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import OUTPUTS_DIR, ensure_dir, project_relative, write_json
from ml.alignment_reliability import DECISIONS, GROUPINGS, apply_targets, cross_validate, expected_calibration_error, load_rows


DEFAULT_OUTPUT = OUTPUTS_DIR / "ml_alignment_v2" / "strict_targets"
DEFAULT_FEATURES = OUTPUTS_DIR / "ml_alignment_v2" / "alignment_features.csv.gz"
MIN_ABSTENTION_SHARE = 0.001
COMPARED_THRESHOLDS = (10, 20, 30, 40, 50, 60)
TARGETS = (
    {"name": "error_0.5pct", "max_accept_error": 0.005, "min_reject_precision": 0.95},
    {"name": "error_0.2pct", "max_accept_error": 0.002, "min_reject_precision": 0.99},
    {"name": "error_0.1pct", "max_accept_error": 0.001, "min_reject_precision": 0.99},
    {"name": "error_0.05pct", "max_accept_error": 0.0005, "min_reject_precision": 0.995},
    {"name": "coverage_95pct", "max_accept_error": 0.01, "min_reject_precision": 0.99, "target_coverage": 0.95},
    {"name": "coverage_90pct", "max_accept_error": 0.01, "min_reject_precision": 0.99, "target_coverage": 0.90},
)


def acceptance_weights(prob: np.ndarray, coverage: float) -> np.ndarray:
    """Accept lowest p(incorrect) first; a tier of tied rows that straddles the cutoff is shared evenly."""
    weights = np.zeros(len(prob))
    remaining = coverage * len(prob)
    values, inverse, counts = np.unique(prob, return_inverse=True, return_counts=True)
    fraction = np.zeros(len(values))
    for index, count in enumerate(counts):
        if remaining <= 0:
            break
        fraction[index] = min(1.0, remaining / count)
        remaining -= count
    weights = fraction[inverse]
    return weights


def weighted_error(y: np.ndarray, weights: np.ndarray) -> float | None:
    total = weights.sum()
    return float((weights * y).sum() / total) if total else None


def assess_abstention(y: np.ndarray, decision: np.ndarray, thresholds: list[dict]) -> dict:
    accepted, abstained, rejected = decision == "correct", decision == "insufficient_evidence", decision == "incorrect"

    def rate(mask):
        return float(y[mask].mean()) if mask.any() else None

    folds_with_band = sum(math.isfinite(t["low"]) and t["low"] < t["high"] for t in thresholds)
    accepted_rate, abstained_rate, rejected_rate = rate(accepted), rate(abstained), rate(rejected)
    meaningful = (
        folds_with_band == len(thresholds)
        and abstained.mean() >= MIN_ABSTENTION_SHARE
        and accepted_rate is not None and abstained_rate is not None
        and abstained_rate > accepted_rate
        and (rejected_rate is None or abstained_rate < rejected_rate)
    )
    return {"abstentions": int(abstained.sum()), "abstention_share": float(abstained.mean()), "folds_with_band": folds_with_band,
            "folds": len(thresholds), "accepted_error_rate": accepted_rate, "abstained_error_rate": abstained_rate,
            "rejected_error_rate": rejected_rate, "meaningful": bool(meaningful)}


def _split_counts(y: np.ndarray, decision: np.ndarray) -> dict:
    return {
        "rows": int(len(y)), "incorrect_rows": int(y.sum()),
        "decisions": {name: int((decision == name).sum()) for name in DECISIONS},
        "accepted_error_rate": float(y[decision == "correct"].mean()) if (decision == "correct").any() else None,
        "incorrect_accepted": int(((decision == "correct") & (y == 1)).sum()),
        "incorrect_abstained": int(((decision == "insufficient_evidence") & (y == 1)).sum()),
        "incorrect_rejected": int(((decision == "incorrect") & (y == 1)).sum()),
        "correct_rejected": int(((decision == "incorrect") & (y == 0)).sum()),
        "correct_abstained": int(((decision == "insufficient_evidence") & (y == 0)).sum()),
    }


def equal_coverage_comparison(y: np.ndarray, prob: np.ndarray, mapq: np.ndarray, decoy: np.ndarray, groups: np.ndarray, *,
                              thresholds=COMPARED_THRESHOLDS, resamples: int = 200, seed: int = 0) -> list[dict]:
    """Model minus MAPQ error at each MAPQ threshold's accepted share, with a window-resampled 95% interval."""
    rng = np.random.default_rng(seed)
    members = {group: np.flatnonzero(groups == group) for group in sorted(set(groups))}
    names = list(members)
    samples = [np.concatenate([members[names[i]] for i in rng.integers(0, len(names), len(names))]) for _ in range(resamples)]
    rows = []
    for threshold in thresholds:
        accepted = mapq >= threshold
        coverage = float(accepted.mean())
        weights = acceptance_weights(prob, coverage)
        differences = []
        for index in samples:
            share = (mapq[index] >= threshold)
            if share.any():
                model = weighted_error(y[index], acceptance_weights(prob[index], float(share.mean())))
                differences.append(model - float(y[index][share].mean()))
        rows.append({
            "threshold": int(threshold), "coverage": coverage,
            "mapq_error": float(y[accepted].mean()) if accepted.any() else None,
            "model_error": weighted_error(y, weights),
            "mapq_accepted_within_locus": int((accepted & (y == 1) & (decoy == 0)).sum()),
            "mapq_accepted_off_target": int((accepted & (y == 1) & (decoy == 1)).sum()),
            "model_accepted_within_locus": float((weights * y * (decoy == 0)).sum()),
            "model_accepted_off_target": float((weights * y * (decoy == 1)).sum()),
            "difference_ci95": (float(np.percentile(differences, 2.5)), float(np.percentile(differences, 97.5))) if differences else (None, None),
        })
    return rows


def mapq_verdict(rows: list[dict]) -> str:
    usable = [row for row in rows if row["difference_ci95"][1] is not None]
    if usable and all(row["difference_ci95"][1] < 0 for row in usable):
        return ("Gradient boosting has lower error than the plain MAPQ threshold at every compared coverage (95% interval below zero). "
                "This holds only for these simulated, locus-scale reads.")
    return ("Gradient boosting does not have reliably lower error than the plain MAPQ threshold, not at every compared coverage, "
            "so these results do not support using it instead of MAPQ.")


def _plot(y, prob, mapq, decoy, output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    subsets = (("All errors", np.ones(len(y), dtype=bool)), ("Within-locus errors", decoy == 0), ("Off-target errors", (y == 0) | (decoy == 1)))
    figure, axes = plt.subplots(1, 3, figsize=(13, 4))
    for axis, (title, keep) in zip(axes, subsets):
        grid = np.linspace(0.5, 1.0, 51)
        curve = [weighted_error(y[keep], acceptance_weights(prob[keep], value)) for value in grid]
        axis.plot(grid, curve, label="gradient boosting")
        points = [((mapq[keep] >= q).mean(), y[keep][mapq[keep] >= q].mean()) for q in range(0, 61, 10) if (mapq[keep] >= q).any()]
        axis.plot(*zip(*points), marker="o", linestyle="--", label="plain MAPQ threshold")
        axis.set(title=title, xlabel="Share of alignments accepted", ylabel="Error among accepted")
        axis.legend(fontsize=8); axis.grid(alpha=0.3)
    figure.tight_layout(); figure.savefig(output_dir / "risk_coverage_split.png", dpi=150); plt.close(figure)

    figure, axis = plt.subplots(figsize=(5.2, 5.2))
    edges = np.array([0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    index = np.clip(np.digitize(prob, edges[1:-1], right=True), 0, len(edges) - 2)
    points = [(prob[index == b].mean(), y[index == b].mean()) for b in range(len(edges) - 1) if (index == b).any()]
    axis.plot(*zip(*points), marker="o", label="gradient boosting")
    axis.plot([1e-4, 1], [1e-4, 1], color="grey", linestyle="--", linewidth=1)
    axis.set(xscale="log", yscale="log", xlabel="Predicted p(incorrect)", ylabel="Observed fraction incorrect", title="Calibration on held-out folds")
    axis.legend(fontsize=8); axis.grid(alpha=0.3)
    figure.tight_layout(); figure.savefig(output_dir / "reliability.png", dpi=150); plt.close(figure)


def _pct(value, digits=2) -> str:
    return "–" if value is None else f"{value * 100:.{digits}f}%"


def _report(results: dict) -> str:
    lines = [
        f"# Model A strict targets: {results['grouping']}", "",
        f"Grouping: {GROUPINGS[results['grouping']]}. Thresholds come only from a validation part of each training fold, disjoint by group from the rows that fitted the model, the rows that fitted the calibrator, and the test fold.", "",
        f"Held-out calibration of gradient boosting: ECE {results['calibration']['ece']:.4f}, Brier {results['calibration']['brier']:.4f}, mean predicted {_pct(results['calibration']['mean_predicted'])} against observed {_pct(results['calibration']['observed'])}.", "",
        "## Targets and abstention", "",
        "A band is meaningful only if it opens in every fold, covers at least "
        f"{MIN_ABSTENTION_SHARE:.1%} of alignments, and its error rate lies between accepted and rejected.", "",
        "| Target | Model | Accepted | Abstentions | Rejected | Error among accepted | Error among abstained | Precision of rejections | Folds with band | Meaningful |", "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in results["targets"]:
        a = row["abstention"]; d = row["overall"]["decisions"]
        lines.append(f"| {row['target']} | {row['model']} | {d['correct']} | {d['insufficient_evidence']} | {d['incorrect']} | {_pct(a['accepted_error_rate'], 3)} | "
                     f"{_pct(a['abstained_error_rate'])} | {_pct(a['rejected_error_rate'], 1)} | {a['folds_with_band']}/{a['folds']} | {'yes' if a['meaningful'] else 'no'} |")
    lines += ["", "## Within-locus and off-target separately (gradient boosting)", "",
              "| Target | Within-locus errors accepted / abstained / rejected | Within-locus correct abstained / rejected | Off-target errors accepted / abstained / rejected |", "|---|---|---|---|"]
    for row in results["targets"]:
        if row["model"] == "gradient_boosting":
            w, o = row["within_locus"], row["off_target"]
            lines.append(f"| {row['target']} | {w['incorrect_accepted']} / {w['incorrect_abstained']} / {w['incorrect_rejected']} | {w['correct_abstained']} / {w['correct_rejected']} | "
                         f"{o['incorrect_accepted']} / {o['incorrect_abstained']} / {o['incorrect_rejected']} |")
    lines += ["", "## Gradient boosting against plain MAPQ thresholds at equal coverage", "",
              "Difference is model minus MAPQ error, with a 95% interval resampling whole 20 kb windows. Accepted errors are split into within-locus and off-target.", "",
              "| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) | MAPQ accepted errors: within / off-target | Model accepted errors: within / off-target |", "|---|---|---|---|---|---|---|"]
    for row in results["mapq_comparison"]:
        low, high = row["difference_ci95"]
        lines.append(f"| {row['threshold']} | {_pct(row['coverage'], 1)} | {_pct(row['mapq_error'], 3)} | {_pct(row['model_error'], 3)} | "
                     f"{row['model_error'] - row['mapq_error']:+.5f} ({low:+.5f} to {high:+.5f}) | {row['mapq_accepted_within_locus']} / {row['mapq_accepted_off_target']} | "
                     f"{row['model_accepted_within_locus']:.0f} / {row['model_accepted_off_target']:.0f} |")
    lines += ["", f"**{results['mapq_comparison_verdict']}**", ""]
    for aligner, rows in results["mapq_comparison_by_aligner"].items():
        lines += [f"`{aligner}` only:", "", "| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) |", "|---|---|---|---|---|"]
        for row in rows:
            low, high = row["difference_ci95"]
            lines.append(f"| {row['threshold']} | {_pct(row['coverage'], 1)} | {_pct(row['mapq_error'], 3)} | {_pct(row['model_error'], 3)} | {low:+.5f} to {high:+.5f} |")
        lines.append("")
    lines += ["## Limits", "",
              "- Simulated PBSIM2 reads at 75% accuracy on three loci; real data and whole-genome graphs are not covered.",
              "- Errors from the same regions repeat, so the effective sample is smaller than the row counts, and there are few independent windows.",
              "- Rankings pool probabilities from different folds' calibrators.",
              "- Targets are compared across this one dataset; no separate untouched test set was held back for choosing among them.", ""]
    return "\n".join(lines)


def run_grouping(rows: list[dict], grouping: str, output_dir: Path, *, n_splits: int = 5, resamples: int = 200, seed: int = 0) -> dict:
    output_dir = ensure_dir(Path(output_dir))
    if (output_dir / "strict_results.json").exists():
        raise FileExistsError(f"{output_dir} already has results; choose another --output-dir.")
    fitted = {model: cross_validate(rows, model, n_splits=n_splits, seed=seed, grouping=grouping, separate_threshold_set=True)
              for model in ("gradient_boosting", "mapq_only")}
    gb = fitted["gradient_boosting"]
    y, prob, labelled, groups = gb["y"], gb["prob_incorrect"], gb["rows"], gb["groups"]
    decoy = np.array([int(row["decoy"]) for row in labelled])
    mapq = np.array([float(row["mapq"]) for row in labelled])
    aligners = np.array([row["aligner"] for row in labelled])
    target_rows = []
    for target in TARGETS:
        settings = {key: target[key] for key in ("max_accept_error", "min_reject_precision")}
        settings["target_coverage"] = target.get("target_coverage")
        for model, result in fitted.items():
            decision, thresholds = apply_targets(result, **settings)
            target_rows.append({
                "target": target["name"], "model": model, **settings, "thresholds": thresholds,
                "abstention": assess_abstention(result["y"], decision, thresholds),
                "overall": _split_counts(result["y"], decision),
                "within_locus": _split_counts(result["y"][decoy == 0], decision[decoy == 0]),
                "off_target": _split_counts(result["y"][decoy == 1], decision[decoy == 1]),
            })
    comparison = equal_coverage_comparison(y, prob, mapq, decoy, groups, resamples=resamples, seed=seed)
    by_aligner = {aligner: equal_coverage_comparison(y[aligners == aligner], prob[aligners == aligner], mapq[aligners == aligner], decoy[aligners == aligner],
                                                     groups[aligners == aligner], resamples=resamples, seed=seed)
                  for aligner in sorted(set(aligners))}
    results = {
        "grouping": grouping, "evaluation": GROUPINGS[grouping], "seed": seed, "resamples": resamples, "rows": len(labelled),
        "incorrect_rows": int(y.sum()), "within_locus_incorrect": int((y * (decoy == 0)).sum()), "off_target_incorrect": int((y * (decoy == 1)).sum()),
        "threshold_selection": "fit 60%, calibration 20%, thresholds 20% of each training fold, disjoint by group; test fold never used",
        "calibration": {"ece": expected_calibration_error(y, prob), "brier": float(np.mean((prob - y) ** 2)), "mean_predicted": float(prob.mean()), "observed": float(y.mean())},
        "targets": target_rows, "mapq_comparison": comparison, "mapq_comparison_by_aligner": by_aligner,
        "mapq_comparison_verdict": mapq_verdict(comparison),
        "fold_groups": [{key: fold[key] for key in ("fold", "test_groups", "calibration_groups", "threshold_groups")} for fold in gb["thresholds"]],
    }
    _plot(y, prob, mapq, decoy, output_dir)
    write_json(output_dir / "strict_results.json", results)
    (output_dir / "REPORT.md").write_text(_report(results), encoding="utf-8")
    return results


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Model A strict-target evaluation with validation-only thresholds.")
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--groupings", nargs="+", choices=("window", "leave_one_gene_out"), default=["window", "leave_one_gene_out"])
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--resamples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    rows = load_rows(args.features)
    summary = {}
    for grouping in args.groupings:
        results = run_grouping(rows, grouping, args.output_dir.resolve() / grouping, n_splits=args.folds, resamples=args.resamples, seed=args.seed)
        summary[grouping] = results["mapq_comparison_verdict"]
    print(json.dumps({"command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]), **summary}, indent=2))


if __name__ == "__main__":
    main()
