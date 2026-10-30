"""Model B feasibility study: can the healthy/unhealthy labels be recovered without leakage?

Research label recovery only, not disease association or diagnosis. The source,
version and criteria of the healthy/unhealthy labels are unknown (user decision,
2026-10-30). Each of the 220 labelled sequences is its locus' GRCh38 reference with
one inserted variant, so each example is reduced to that variant: gene, position,
alleles, type, transition/transversion, local GC and CpG context.

Evaluation is StratifiedGroupKFold over 1 kb genomic windows, so neighbouring variants
never sit on both sides of a split, repeated over several seeds. Every model is
compared with:
  majority            the training-fold majority class
  nearest_neighbour   label of the closest training variant in the same gene
  position_only       logistic regression on gene + relative position
and two logistic models, context_only (no position) and context_and_position, which
show how much of any signal is the position confound. Logistic models answer
`uncertain` when |p - 0.5| is below a margin chosen by inner grouped CV on the
training fold (smallest margin reaching `target_accuracy`; otherwise everything is
uncertain).
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, brier_score_loss, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# Runnable as a script: put scripts/ on the import path for common, truth_eval, ml.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, project_relative, read_json, resolve_project_path, write_json
from simulate_truth_reads import read_fasta


GENES = ("APP", "PSEN1", "PSEN2")
WINDOW = 1000
GENE_COLUMNS = [f"gene_{gene}" for gene in GENES]
CONTEXT_COLUMNS = ["is_indel", "abs_length_change", "transition", "gc_fraction", "cpg"]
FEATURE_SETS = {
    "position_only": [*GENE_COLUMNS, "position_fraction"],
    "context_only": [*GENE_COLUMNS, *CONTEXT_COLUMNS],
    "context_and_position": [*GENE_COLUMNS, *CONTEXT_COLUMNS, "position_fraction"],
}
MARGINS = [round(value, 2) for value in np.arange(0.0, 0.5, 0.05)]
DEFAULT_OUTPUT = OUTPUTS_DIR / "ml_variant"
LABEL_SOURCE_NOTE = "Label source, version and criteria are unknown; results describe recovery of these labels, not disease association."


def extract_variant(reference: str, haplotype: str) -> dict:
    if len(reference) == len(haplotype):
        diffs = [index for index, (a, b) in enumerate(zip(reference, haplotype)) if a != b]
        if not diffs:
            raise ValueError("haplotype is identical to the reference")
        first = diffs[0]
        kind = "SNV" if len(diffs) == 1 else "multi_SNV"
        return {"position": first, "ref": reference[first], "alt": haplotype[first], "kind": kind, "length_change": 0, "snv_count": len(diffs)}
    limit = min(len(reference), len(haplotype))
    prefix = 0
    while prefix < limit and reference[prefix] == haplotype[prefix]:
        prefix += 1
    suffix = 0
    while suffix < limit - prefix and reference[-1 - suffix] == haplotype[-1 - suffix]:
        suffix += 1
    ref_allele, alt_allele = reference[prefix:len(reference) - suffix], haplotype[prefix:len(haplotype) - suffix]
    change = len(haplotype) - len(reference)
    kind = "insertion" if not ref_allele else "deletion" if not alt_allele else "complex"
    return {"position": prefix, "ref": ref_allele[:20], "alt": alt_allele[:20], "kind": kind, "length_change": change, "snv_count": 0}


def is_transition(ref: str, alt: str) -> int:
    return int(len(ref) == 1 and len(alt) == 1 and {ref, alt} in ({"A", "G"}, {"C", "T"}))


def in_cpg(reference: str, position: int) -> int:
    base = reference[position]
    return int((base == "C" and reference[position + 1:position + 2] == "G") or (base == "G" and position > 0 and reference[position - 1] == "C"))


def gc_fraction(reference: str, position: int, flank: int = 50) -> float:
    window = reference[max(0, position - flank):position + flank + 1]
    return sum(base in "GC" for base in window) / len(window) if window else 0.0


def window_group(gene: str, position: int, window: int = WINDOW) -> str:
    return f"{gene}:{position // window}"


def build_variant_table() -> list[dict]:
    manifest = read_json(METADATA_DIR / "preprocess_manifest.json")
    test_manifest = read_json(METADATA_DIR / "test_manifest.json")
    test_ids = {entry["id"] for payload in test_manifest["genes"].values() for bucket in ("healthy", "unhealthy") for entry in payload[bucket]}
    table = []
    for gene, payload in manifest["genes"].items():
        reference = read_fasta(resolve_project_path(payload["reference"]["path"]))[0][1]
        for bucket in ("healthy", "unhealthy"):
            for entry in payload[bucket]:
                variant = extract_variant(reference, read_fasta(resolve_project_path(entry["path"]))[0][1])
                position = variant["position"]
                table.append({
                    "id": entry["id"], "gene": gene, "label": bucket.upper(), "split": "test" if entry["id"] in test_ids else "train",
                    **variant, "position_fraction": position / len(reference),
                    "transition": is_transition(variant["ref"], variant["alt"]),
                    "gc_fraction": gc_fraction(reference, min(position, len(reference) - 1)),
                    "cpg": in_cpg(reference, min(position, len(reference) - 1)),
                    "group": window_group(gene, position),
                })
    return table


def feature_matrix(rows: list[dict], names: list[str]) -> np.ndarray:
    def value(row: dict, name: str) -> float:
        if name.startswith("gene_"):
            return float(row["gene"] == name[5:])
        if name == "is_indel":
            return float(row["kind"] in ("insertion", "deletion", "complex"))
        if name == "abs_length_change":
            return float(abs(int(row["length_change"])))
        return float(row[name])
    return np.array([[value(row, name) for name in names] for row in rows], dtype=float)


def grouped_folds(y: np.ndarray, groups: np.ndarray, n_splits: int, seed: int):
    return StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups)


def nearest_neighbour_predict(train: list[dict], test: list[dict]) -> list[int]:
    predictions = []
    for row in test:
        candidates = [item for item in train if item["gene"] == row["gene"]] or train
        nearest = min(candidates, key=lambda item: abs(int(item["position"]) - int(row["position"])))
        predictions.append(int(nearest["label"] == "UNHEALTHY"))
    return predictions


def choose_margin(prob: np.ndarray, y: np.ndarray, *, target_accuracy: float) -> float:
    distance = np.round(np.abs(prob - 0.5), 10)
    for margin in MARGINS:
        confident = distance >= margin
        if confident.any() and ((prob[confident] >= 0.5).astype(int) == y[confident]).mean() >= target_accuracy:
            return margin
    return 0.5


def decide(prob: np.ndarray, margin: float) -> np.ndarray:
    decision = np.full(len(prob), "uncertain", dtype=object)
    if margin < 0.5:
        confident = np.round(np.abs(prob - 0.5), 10) >= margin
        decision[confident & (prob >= 0.5)] = "disease_associated"
        decision[confident & (prob < 0.5)] = "healthy_associated"
    return decision


def _logistic():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))


def _fold_logistic(rows, X, y, groups, train, test, seed, target_accuracy):
    inner_prob = np.zeros(len(train))
    inner_splits = min(3, len(set(groups[train])))
    for inner_train, inner_test in grouped_folds(y[train], groups[train], inner_splits, seed):
        model = _logistic().fit(X[train][inner_train], y[train][inner_train])
        inner_prob[inner_test] = model.predict_proba(X[train][inner_test])[:, 1]
    margin = choose_margin(inner_prob, y[train], target_accuracy=target_accuracy)
    prob = _logistic().fit(X[train], y[train]).predict_proba(X[test])[:, 1]
    return prob, margin


def run_study(table: list[dict], *, seeds: list[int], n_splits: int = 5, target_accuracy: float = 0.8) -> dict:
    y = np.array([int(row["label"] == "UNHEALTHY") for row in table])
    groups = np.array([row["group"] for row in table])
    per_seed: dict[str, list[dict]] = {name: [] for name in ("majority", "nearest_neighbour", *FEATURE_SETS)}
    for seed in seeds:
        hard = {name: np.zeros(len(y), dtype=int) for name in per_seed}
        prob = {name: np.zeros(len(y)) for name in FEATURE_SETS}
        decision = {name: np.empty(len(y), dtype=object) for name in FEATURE_SETS}
        for train, test in grouped_folds(y, groups, n_splits, seed):
            hard["majority"][test] = int(y[train].mean() >= 0.5)
            hard["nearest_neighbour"][test] = nearest_neighbour_predict([table[i] for i in train], [table[i] for i in test])
            for name, columns in FEATURE_SETS.items():
                X = feature_matrix(table, columns)
                fold_prob, margin = _fold_logistic(table, X, y, groups, train, test, seed, target_accuracy)
                prob[name][test] = fold_prob
                hard[name][test] = (fold_prob >= 0.5).astype(int)
                decision[name][test] = decide(fold_prob, margin)
        for name in per_seed:
            metrics = {"balanced_accuracy": balanced_accuracy_score(y, hard[name]), "mcc": matthews_corrcoef(y, hard[name])}
            if name in FEATURE_SETS:
                confident = decision[name] != "uncertain"
                predicted = decision[name] == "disease_associated"
                metrics.update({
                    "roc_auc": roc_auc_score(y, prob[name]), "brier": brier_score_loss(y, prob[name]),
                    "coverage": float(confident.mean()),
                    "confident_rows": int(confident.sum()),
                    "confident_accuracy": float((predicted[confident] == y[confident].astype(bool)).mean()) if confident.any() else None,
                })
            per_seed[name].append(metrics)

    def aggregate(records: list[dict]) -> dict:
        summary = {}
        for key in records[0]:
            values = [record[key] for record in records if record[key] is not None]
            summary[f"{key}_mean"] = statistics.fmean(values) if values else None
            summary[f"{key}_sd"] = statistics.stdev(values) if len(values) > 1 else 0.0
        return summary

    return {"seeds": seeds, "n_splits": n_splits, "target_accuracy": target_accuracy, "window_bp": WINDOW,
            "models": {name: aggregate(records) for name, records in per_seed.items()}}


def random_split_nearest_neighbour(table: list[dict]) -> dict:
    """The original seed-42 per-sequence split, to show how much it leaks through position."""
    train = [row for row in table if row["split"] == "train"]
    test = [row for row in table if row["split"] == "test"]
    predictions = nearest_neighbour_predict(train, test)
    y = [int(row["label"] == "UNHEALTHY") for row in test]
    return {"test_rows": len(test), "correct": int(sum(p == t for p, t in zip(predictions, y))),
            "balanced_accuracy": balanced_accuracy_score(y, predictions)}


def write_report(path: Path, table: list[dict], results: dict, leak: dict) -> None:
    def row(name: str) -> str:
        m = results["models"][name]
        fmt = lambda key: "–" if m.get(f"{key}_mean") is None else f"{m[f'{key}_mean']:.3f} ± {m[f'{key}_sd']:.3f}"
        rows = "–" if m.get("confident_rows_mean") is None else f"{m['confident_rows_mean']:.1f} of {len(table)}"
        return f"| {name} | {fmt('balanced_accuracy')} | {fmt('mcc')} | {fmt('roc_auc')} | {rows} | {fmt('confident_accuracy')} |"

    kinds = Counter(item["kind"] for item in table)
    lines = [
        "# Model B feasibility study (research label recovery, not diagnosis)", "",
        f"**{LABEL_SOURCE_NOTE}**", "",
        f"Data: {len(table)} sequences, each the GRCh38 locus reference with one inserted variant ({dict(kinds)}); "
        f"{sum(item['label'] == 'UNHEALTHY' for item in table)} labelled unhealthy.", "",
        f"Evaluation: StratifiedGroupKFold({results['n_splits']}) over {results['window_bp']} bp genomic windows, seeds {results['seeds']}; "
        f"mean ± SD over seeds. `uncertain` margin chosen by inner grouped CV for target accuracy {results['target_accuracy']}.", "",
        "| Model | Balanced accuracy | MCC | ROC-AUC | Confident predictions (mean) | Accuracy when confident |", "|---|---|---|---|---|---|",
        *[row(name) for name in results["models"]], "",
        f"Original seed-42 per-sequence split, nearest training variant: {leak['correct']}/{leak['test_rows']} correct "
        f"(balanced accuracy {leak['balanced_accuracy']:.3f}); the random split leaks through position.", "",
        "Read `context_only` against `position_only` and `nearest_neighbour`: anything the position models achieve "
        "reflects where the labelled variants were chosen, not a property learned from sequence context.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Model B feasibility study (research label recovery, not diagnosis).")
    parser.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--target-accuracy", type=float, default=0.8)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    output_dir = ensure_dir(args.output_dir.resolve())
    table = build_variant_table()
    with (output_dir / "variant_table.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader(); writer.writerows(table)
    results = run_study(table, seeds=args.seeds, n_splits=args.folds, target_accuracy=args.target_accuracy)
    leak = random_split_nearest_neighbour(table)
    write_json(output_dir / "feasibility_results.json", {
        "command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "label_source": LABEL_SOURCE_NOTE, "rows": len(table), "random_split_nearest_neighbour": leak, **results,
    })
    write_report(output_dir / "REPORT.md", table, results, leak)
    print(json.dumps({name: round(m["balanced_accuracy_mean"], 3) for name, m in results["models"].items()}, indent=2))


if __name__ == "__main__":
    main()
