"""Train and evaluate the replacement Model B under the protocol fixed in the validation snapshot.

The target is the direction of a ClinVar classification in APP, PSEN1 and PSEN2. Results describe recovery of that
annotation, not diagnosis, individual risk, Alzheimer's association or VUS reclassification.

Primary comparison: missense logistic model minus the best baseline, by balanced accuracy, with folds and clusters
from the snapshot. Models fit on train only. Validation sets settings and the threshold. Test is only scored.

Usage: variant_model_v2.py --output <new folder> [--control permuted]
With --control permuted, labels are shuffled between clusters: a dry run that exposes no real result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import warnings
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from ml.variant_dataset_validation import choose_validation_fold, forbidden_feature_columns  # also sets the pandas string option

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "outputs" / "model_b_dataset_validation" / "20261030_validation_02"
FEATURES = ROOT / "outputs" / "model_b_features" / "20261030_features_01" / "data"
GC_COLUMNS = ("gc_50", "gc_200")
ARMS = ("missense", "missense_no_gc", "amino_only", "main", "main_no_gc")
LOGISTIC_GRID = ({"C": 0.01}, {"C": 0.1}, {"C": 1.0}, {"C": 10.0})
BOOSTING_GRID = ({"max_depth": 2, "max_iter": 50}, {"max_depth": 2, "max_iter": 150}, {"max_depth": 3, "max_iter": 50}, {"max_depth": 3, "max_iter": 150})
N_RESAMPLES = 2000
BOOTSTRAP_SEED = 20261031
BASELINES = {"gene": "gene", "position_vote": "position_vote", "nearest_variant": "nearest_variant"}
SCHEME_PLAN = {
    "missense_primary": {"arms": ("missense", "missense_no_gc", "amino_only"), "models": ("logistic", "gradient_boosting"), "consequence_baseline": False},
    "all_consequences": {"arms": ("main", "main_no_gc"), "models": ("logistic", "gradient_boosting"), "consequence_baseline": True},
    "two_star_all": {"arms": ("main",), "models": ("logistic",), "consequence_baseline": True},
    "psen1_all": {"arms": ("main",), "models": ("logistic",), "consequence_baseline": True},
    "gene_held_out_all": {"arms": ("main",), "models": ("logistic",), "consequence_baseline": True},
    "gene_held_out_missense": {"arms": ("missense",), "models": ("logistic",), "consequence_baseline": False},
    "time_2024_not_region_safe": {"arms": ("main",), "models": ("logistic",), "consequence_baseline": True},
}


def balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    positives, negatives = y == 1, y == 0
    if not positives.any() or not negatives.any():
        return math.nan
    return float(((pred[positives] == 1).mean() + (pred[negatives] == 0).mean()) / 2)


def choose_threshold(scores: np.ndarray, y: np.ndarray) -> float:
    """The cut with the highest balanced accuracy. Ties go to the value nearest 0.5."""
    unique = np.unique(scores)
    candidates = [0.5, float(unique[0]) - 1e-6, float(unique[-1]) + 1e-6, *((unique[:-1] + unique[1:]) / 2).tolist()]
    if len(unique) == 1:
        candidates.append(float(unique[0]))
    best, best_key = 0.5, None
    for candidate in candidates:
        value = balanced_accuracy(y, (scores >= candidate).astype(int))
        key = (round(value, 12), -abs(candidate - 0.5), -candidate)
        if best_key is None or key > best_key:
            best, best_key = candidate, key
    return float(best)


def group_rate_scores(train: pd.DataFrame, target: pd.DataFrame, column: str) -> np.ndarray:
    """Train pathogenic rate per group, smoothed by one half; 0.5 for a group unseen in train."""
    grouped = train.groupby(column)["y"].agg(["sum", "count"])
    rates = (grouped["sum"] + 0.5) / (grouped["count"] + 1)
    return target[column].map(rates).fillna(0.5).to_numpy(dtype=float)


def _gene_lookup(train: pd.DataFrame) -> dict:
    return {gene: (part["pos"].to_numpy(dtype=np.int64), part["y"].to_numpy(dtype=float)) for gene, part in train.groupby("gene")}


def nearest_variant_scores(train: pd.DataFrame, target: pd.DataFrame) -> np.ndarray:
    """Label of the nearest train variant in the same gene; the train rate if the gene has none."""
    lookup, fallback = _gene_lookup(train), float(train["y"].mean())
    scores = []
    for gene, position in zip(target["gene"], target["pos"].astype(np.int64)):
        if gene not in lookup:
            scores.append(fallback)
        else:
            positions, labels = lookup[gene]
            scores.append(float(labels[np.argmin(np.abs(positions - position))]))
    return np.array(scores)


def position_vote_scores(train: pd.DataFrame, target: pd.DataFrame, k: int = 5) -> np.ndarray:
    """Pathogenic share of the k nearest train variants in the same gene; the train rate if the gene has none."""
    lookup, fallback = _gene_lookup(train), float(train["y"].mean())
    scores = []
    for gene, position in zip(target["gene"], target["pos"].astype(np.int64)):
        if gene not in lookup:
            scores.append(fallback)
        else:
            positions, labels = lookup[gene]
            scores.append(float(labels[np.argsort(np.abs(positions - position), kind="stable")[:k]].mean()))
    return np.array(scores)


def arm_columns(arm: str, main_columns: list[str], amino_columns: list[str]) -> list[str]:
    no_gc = [column for column in main_columns if column not in GC_COLUMNS]
    columns = {"missense": [*main_columns, *amino_columns], "missense_no_gc": [*no_gc, *amino_columns], "amino_only": list(amino_columns), "main": list(main_columns), "main_no_gc": no_gc}[arm]
    bad = forbidden_feature_columns(columns)
    if bad:
        raise ValueError(f"forbidden feature columns: {bad}")
    return columns


def _pipeline(frame: pd.DataFrame, columns: list[str], model: str, params: dict, seed: int) -> Pipeline:
    numeric = [column for column in columns if pd.api.types.is_numeric_dtype(frame[column])]
    categorical = [column for column in columns if column not in numeric]
    prepare = ColumnTransformer([
        ("numeric", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric),
        ("categorical", OneHotEncoder(handle_unknown="ignore", sparse_output=False), categorical),
    ])
    if model == "logistic":
        estimator = LogisticRegression(max_iter=5000, class_weight="balanced", **params)
    elif model == "gradient_boosting":
        estimator = HistGradientBoostingClassifier(class_weight="balanced", learning_rate=0.05, min_samples_leaf=5, early_stopping=False, random_state=seed, **params)
    else:
        raise ValueError(model)
    return Pipeline([("prepare", prepare), ("model", estimator)])


def fit_fold(data: pd.DataFrame, role: pd.Series, columns: list[str], model: str, seed: int) -> dict:
    """Fit on train; choose settings and threshold on validation; score test. Test labels are never read."""
    grid = LOGISTIC_GRID if model == "logistic" else BOOSTING_GRID
    train = data[role == "train"].sort_values("variation_id").sample(frac=1, random_state=seed)
    validation, test = data[role == "validation"], data[role == "test"]
    y_validation = validation["y"].to_numpy()
    best = None
    for params in grid:
        pipeline = _pipeline(data, columns, model, params, seed)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            pipeline.fit(train[columns], train["y"].to_numpy())
        scores = pipeline.predict_proba(validation[columns])[:, 1]
        threshold = choose_threshold(scores, y_validation)
        value = balanced_accuracy(y_validation, (scores >= threshold).astype(int))
        if best is None or value > best["validation_ba"] + 1e-12:
            best = {"validation_ba": value, "threshold": threshold, "params": params, "pipeline": pipeline, "validation_scores": scores}
    best["test_scores"] = best["pipeline"].predict_proba(test[columns])[:, 1]
    return best


def choose_baseline(validation_values: dict[str, float]) -> str:
    return sorted(validation_values, key=lambda name: (-validation_values[name], name))[0]


def _cluster_counts(y: np.ndarray, pred: np.ndarray, clusters: np.ndarray, order: np.ndarray) -> np.ndarray:
    counts = np.zeros((len(order), 4))
    index = np.searchsorted(order, clusters)
    np.add.at(counts[:, 0], index, (y == 1) & (pred == 1))
    np.add.at(counts[:, 1], index, y == 1)
    np.add.at(counts[:, 2], index, (y == 0) & (pred == 0))
    np.add.at(counts[:, 3], index, y == 0)
    return counts


def cluster_bootstrap_difference(y: np.ndarray, pred_a: np.ndarray, pred_b: np.ndarray, clusters: np.ndarray, n_resamples: int = N_RESAMPLES, seed: int = BOOTSTRAP_SEED) -> dict:
    """Balanced accuracy of a minus b, with a percentile interval from resampling whole clusters. The same clusters serve both."""
    order = np.unique(clusters)
    a, b = _cluster_counts(y, pred_a, clusters, order), _cluster_counts(y, pred_b, clusters, order)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(order), size=(n_resamples, len(order)))

    def ba(counts: np.ndarray) -> np.ndarray:
        total = counts[draws].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return (total[:, 0] / total[:, 1] + total[:, 2] / total[:, 3]) / 2

    differences = ba(a) - ba(b)
    usable = differences[~np.isnan(differences)]
    return {"difference": balanced_accuracy(y, pred_a) - balanced_accuracy(y, pred_b), "low": float(np.percentile(usable, 2.5)), "high": float(np.percentile(usable, 97.5)), "usable_resamples": int(len(usable)), "clusters": int(len(order))}


def permute_labels_by_cluster(data: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Dry-run control: lay labels back down in a shuffled cluster order. Keeps the class total and label clumping, breaks any tie to features."""
    ordered = data.sort_values(["cluster", "variation_id"])
    clusters = ordered["cluster"].unique()
    shuffled = np.random.default_rng(seed).permutation(clusters)
    labels = np.concatenate([ordered.loc[ordered["cluster"] == cluster, "y"].to_numpy() for cluster in shuffled])
    out = data.copy()
    out.loc[ordered.index, "y"] = labels
    return out


def calibration_error(y: np.ndarray, scores: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0, 1, bins + 1)
    total = 0.0
    for low, high in zip(edges[:-1], edges[1:]):
        chosen = (scores >= low) & (scores <= high if high == 1 else scores < high)
        if chosen.any():
            total += chosen.mean() * abs(scores[chosen].mean() - y[chosen].mean())
    return float(total)


def metric_row(y: np.ndarray, scores: np.ndarray, pred: np.ndarray) -> dict:
    two_classes = len(set(y.tolist())) == 2
    tp, fn, tn, fp = int(((y == 1) & (pred == 1)).sum()), int(((y == 1) & (pred == 0)).sum()), int(((y == 0) & (pred == 0)).sum()), int(((y == 0) & (pred == 1)).sum())
    return {"n": len(y), "n_pathogenic": int((y == 1).sum()), "n_benign": int((y == 0).sum()), "balanced_accuracy": balanced_accuracy(y, pred),
            "sensitivity": tp / (tp + fn) if tp + fn else math.nan, "specificity": tn / (tn + fp) if tn + fp else math.nan, "mcc": float(matthews_corrcoef(y, pred)) if two_classes else math.nan,
            "roc_auc": float(roc_auc_score(y, scores)) if two_classes and len(set(scores.tolist())) > 1 else math.nan,
            "pr_auc_pathogenic": float(average_precision_score(y, scores)) if two_classes else math.nan, "pr_auc_benign": float(average_precision_score(1 - y, -scores)) if two_classes else math.nan,
            "brier": float(brier_score_loss(y, np.clip(scores, 0, 1))), "calibration_error": calibration_error(y, np.clip(scores, 0, 1))}


def _with_time_validation(role: pd.Series, data: pd.DataFrame, seed: int) -> pd.Series:
    """The time scheme has no validation role: one grouped fold of its train rows takes that job."""
    train = data[role == "train"]
    fold_of = np.empty(len(train), dtype=int)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="The least populated class")
        for index, (_, test) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(np.zeros(len(train)), train["y"], train["cluster"])):
            fold_of[test] = index
    chosen = choose_validation_fold(fold_of, train["y"].to_numpy())
    out = role.copy()
    out.loc[train.index[fold_of == chosen]] = "validation"
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(frame: pd.DataFrame, path: Path) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, index=False, lineterminator="\n")


def run(output: Path, control: str | None, seeds: list[int] | None) -> None:
    output.mkdir(parents=True, exist_ok=False)
    variants = pd.read_csv(SNAPSHOT / "CORRECTED_VARIANTS.csv", dtype=str, keep_default_na=False)
    main_table = pd.read_csv(FEATURES / "features_main.csv", dtype={"variation_id": str})
    amino_table = pd.read_csv(FEATURES / "features_missense.csv", dtype={"variation_id": str})
    manifest = pd.read_csv(SNAPSHOT / "SPLIT_MANIFEST.csv", dtype={"variation_id": str, "held_out_gene": str}, keep_default_na=False)
    main_columns = [column for column in main_table.columns if column != "variation_id"]
    amino_columns = [column for column in amino_table.columns if column != "variation_id"]
    base = variants[["variation_id", "gene", "pos", "cluster", "consequence_group", "label"]].copy()
    base["pos"] = base["pos"].astype(np.int64)
    base["y"] = (base["label"] == "pathogenic").astype(int)
    if control == "permuted":
        base = permute_labels_by_cluster(base, seed=BOOTSTRAP_SEED)
    base = base.drop(columns="label").merge(main_table, on="variation_id", validate="one_to_one").merge(amino_table, on="variation_id", how="left", validate="one_to_one")

    predictions, choices, coefficients = [], [], []
    for scheme, plan in SCHEME_PLAN.items():
        scheme_rows = manifest[manifest["scheme"] == scheme]
        if seeds is not None:
            scheme_rows = scheme_rows[scheme_rows["seed"].isin(seeds)]
        for (seed, fold, held), part in scheme_rows.groupby(["seed", "fold", "held_out_gene"]):
            data = base.merge(part[["variation_id", "role"]], on="variation_id")
            data = data[data["role"].isin(("train", "validation", "test"))].reset_index(drop=True)
            role = data["role"].copy()
            if "validation" not in set(role):
                role = _with_time_validation(role, data, int(seed))
            train, validation, test = data[role == "train"], data[role == "validation"], data[role == "test"]
            y_validation, y_test = validation["y"].to_numpy(), test["y"].to_numpy()
            ids = test["variation_id"].to_numpy()

            def record(arm: str, model: str, scores: np.ndarray, threshold: float, params: str, validation_ba: float) -> None:
                predictions.append(pd.DataFrame({"scheme": scheme, "seed": seed, "fold": fold, "held_out_gene": held, "arm": arm, "model": model, "variation_id": ids, "score": scores, "threshold": threshold, "prediction": (scores >= threshold).astype(int)}))
                choices.append({"scheme": scheme, "seed": seed, "fold": fold, "held_out_gene": held, "arm": arm, "model": model, "settings": params, "threshold": threshold, "validation_balanced_accuracy": validation_ba})

            functions = {"gene": lambda t, q: group_rate_scores(t, q, "gene"), "position_vote": position_vote_scores, "nearest_variant": nearest_variant_scores}
            if plan["consequence_baseline"]:
                functions["consequence"] = lambda t, q: group_rate_scores(t, q, "consequence_group")
            fitted, values = {}, {}
            for name, function in functions.items():
                validation_scores = function(train, validation)
                threshold = choose_threshold(validation_scores, y_validation)
                values[name] = balanced_accuracy(y_validation, (validation_scores >= threshold).astype(int))
                fitted[name] = (function(train, test), threshold)
                record("baseline", name, fitted[name][0], threshold, "", values[name])
            chosen = choose_baseline(values)
            record("baseline", "best_baseline", fitted[chosen][0], fitted[chosen][1], chosen, values[chosen])

            for arm in plan["arms"]:
                columns = arm_columns(arm, main_columns, amino_columns)
                for model in plan["models"]:
                    result = fit_fold(data, role, columns, model, int(seed))
                    record(arm, model, result["test_scores"], result["threshold"], json.dumps(result["params"], sort_keys=True), result["validation_ba"])
                    if scheme == "missense_primary" and model == "logistic":
                        names = result["pipeline"].named_steps["prepare"].get_feature_names_out()
                        for name, value in zip(names, result["pipeline"].named_steps["model"].coef_[0]):
                            coefficients.append({"arm": arm, "seed": seed, "fold": fold, "term": name, "coefficient": float(value)})
    prediction_table = pd.concat(predictions, ignore_index=True).merge(base[["variation_id", "gene", "cluster", "y"]], on="variation_id")

    # Metrics.
    metric_rows = []
    keys = ["scheme", "arm", "model", "seed", "held_out_gene"]
    for key, part in prediction_table.groupby(keys):
        metric_rows.append({**dict(zip(keys, key)), "group": "all", **metric_row(part["y"].to_numpy(), part["score"].to_numpy(), part["prediction"].to_numpy())})
        for gene, gene_part in part.groupby("gene"):
            metric_rows.append({**dict(zip(keys, key)), "group": gene, **metric_row(gene_part["y"].to_numpy(), gene_part["score"].to_numpy(), gene_part["prediction"].to_numpy())})
    metrics = pd.DataFrame(metric_rows)
    numeric = ["balanced_accuracy", "sensitivity", "specificity", "mcc", "roc_auc", "pr_auc_pathogenic", "pr_auc_benign", "brier", "calibration_error"]
    summary = metrics.groupby(["scheme", "arm", "model", "held_out_gene", "group"])[numeric + ["n", "n_pathogenic", "n_benign"]].mean().round(4).reset_index()
    summary["seeds"] = metrics.groupby(["scheme", "arm", "model", "held_out_gene", "group"])["seed"].nunique().to_numpy()
    summary["balanced_accuracy_sd_across_seeds"] = metrics.groupby(["scheme", "arm", "model", "held_out_gene", "group"])["balanced_accuracy"].std().round(4).to_numpy()

    # Comparisons.
    def pair(scheme: str, arm_a: str, model_a: str, arm_b: str, model_b: str, seed: int) -> dict:
        part = prediction_table[(prediction_table["scheme"] == scheme) & (prediction_table["seed"] == seed)]
        a = part[(part["arm"] == arm_a) & (part["model"] == model_a)].set_index("variation_id").sort_index()
        b = part[(part["arm"] == arm_b) & (part["model"] == model_b)].set_index("variation_id").sort_index()
        if not a.index.equals(b.index):
            raise ValueError("comparison rows differ")
        return cluster_bootstrap_difference(a["y"].to_numpy(), a["prediction"].to_numpy(), b["prediction"].to_numpy(), a["cluster"].to_numpy())

    available_seeds = sorted(prediction_table["seed"].unique().tolist())
    primary = {str(seed): pair("missense_primary", "missense", "logistic", "baseline", "best_baseline", seed) for seed in available_seeds}
    headline_seed = available_seeds[0]
    secondary = {}
    for name, (scheme, arm_a, model_a, arm_b, model_b) in {
        "missense_boosting_minus_best_baseline": ("missense_primary", "missense", "gradient_boosting", "baseline", "best_baseline"),
        "missense_no_gc_logistic_minus_best_baseline": ("missense_primary", "missense_no_gc", "logistic", "baseline", "best_baseline"),
        "missense_no_gc_boosting_minus_best_baseline": ("missense_primary", "missense_no_gc", "gradient_boosting", "baseline", "best_baseline"),
        "missense_logistic_minus_no_gc": ("missense_primary", "missense", "logistic", "missense_no_gc", "logistic"),
        "amino_only_logistic_minus_best_baseline": ("missense_primary", "amino_only", "logistic", "baseline", "best_baseline"),
        "main_logistic_minus_best_baseline": ("all_consequences", "main", "logistic", "baseline", "best_baseline"),
        "main_logistic_minus_consequence": ("all_consequences", "main", "logistic", "baseline", "consequence"),
        "main_boosting_minus_consequence": ("all_consequences", "main", "gradient_boosting", "baseline", "consequence"),
        "main_no_gc_logistic_minus_consequence": ("all_consequences", "main_no_gc", "logistic", "baseline", "consequence"),
    }.items():
        secondary[name] = pair(scheme, arm_a, model_a, arm_b, model_b, headline_seed)
    primary_mean = float(metrics[(metrics["scheme"] == "missense_primary") & (metrics["arm"] == "missense") & (metrics["model"] == "logistic") & (metrics["group"] == "all")]["balanced_accuracy"].mean())
    baseline_mean = float(metrics[(metrics["scheme"] == "missense_primary") & (metrics["model"] == "best_baseline") & (metrics["group"] == "all")]["balanced_accuracy"].mean())

    write_new(prediction_table.drop(columns=["y"]).assign(score=lambda d: d["score"].round(6), threshold=lambda d: d["threshold"].round(6)), output / "predictions.csv")
    write_new(pd.DataFrame(choices).round(6), output / "fold_choices.csv")
    write_new(metrics.round(6), output / "metrics_by_seed.csv")
    write_new(summary, output / "metrics_summary.csv")
    coefficient_table = pd.DataFrame(coefficients)
    coefficient_summary = coefficient_table.groupby(["arm", "term"])["coefficient"].agg(mean="mean", sd="std", share_positive=lambda s: float((s > 0).mean())).round(4).reset_index()
    write_new(coefficient_summary, output / "logistic_coefficients_missense.csv")
    results = {
        "control": control, "models_trained_on": "train role only", "primary": {"comparison": "missense logistic minus best baseline, balanced accuracy", "model_mean_across_seeds": primary_mean, "best_baseline_mean_across_seeds": baseline_mean, "difference_of_means": primary_mean - baseline_mean, "by_seed": primary, "headline_seed": headline_seed},
        "secondary_headline_seed": secondary, "parameters": {"logistic_grid": LOGISTIC_GRID, "boosting_grid": BOOSTING_GRID, "resamples": N_RESAMPLES, "bootstrap_seed": BOOTSTRAP_SEED, "learning_rate": 0.05, "min_samples_leaf": 5},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scikit_learn": sklearn.__version__},
        "inputs_sha256": {path.name: sha256(path) for path in (SNAPSHOT / "CORRECTED_VARIANTS.csv", SNAPSHOT / "SPLIT_MANIFEST.csv", FEATURES / "features_main.csv", FEATURES / "features_missense.csv")},
    }
    results["outputs_sha256"] = {path.name: sha256(path) for path in sorted(output.glob("*.csv"))}
    with (output / "results.json").open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(results, handle, indent=1)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True, help="A new folder; existing folders are refused.")
    parser.add_argument("--control", choices=("permuted",))
    parser.add_argument("--seeds", type=int, nargs="+", help="Limit seeds; for dry runs only.")
    args = parser.parse_args(argv)
    run(args.output, args.control, args.seeds)


if __name__ == "__main__":
    main()
