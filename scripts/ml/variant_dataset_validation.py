"""Corrected Model B dataset snapshot: exclusions, genomic clusters, split manifests, leakage checks, feature audit.

Reads the frozen label table and the built features; both are left untouched. No model is fitted. The feature
audit does read labels, but only to report screening associations: nothing is selected or removed because of them.

Genomic clusters replace the freeze's fixed windows. Two variants within MAX_GAP bases, or on the same protein
residue, are always in one cluster, so no fixed-window edge can separate neighbours and no rows need dropping.

Usage: variant_dataset_validation.py --output <new folder>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.model_selection import StratifiedGroupKFold

# Pandas 3 stores strings with pyarrow by default, which crashes this Windows Python once the full test suite's
# native libraries are loaded. Plain Python strings give the same results here.
pd.options.mode.string_storage = "python"

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "outputs" / "model_b_task1_freeze" / "20261030_freeze_01" / "data"
FEATURES = ROOT / "outputs" / "model_b_features" / "20261030_features_01" / "data"
AUDIT_DERIVED = ROOT / "outputs" / "dataset_replacement_feasibility" / "20261030_audit_01" / "derived"
GENES = ("APP", "PSEN1", "PSEN2")
LABELS = ("pathogenic", "benign")
MAX_GAP = 1000
SEEDS = (0, 1, 2, 3, 4)
N_FOLDS = 5
PROXY_THRESHOLD = 0.5
SPLIT_ROLES = ("train", "validation", "test")
LABEL_FOR_CLASS = {"Pathogenic": "pathogenic", "Likely pathogenic": "pathogenic", "Pathogenic/Likely pathogenic": "pathogenic",
                   "Benign": "benign", "Likely benign": "benign", "Benign/Likely benign": "benign"}
FORBIDDEN_FEATURES = {
    "label", "label_code", "aggregate_class", "review_status", "stars", "number_of_submitters", "sole_submitter", "opposing_scv", "vus_scv",
    "last_evaluated", "first_seen_release", "name", "vcv", "vcv_version", "rcvs", "canonical_spdi", "condition", "conditions", "submitter", "in_legacy_220", "clinvar_snapshot",
    "gene", "chrom", "pos", "relative_position", "consequence", "consequence_group", "vep_consequence", "protein_position", "residue", "window_1kb", "window_20kb", "residue_group", "cluster", "group_id",
}


def invalid_reasons(row: dict) -> list[str]:
    reasons = []
    ref, alt = row.get("ref", ""), row.get("alt", "")
    for name, allele in (("REF", ref), ("ALT", alt)):
        if not allele or set(allele) - set("ACGT"):
            reasons.append(f"{name} is not a plain nucleotide allele")
    if ref and ref == alt:
        reasons.append("REF equals ALT: reference-allele record, not a variant")
    if row.get("label") not in LABELS:
        reasons.append("label is not pathogenic or benign")
    if row.get("gene") not in GENES:
        reasons.append("gene is not APP, PSEN1 or PSEN2")
    vcv = row.get("vcv", "")
    if vcv[3:].lstrip("0") != row.get("variation_id"):
        reasons.append("accession does not match Variation ID")
    return reasons


def cluster_variants(frame: pd.DataFrame, max_gap: int = MAX_GAP) -> pd.Series:
    """Cluster IDs aligned to `frame`. Named by gene and first position, so they do not depend on row order."""
    parent: dict[str, str] = {identifier: identifier for identifier in frame["variation_id"]}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(a: str, b: str) -> None:
        parent[find(a)] = find(b)

    starts = frame["pos"].astype(int)
    ends = starts + frame["ref"].str.len() - 1
    table = frame.assign(var_start=starts, var_end=ends)
    for _, part in table.groupby("gene"):
        previous_id, furthest = None, None
        for row in part.sort_values(["var_start", "var_end", "variation_id"]).itertuples():
            if previous_id is not None and row.var_start - furthest <= max_gap:
                union(previous_id, row.variation_id)
            furthest = row.var_end if furthest is None else max(furthest, row.var_end)
            previous_id = row.variation_id
    with_residue = table[table["residue"] != ""]
    for _, part in with_residue.groupby(["gene", "residue"]):
        identifiers = part["variation_id"].tolist()
        for other in identifiers[1:]:
            union(identifiers[0], other)
    roots = table["variation_id"].map(find)
    first = table.groupby(roots)["var_start"].transform("min")
    return table["gene"] + ":" + first.astype(str)


def _stratify_key(frame: pd.DataFrame) -> pd.Series:
    return frame["gene"] + "|" + frame["label"]


def _folds(frame: pd.DataFrame, seed: int, n_folds: int) -> np.ndarray:
    fold_of = np.empty(len(frame), dtype=int)
    splitter = StratifiedGroupKFold(n_folds, shuffle=True, random_state=seed)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="The least populated class")
        for index, (_, test) in enumerate(splitter.split(np.zeros(len(frame)), _stratify_key(frame), frame["cluster"])):
            fold_of[test] = index
    return fold_of


def fold_manifest(frame: pd.DataFrame, scheme: str, seeds, n_folds: int = N_FOLDS) -> pd.DataFrame:
    """Each fold is test once; the next fold validates; the rest train."""
    rows = []
    for seed in seeds:
        fold_of = _folds(frame, seed, n_folds)
        for fold in range(n_folds):
            roles = np.full(len(frame), "train", dtype=object)
            roles[fold_of == (fold + 1) % n_folds] = "validation"
            roles[fold_of == fold] = "test"
            rows.append(pd.DataFrame({"scheme": scheme, "seed": seed, "fold": fold, "held_out_gene": "", "variation_id": frame["variation_id"].to_numpy(), "role": roles, "group_id": frame["cluster"].to_numpy()}))
    return pd.concat(rows, ignore_index=True)


def choose_validation_fold(fold_of: np.ndarray, labels: np.ndarray) -> int:
    """The first fold holding both classes, since a one-class validation set cannot set a threshold. Fold 0 if none does."""
    for fold in sorted(set(fold_of.tolist())):
        if len(set(labels[fold_of == fold].tolist())) == 2:
            return fold
    return 0


def gene_held_out_manifest(frame: pd.DataFrame, scheme: str, seeds, n_folds: int = N_FOLDS) -> pd.DataFrame:
    """One gene is the test set. One grouped fold of the remaining genes validates, and a new seed redraws it."""
    rows = []
    for held_index, gene in enumerate(sorted(frame["gene"].unique())):
        held, rest = frame[frame["gene"] == gene], frame[frame["gene"] != gene]
        for seed in seeds:
            fold_of = _folds(rest, seed, n_folds)
            roles = np.where(fold_of == choose_validation_fold(fold_of, rest["label"].to_numpy()), "validation", "train")
            for part, part_roles in ((held, np.full(len(held), "test", dtype=object)), (rest, roles)):
                rows.append(pd.DataFrame({"scheme": scheme, "seed": seed, "fold": held_index, "held_out_gene": gene, "variation_id": part["variation_id"].to_numpy(), "role": part_roles, "group_id": part["cluster"].to_numpy()}))
    return pd.concat(rows, ignore_index=True)


def _fold_keys(manifest: pd.DataFrame) -> list[str]:
    return [column for column in ("scheme", "seed", "fold", "held_out_gene") if column in manifest]


def nearby_crossings(frame: pd.DataFrame, part: pd.DataFrame, max_gap: int) -> int:
    """Pairs adjacent in genomic order, within max_gap, in different split roles. Dropped rows do not count."""
    used = part[part["role"].isin(SPLIT_ROLES)].merge(frame[["variation_id", "gene", "pos", "ref"]], on="variation_id")
    count = 0
    for _, gene_part in used.groupby("gene"):
        starts = gene_part["pos"].astype(int)
        ordered = gene_part.assign(var_s=starts, var_e=starts + gene_part["ref"].str.len() - 1).sort_values(["var_s", "var_e"])
        previous_role, furthest = None, None
        for row in ordered.itertuples():
            if previous_role is not None and row.var_s - furthest <= max_gap and row.role != previous_role:
                count += 1
            furthest = row.var_e if furthest is None else max(furthest, row.var_e)
            previous_role = row.role
    return count


def audit_manifest(frame: pd.DataFrame, manifest: pd.DataFrame, max_gap: int = MAX_GAP) -> list[str]:
    """Every leakage problem in one manifest, as text. Empty means clean."""
    failures = []
    unknown = set(manifest["variation_id"]) - set(frame["variation_id"])
    if unknown:
        failures.append(f"{len(unknown)} variants are not in the dataset, e.g. {sorted(unknown)[:3]}")
    known = manifest[manifest["variation_id"].isin(frame["variation_id"])]
    for key, part in known.groupby(_fold_keys(known)):
        duplicates = part["variation_id"].duplicated().sum()
        if duplicates:
            failures.append(f"{key}: {duplicates} variants appear more than once")
        used = part[part["role"].isin(SPLIT_ROLES)]
        if "cluster" in frame:
            spread = used.merge(frame[["variation_id", "cluster"]], on="variation_id").groupby("cluster")["role"].nunique()
            if (spread > 1).any():
                failures.append(f"{key}: {(spread > 1).sum()} clusters cross roles")
        residues = used.merge(frame[["variation_id", "gene", "residue"]], on="variation_id")
        residues = residues[residues["residue"] != ""]
        crossing = residues.groupby(["gene", "residue"])["role"].nunique()
        if (crossing > 1).any():
            failures.append(f"{key}: {(crossing > 1).sum()} protein residues cross roles")
        pairs = nearby_crossings(frame, part, max_gap)
        if pairs:
            failures.append(f"{key}: {pairs} neighbouring pairs within {max_gap} bp cross roles")
    return failures


def forbidden_feature_columns(columns) -> list[str]:
    return [column for column in columns if column in FORBIDDEN_FEATURES]


def join_problems(expected, found, table: str, subset: bool = False) -> list[str]:
    """Problems joining a feature table to the variant list one-to-one. `subset` allows absent variants."""
    problems = []
    found = list(found)
    if len(found) != len(set(found)):
        problems.append(f"{table}: duplicate rows")
    missing, extra = set(expected) - set(found), set(found) - set(expected)
    if missing and not subset:
        problems.append(f"{table}: {len(missing)} variants have no row")
    if extra:
        problems.append(f"{table}: {len(extra)} rows have no variant")
    return problems


def _auc(values: pd.Series, positive: pd.Series) -> float:
    """Rank AUC: how often the positive group has the higher value. NaN if one group is empty."""
    keep = values.notna()
    ranks = values[keep].rank()
    flags = positive[keep].astype(bool)
    n_pos, n_neg = int(flags.sum()), int((~flags).sum())
    if not n_pos or not n_neg:
        return math.nan
    return round(float((ranks[flags].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)), 4)


def _cramers_v(a: pd.Series, b: pd.Series) -> float:
    table = pd.crosstab(a, b)
    if min(table.shape) < 2:
        return math.nan
    observed = table.to_numpy(dtype=float)
    expected = observed.sum(axis=1, keepdims=True) * observed.sum(axis=0, keepdims=True) / observed.sum()
    chi2 = float(((observed - expected) ** 2 / expected).sum())
    return round(math.sqrt(chi2 / (observed.sum() * (min(table.shape) - 1))), 4)


def _eta(values: pd.Series, groups: pd.Series) -> float:
    """Correlation ratio: the share of a numeric feature's spread explained by a grouping, as a root."""
    keep = values.notna()
    values, groups = values[keep].astype(float), groups[keep]
    total = float(((values - values.mean()) ** 2).sum())
    if not total:
        return 0.0
    between = sum(len(part) * (part.mean() - values.mean()) ** 2 for _, part in values.groupby(groups))
    return round(math.sqrt(between / total), 4)


def feature_audit(features: pd.DataFrame, context: pd.DataFrame, table: str) -> pd.DataFrame:
    """Screening associations for each feature. A flag is a prompt to investigate, never a reason to remove."""
    joined = features.merge(context, on="variation_id", validate="one_to_one")
    is_pathogenic = joined["label"] == "pathogenic"
    is_missense = joined["consequence_group"] == "missense"
    rows = []
    for column in features.columns:
        if column == "variation_id":
            continue
        values = joined[column]
        numeric = pd.api.types.is_numeric_dtype(values)
        counts = values.value_counts(dropna=True)
        row = {"table": table, "feature": column, "kind": "numeric" if numeric else "categorical", "n": len(values), "n_missing": int(values.isna().sum()),
               "n_unique": int(values.nunique()), "top_value_share": round(float(counts.iloc[0] / len(values)), 4) if len(counts) else math.nan}
        if numeric:
            label_auc, missense_auc = _auc(values, is_pathogenic), _auc(values, is_missense)
            within = _auc(values[is_missense], is_pathogenic[is_missense])
            psen1 = joined["gene"] == "PSEN1"
            row.update(label_auc=label_auc, missense_auc=missense_auc, label_auc_within_missense=within, label_auc_within_psen1=_auc(values[psen1], is_pathogenic[psen1]),
                       label_association=round(abs(2 * label_auc - 1), 4) if not math.isnan(label_auc) else math.nan,
                       consequence_association=_eta(values, joined["consequence_group"]), gene_association=_eta(values, joined["gene"]))
        else:
            known = values.notna()
            row.update(label_auc=math.nan, missense_auc=math.nan, label_auc_within_missense=math.nan, label_auc_within_psen1=math.nan,
                       label_association=_cramers_v(values[known], joined["label"][known]), consequence_association=_cramers_v(values[known], joined["consequence_group"][known]),
                       gene_association=_cramers_v(values[known], joined["gene"][known]))
        flags = []
        if row["n_unique"] <= 1:
            flags.append("constant")
        if row["n_missing"]:
            flags.append("has_missing")
        if row["n_unique"] > 1 and row["top_value_share"] >= 0.95:
            flags.append("nearly_constant")
        for name, field in (("possible_label_proxy", "label_association"), ("possible_consequence_proxy", "consequence_association"), ("possible_gene_proxy", "gene_association")):
            if not math.isnan(row[field]) and row[field] >= PROXY_THRESHOLD:
                flags.append(name)
        row["flags"] = ";".join(flags)
        rows.append(row)
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(frame: pd.DataFrame, path: Path) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, index=False, lineterminator="\n")


def split_summary(manifest: pd.DataFrame, frame: pd.DataFrame) -> pd.DataFrame:
    merged = manifest.merge(frame[["variation_id", "gene", "label"]], on="variation_id")
    counts = merged.groupby(["scheme", "seed", "fold", "held_out_gene", "role", "gene", "label"]).size().rename("n").reset_index()
    return counts


def fold_usability(manifest: pd.DataFrame, frame: pd.DataFrame) -> pd.DataFrame:
    """Per fold, the smallest class count in each role; a fold is usable only if every role has both classes."""
    merged = manifest.merge(frame[["variation_id", "label"]], on="variation_id")
    rows = []
    for (scheme, seed, fold, held), part in merged.groupby(["scheme", "seed", "fold", "held_out_gene"]):
        row = {"scheme": scheme, "seed": seed, "fold": fold, "held_out_gene": held}
        usable = True
        for role in ("train", "validation", "test"):
            roles = part[part["role"] == role]
            if roles.empty:
                continue
            pathogenic, benign = int((roles["label"] == "pathogenic").sum()), int((roles["label"] == "benign").sum())
            row.update({f"{role}_pathogenic": pathogenic, f"{role}_benign": benign, f"{role}_groups": roles["group_id"].nunique()})
            usable = usable and pathogenic > 0 and benign > 0
        row["usable"] = usable
        rows.append(row)
    return pd.DataFrame(rows)


def build(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=False)
    read = lambda path: pd.read_csv(path, dtype=str, keep_default_na=False)
    freeze = read(FREEZE / "task1_labels.csv")
    results: dict = {"freeze_rows": len(freeze)}

    # Exclusions: found by rule, then required to match the feature build's list.
    excluded_rows, reasons_by_id = [], {}
    for row in freeze.to_dict("records"):
        reasons = invalid_reasons(row)
        wrong_label = LABEL_FOR_CLASS.get(row["aggregate_class"]) != row["label"]
        if wrong_label:
            reasons.append("label does not follow from the aggregate class")
        if reasons:
            reasons_by_id[row["variation_id"]] = "; ".join(reasons)
            excluded_rows.append({"variation_id": row["variation_id"], "vcv": row["vcv"], "name": row["name"], "gene": row["gene"], "chrom": row["chrom"], "pos": row["pos"], "ref": row["ref"], "alt": row["alt"],
                                  "original_label": row["label"], "stars": row["stars"], "in_legacy_220": row["in_legacy_220"], "reason": reasons_by_id[row["variation_id"]], "source": "task1_labels.csv, freeze 20261030_freeze_01"})
    feature_excluded = set(read(FEATURES / "excluded_variants.csv")["variation_id"])
    if feature_excluded != set(reasons_by_id):
        raise SystemExit(f"exclusions differ: rule {sorted(reasons_by_id)}, feature build {sorted(feature_excluded)}")

    corrected = freeze[~freeze["variation_id"].isin(reasons_by_id)].copy()
    corrected["cluster"] = cluster_variants(corrected, MAX_GAP)
    meta = read(AUDIT_DERIVED / "record_submitter_meta.csv.gz")[["variation_id", "germline_last_evaluated"]]
    corrected = corrected.merge(meta, on="variation_id", how="left", validate="one_to_one")
    corrected["data_note"] = np.where(corrected["consequence"] == "", "ClinVar lists no molecular consequence", "")
    write_new(corrected, output / "CORRECTED_VARIANTS.csv")
    write_new(pd.DataFrame(excluded_rows), output / "EXCLUDED_VARIANTS.csv")

    # Joins.
    main_table = pd.read_csv(FEATURES / "features_main.csv", dtype={"variation_id": str})
    missense_table = pd.read_csv(FEATURES / "features_missense.csv", dtype={"variation_id": str})
    baselines = read(FEATURES / "baselines.csv")
    ids = corrected["variation_id"]
    joins = join_problems(ids, main_table["variation_id"], "features_main") + join_problems(ids, baselines["variation_id"], "baselines") + join_problems(ids, missense_table["variation_id"], "features_missense", subset=True)
    missense_ids = set(corrected.loc[corrected["consequence_group"] == "missense", "variation_id"])
    if set(missense_table["variation_id"]) != missense_ids:
        joins.append("features_missense: rows differ from the missense set")
    results["join_problems"] = joins
    feature_manifest = json.loads((FEATURES / "feature_manifest.json").read_text(encoding="utf-8"))
    results["feature_build_used_this_freeze"] = feature_manifest["freeze_sha256"] == sha256(FREEZE / "task1_labels.csv")

    # Splits.
    missense = corrected[corrected["consequence_group"] == "missense"]
    manifests = [
        fold_manifest(corrected, "all_consequences", SEEDS),
        fold_manifest(missense, "missense_primary", SEEDS),
        fold_manifest(corrected[corrected["stars"] == "2"], "two_star_all", SEEDS),
        fold_manifest(corrected[corrected["gene"] == "PSEN1"], "psen1_all", SEEDS),
        gene_held_out_manifest(corrected, "gene_held_out_all", SEEDS),
        gene_held_out_manifest(missense, "gene_held_out_missense", SEEDS),
    ]
    time_roles = read(FREEZE / "splits" / "time_2024_all.csv")
    time_roles = time_roles[~time_roles["variation_id"].isin(reasons_by_id)]
    clusters = corrected.set_index("variation_id")["cluster"]
    time_manifest = pd.DataFrame({"scheme": "time_2024_not_region_safe", "seed": 0, "fold": 0, "held_out_gene": "", "variation_id": time_roles["variation_id"].to_numpy(), "role": time_roles["role"].to_numpy(), "group_id": time_roles["variation_id"].map(clusters).to_numpy()})
    manifest = pd.concat([*manifests, time_manifest], ignore_index=True)
    write_new(manifest, output / "SPLIT_MANIFEST.csv")
    write_new(split_summary(manifest, corrected), output / "split_counts.csv")
    usability = fold_usability(manifest, corrected)
    write_new(usability, output / "fold_usability.csv")

    leak = {}
    for scheme, part in manifest.groupby("scheme"):
        leak[scheme] = audit_manifest(corrected, part, MAX_GAP)
    results["leakage_failures"] = {k: v for k, v in leak.items() if k != "time_2024_not_region_safe"}
    results["time_scheme_known_crossings"] = leak["time_2024_not_region_safe"]
    results["usable_folds"] = {k: {"usable": int(v["usable"].sum()), "total": len(v)} for k, v in usability.groupby("scheme")}

    # Near duplicates, all of which must be inside one cluster.
    protein = corrected["name"].str.extract(r"\(p\.([^)]+)\)")[0]
    same_protein = corrected[protein.notna() & corrected.assign(_p=protein).duplicated(["gene", "_p"], keep=False)]
    same_position = corrected[corrected.duplicated(["chrom", "pos"], keep=False)]
    residue_part = corrected[corrected["residue"] != ""]
    results["near_duplicates"] = {
        "unique_ids": int(corrected["variation_id"].nunique()), "duplicate_alleles": int(corrected.duplicated(["chrom", "pos", "ref", "alt"]).sum()), "duplicate_spdi": int(corrected["canonical_spdi"].duplicated().sum()),
        "rows_sharing_a_protein_change": len(same_protein), "rows_sharing_a_position": len(same_position), "positions_with_both_classes": int((same_position.groupby(["chrom", "pos"])["label"].nunique() > 1).sum()),
        "residues_with_both_classes": int((residue_part.groupby(["gene", "residue"])["label"].nunique() > 1).sum()),
        "protein_changes_in_one_cluster": bool((same_protein.assign(_p=protein).groupby(["gene", "_p"])["cluster"].nunique() == 1).all()),
        "positions_in_one_cluster": bool((same_position.groupby(["chrom", "pos"])["cluster"].nunique() == 1).all()),
    }
    cluster_stats = corrected.groupby("cluster").agg(gene=("gene", "first"), n=("variation_id", "size"), pathogenic=("label", lambda s: int((s == "pathogenic").sum())), benign=("label", lambda s: int((s == "benign").sum())))
    write_new(cluster_stats.reset_index(), output / "cluster_summary.csv")
    results["clusters"] = {"total": len(cluster_stats), "largest": int(cluster_stats["n"].max()), "with_pathogenic": int((cluster_stats["pathogenic"] > 0).sum()), "with_benign": int((cluster_stats["benign"] > 0).sum()), "mixed": int(((cluster_stats["pathogenic"] > 0) & (cluster_stats["benign"] > 0)).sum()),
                           "by_gene": {g: {"clusters": int((cluster_stats["gene"] == g).sum()), "with_pathogenic": int(((cluster_stats["gene"] == g) & (cluster_stats["pathogenic"] > 0)).sum()), "with_benign": int(((cluster_stats["gene"] == g) & (cluster_stats["benign"] > 0)).sum())} for g in GENES},
                           "largest_pathogenic_share": round(float(cluster_stats["pathogenic"].max() / cluster_stats["pathogenic"].sum()), 4),
                           "missense_clusters": int(missense["cluster"].nunique()), "missense_with_pathogenic": int(missense[missense["label"] == "pathogenic"]["cluster"].nunique()), "missense_with_benign": int(missense[missense["label"] == "benign"]["cluster"].nunique())}

    # Features.
    results["forbidden_feature_columns"] = {"main": forbidden_feature_columns(main_table.columns), "missense": forbidden_feature_columns(missense_table.columns)}
    feature_ids_ordered = main_table["variation_id"].tolist() == corrected.sort_values(["gene"], kind="stable")["variation_id"].tolist()
    results["features_main_row_order_follows_gene_and_position"] = bool(feature_ids_ordered)
    context = corrected[["variation_id", "label", "consequence_group", "gene"]]
    audit = pd.concat([feature_audit(main_table, context, "features_main"), feature_audit(missense_table, context, "features_missense")], ignore_index=True)
    write_new(audit, output / "FEATURE_AUDIT.csv")

    results["counts"] = {
        "variants": len(corrected), "by_label": corrected["label"].value_counts().to_dict(), "by_gene_label": {f"{g}|{l}": int(n) for (g, l), n in corrected.groupby(["gene", "label"]).size().items()},
        "two_star_by_gene_label": {f"{g}|{l}": int(n) for (g, l), n in corrected[corrected["stars"] == "2"].groupby(["gene", "label"]).size().items()},
        "missense_by_gene_label": {f"{g}|{l}": int(n) for (g, l), n in missense.groupby(["gene", "label"]).size().items()}, "missense": len(missense),
        "variant_types": main_table["variant_type"].value_counts().to_dict(), "empty_fields": {c: int((corrected[c] == "").sum()) for c in corrected.columns if (corrected[c] == "").any()},
        "rows_noted": corrected.loc[corrected["data_note"] != "", ["vcv", "name", "consequence_group"]].to_dict("records"),
    }
    results["environment"] = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scikit_learn": sklearn.__version__}
    results["parameters"] = {"max_gap": MAX_GAP, "seeds": list(SEEDS), "n_folds": N_FOLDS, "stratify": "gene and label", "proxy_screen_threshold": PROXY_THRESHOLD}
    inputs = [FREEZE / "task1_labels.csv", FREEZE / "splits" / "time_2024_all.csv", FEATURES / "features_main.csv", FEATURES / "features_missense.csv", FEATURES / "baselines.csv", FEATURES / "excluded_variants.csv", AUDIT_DERIVED / "record_submitter_meta.csv.gz"]
    results["inputs_sha256"] = {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in inputs}
    results["outputs_sha256"] = {path.name: sha256(path) for path in sorted(output.glob("*.csv"))}
    results["models_fitted"] = 0
    with (output / "validation_results.json").open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(results, handle, indent=1, default=str)
    print(json.dumps({key: results[key] for key in ("counts", "clusters", "near_duplicates", "join_problems", "usable_folds", "forbidden_feature_columns", "feature_build_used_this_freeze", "features_main_row_order_follows_gene_and_position")}, indent=0, default=str))
    print("leakage failures:", {k: len(v) for k, v in results["leakage_failures"].items()}, "time crossings folds:", len(results["time_scheme_known_crossings"]), results["time_scheme_known_crossings"][:2])
    print(audit.drop(columns=["table", "n"]).to_string())


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True, help="A new folder; existing folders are refused.")
    args = parser.parse_args(argv)
    build(args.output)


if __name__ == "__main__":
    main()
