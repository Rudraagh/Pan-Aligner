"""Descriptive confound checks for the Model B labels. Fits and retrains no model; reads existing tables."""
from __future__ import annotations

import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = HERE.parents[2]


def main() -> None:
    table = list(csv.DictReader((ROOT / "outputs" / "ml_variant" / "variant_table.csv").open(encoding="utf-8")))
    prov = list(csv.DictReader((OUT / "variant_provenance.csv").open(encoding="utf-8")))
    result: dict = {}
    snv = {"HEALTHY": 0, "UNHEALTHY": 0}
    for r in prov:
        snv[r["original_label"]] += r["variant_type"] == "SNV"
    result["snv_by_label"] = snv
    result["snv_vs_other_fisher_p"] = fisher_exact([[snv["HEALTHY"], 110 - snv["HEALTHY"]], [snv["UNHEALTHY"], 110 - snv["UNHEALTHY"]]])[1]
    protein = {"missense_variant", "frameshift_variant", "start_lost", "inframe_insertion", "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant"}
    counts = {label: sum(r["original_label"] == label and r["ensembl_consequence"] in protein for r in prov) for label in snv}
    known = {label: sum(r["original_label"] == label and not r["ensembl_consequence"].startswith("not queried") for r in prov) for label in snv}
    result["protein_or_splice_altering"] = {"counts": counts, "with_consequence": known, "fisher_p": fisher_exact([[counts["HEALTHY"], known["HEALTHY"] - counts["HEALTHY"]], [counts["UNHEALTHY"], known["UNHEALTHY"] - counts["UNHEALTHY"]]])[1]}
    afs = {label: [float(r["gnomad_r4_af"]) for r in prov if r["original_label"] == label and r["gnomad_r4_af"]] for label in snv}
    result["gnomad"] = {label: {"present": len(v), "absent_or_unqueried": 110 - len(v), "median": statistics.median(v), "af_ge_1pct": sum(a >= 0.01 for a in v), "af_ge_0.1pct": sum(a >= 0.001 for a in v)} for label, v in afs.items()}

    by_gene = defaultdict(list)
    for r in table:
        by_gene[r["gene"]].append(r)
    neighbour, shared = {}, {}
    test_distance = Counter()
    test_total = same_label_nearest_test = 0
    for gene, rows in by_gene.items():
        same = 0
        distances = []
        positions = defaultdict(list)
        for r in rows:
            positions[int(r["position"])].append(r["label"])
            others = [o for o in rows if o is not r]
            nearest = min(others, key=lambda o: abs(int(o["position"]) - int(r["position"])))
            distances.append(abs(int(nearest["position"]) - int(r["position"])))
            same += nearest["label"] == r["label"]
            if r["split"] == "test":
                train = [o for o in rows if o["split"] == "train"]
                near = min(train, key=lambda o: abs(int(o["position"]) - int(r["position"])))
                gap = abs(int(near["position"]) - int(r["position"]))
                test_total += 1
                same_label_nearest_test += near["label"] == r["label"]
                for limit in (10, 100, 1000):
                    test_distance[f"within_{limit}_bp"] += gap <= limit
        neighbour[gene] = {"n": len(rows), "nearest_same_label": same, "fraction": round(same / len(rows), 3), "median_gap_bp": statistics.median(distances)}
        shared[gene] = {"positions_with_multiple_records": sum(len(v) > 1 for v in positions.values()), "positions_with_both_labels": sum(len(set(v)) > 1 for v in positions.values())}
        result.setdefault("position_range", {})[gene] = {label: (min(int(r["position"]) for r in rows if r["label"] == label), statistics.median(int(r["position"]) for r in rows if r["label"] == label), max(int(r["position"]) for r in rows if r["label"] == label)) for label in {r["label"] for r in rows}}
    result["nearest_neighbour_all_variants"] = neighbour
    result["shared_positions"] = shared
    result["seed42_split"] = {"split_counts": Counter(r["split"] for r in table), "test": test_total, "nearest_train_same_label": same_label_nearest_test, **test_distance}
    (OUT / "raw" / "confound_summary.json").write_text(json.dumps(result, indent=1, default=dict), encoding="utf-8")
    print(json.dumps(result, default=dict))


if __name__ == "__main__":
    main()
