"""Write the per-variant comparison table and the per-submission table; print condition counts and input checksums."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location("tables", HERE / "04_build_tables.py")
tables = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tables)


def main() -> None:
    rows = list(csv.DictReader((OUT / "variant_provenance.csv").open(encoding="utf-8")))
    by_vcv = {r["vcv"]: r for r in tables.parse_vcv().values()}
    lines = ["# Provenance comparison table", "", "Original labels are copied unchanged. ClinVar retrieved 2026-10-30. Evidence strength is the ClinVar review status in stars.", "",
             "| Variant | Original label | External annotation | Condition | Evidence strength | Label agreement | Category and notes |", "|---|---|---|---|---|---|---|"]
    submissions = []
    condition_counts = Counter()
    for r in rows:
        vcv = r["clinvar_vcv"] or r["related_vcv"]
        name = r["clinvar_classification"] or (f"Related record only: {r['related_classification']}" if r["related_vcv"] else "None")
        stars = f"{r['clinvar_stars']} star" if r["clinvar_vcv"] else "Not applicable"
        condition = r["ad_rcv_classifications"] or r["related_conditions"] or "None"
        note = r["provenance_category"] + (f". {r['related_note']}" if r["related_note"] else "")
        lines.append(f"| {r['sequence_id']} | {r['original_label']} | {vcv} {name} | {condition} | {stars} | {r['label_agreement']} | {note} |".replace("  ", " "))
        record = by_vcv.get(r["clinvar_vcv"])
        if record:
            for s in record["submissions"]:
                submissions.append({"sequence_id": r["sequence_id"], "original_label": r["original_label"], "vcv": r["clinvar_vcv"], "scv": s["scv"], "submitter": s["submitter"], "submission_date": s["submitted"],
                                    "classification": s["classification"], "last_evaluated": s["last_evaluated"], "review_status": s["review"], "conditions": "; ".join(s["assertion_traits"])[:300]})
            for rcv in record["rcvs"]:
                for c in rcv["conditions"]:
                    condition_counts[(r["original_label"], c)] += 1
    (OUT / "provenance_comparison_table.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    with (OUT / "clinvar_submissions.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(submissions[0]))
        writer.writeheader(); writer.writerows(submissions)
    print("submissions", len(submissions))
    for label in ("UNHEALTHY", "HEALTHY"):
        print(label, [(c, n) for (l, c), n in condition_counts.most_common() if l == label][:12])
    for path in sorted((ROOT / "data" / "raw").glob("*_all_labeled.fa")) + sorted((ROOT / "data" / "raw").glob("*_reference.fa")) + sorted(ROOT.glob("*_combined.fasta")):
        print(path.relative_to(ROOT).as_posix(), hashlib.sha256(path.read_bytes()).hexdigest())
    print(json.dumps({"submitters": Counter(s["submitter"] for s in submissions).most_common(8)}))


if __name__ == "__main__":
    main()
