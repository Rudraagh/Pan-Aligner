"""Compute the pilot's go/no-go gates, fixed before the run in outputs/model_a_pilot_plan/.

A gate with no evidence supplied is reported as NOT_EVALUATED and counts as not passed, never as a pass.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

MIN_INCORRECT = 500
MIN_WINDOWS = 30
WINDOW = 100_000
MIN_PER_ALIGNER = 100
MAX_KMER_DISAGREE = 0.05
MAX_TRUTH_UNPROJECTABLE = 0.10
MAX_PATH_UNPROJECTABLE = 0.20
EXPECTED_SETTINGS = 6
EXPECTED_ALIGNERS = 2


def compute(rows, kmer_verdicts: Counter | None, runs_complete: bool | None, reproducible: bool | None) -> dict:
    states: Counter = Counter()
    by_group: dict[tuple[str, str], Counter] = defaultdict(Counter)
    per_aligner: Counter = Counter()
    windows: set[tuple[str, int]] = set()
    kinds: Counter = Counter()
    errors: Counter = Counter()
    combos = set()
    no_truth = 0
    for row in rows:
        state, aligner = row["truth_state"], row["aligner"]
        no_truth += row["truth_reason"] != "ok"
        combos.add((row["setting"], aligner))
        states[state] += 1
        by_group[(row["setting"], aligner)][state] += 1
        if row["path_projected"]:
            kinds[row["path_projected"]] += 1
        if state == "incorrect":
            per_aligner[aligner] += 1
            windows.add((row["truth_ref_chrom"], int(row["truth_ref_start"]) // WINDOW))
            errors[row["error_type"]] += 1
    total = sum(states.values())
    aligned = total - states["unaligned"]
    judged = states["correct"] + states["incorrect"]

    gates: dict[str, dict] = {}

    def gate(name: str, passed: bool | None, value, rule: str) -> None:
        gates[name] = {"result": "NOT_EVALUATED" if passed is None else "PASS" if passed else "FAIL", "value": value, "rule": rule}

    all_labelled = len(combos) == EXPECTED_SETTINGS * EXPECTED_ALIGNERS
    gate("runs", None if runs_complete is None else runs_complete and all_labelled, len(combos), f"{EXPECTED_SETTINGS} settings x {EXPECTED_ALIGNERS} aligners labelled, and jobs finished")
    gate("errors", states["incorrect"] >= MIN_INCORRECT, states["incorrect"], f">= {MIN_INCORRECT} incorrect rows")
    gate("spread", len(windows) >= MIN_WINDOWS, len(windows), f">= {MIN_WINDOWS} separate {WINDOW} bp windows")
    weakest = min((per_aligner[a] for a in {c[1] for c in combos}), default=0)
    gate("per_aligner", weakest >= MIN_PER_ALIGNER, dict(per_aligner), f">= {MIN_PER_ALIGNER} incorrect rows per aligner")
    if kmer_verdicts is None:
        gate("independence", None, None, f"< {MAX_KMER_DISAGREE:.0%} of assessable reads disagree")
    else:
        assessable = kmer_verdicts["agree"] + kmer_verdicts["disagree"]
        share = kmer_verdicts["disagree"] / assessable if assessable else None
        gate("independence", None if share is None else share < MAX_KMER_DISAGREE, {"share": share, **kmer_verdicts}, f"< {MAX_KMER_DISAGREE:.0%} of assessable reads disagree")
    truth_share = no_truth / total if total else None
    gate("truth_coverage", None if truth_share is None else truth_share < MAX_TRUTH_UNPROJECTABLE, truth_share, f"source interval unprojectable for < {MAX_TRUTH_UNPROJECTABLE:.0%} of rows")
    path_share = states["path_unprojectable"] / aligned if aligned else None
    gate("path_coverage", None if path_share is None else path_share < MAX_PATH_UNPROJECTABLE, path_share, f"path_unprojectable < {MAX_PATH_UNPROJECTABLE:.0%} of aligned rows")
    gate("reproducibility", reproducible, reproducible, "same seed gives identical reads")

    return {
        "verdict": "GO" if all(g["result"] == "PASS" for name, g in gates.items() if name != "per_aligner") else "NO_GO",
        "per_aligner_note": "A per_aligner failure does not stop production by itself, but limits comparison between aligners.",
        "gates": gates,
        "states": dict(states),
        "judged_rows": judged,
        "path_projection_kinds": dict(kinds),
        "error_types": dict(errors),
        "by_setting_aligner": {f"{s}|{a}": dict(c) for (s, a), c in sorted(by_group.items())},
        "scope": "chr21 pilot on simulated reads. Says nothing about model performance, real reads or other chromosomes.",
    }


def read_rows(paths: list[Path]):
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            yield from csv.DictReader(handle)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, nargs="+", required=True)
    parser.add_argument("--kmer-check", type=Path)
    parser.add_argument("--runs-complete", choices=("yes", "no"))
    parser.add_argument("--repro-json", type=Path, help='JSON with {"identical": true|false}.')
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    verdicts = None
    if args.kmer_check:
        with gzip.open(args.kmer_check, "rt", encoding="utf-8", newline="") as handle:
            verdicts = Counter(row["verdict"] for row in csv.DictReader(handle))
    reproducible = json.loads(args.repro_json.read_text(encoding="utf-8"))["identical"] if args.repro_json else None
    runs = None if args.runs_complete is None else args.runs_complete == "yes"
    report = compute(read_rows(args.features), verdicts, runs, reproducible)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    print(json.dumps({"verdict": report["verdict"], "gates": {k: v["result"] for k, v in report["gates"].items()}}, indent=2))


if __name__ == "__main__":
    main()
