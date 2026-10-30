"""Truth-backed correctness evaluation of PanAligner and minigraph on simulated reads.

Maps the reads written by simulate_truth_reads.py against each gene's combined train
graph and against the linear reference, then scores them two ways:

  mapeval    `paftools.js mapeval` on the GAF, as in PanAligner/data/evaluate.sh. With
             stable coordinates an alignment counts as correct only if its target is the
             reference name and it overlaps the true interval (overlap/union > 0.1).
             An alignment whose path passes through a non-reference node is counted as
             wrong, whatever its position.
  projected  the same rule, but such paths are first projected back onto the reference
             when both path ends lie on reference segments with one orientation;
             paths that cannot be projected are counted as `unprojectable`, not wrong.

Unaligned reads are simulated reads with no GAF line; mapeval does not count them.
Linux only; on the laptop run it inside WSL.
"""

from __future__ import annotations

import argparse
import csv
import platform
import re
import statistics
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from common import ROOT, ensure_dir, project_relative, read_json, write_json
from simulate_truth_reads import K8_BINARY, PAFTOOLS_JS, TRUTH_DIR, parse_truth_name, sha256_file


ALIGNERS = {
    "panaligner": ROOT / "PanAligner" / "PanAligner",
    "minigraph": ROOT / "minigraph" / "minigraph",
}
MAPQ_THRESHOLDS = (0, 1, 5, 10, 20, 30, 40, 50, 60)
OVERLAP_RATIO = 0.1
STABLE_SEGMENT = re.compile(r"([<>])(.+?):(\d+)-(\d+)(?=[<>]|$)")


def run_aligner(binary: Path, graph: Path, reads: Path, output_gaf: Path, threads: int) -> float:
    ensure_dir(output_gaf.parent)
    start = time.perf_counter()
    with output_gaf.open("w", encoding="utf-8") as out, Path(f"{output_gaf}.log").open("w", encoding="utf-8") as log:
        subprocess.run([str(binary), "-t", str(threads), "-cx", "lr", str(graph), str(reads)], stdout=out, stderr=log, check=True)
    return time.perf_counter() - start


def run_mapeval(gaf: Path) -> str:
    completed = subprocess.run(
        [str(K8_BINARY), str(PAFTOOLS_JS), "mapeval", str(gaf)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True,
    )
    return completed.stdout


def parse_mapeval(text: str) -> list[dict]:
    """Rows of `Q <mapq> <bin reads> <bin wrong> <cumulative error rate> <cumulative reads>`."""
    rows = []
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "Q":
            rows.append({
                "mapq": int(fields[1]), "reads": int(fields[2]), "wrong": int(fields[3]),
                "cumulative_error_rate": float(fields[4]), "cumulative_reads": int(fields[5]),
            })
    return rows


def overlaps(truth: tuple[int, int], start: int, end: int) -> bool:
    """paftools mapeval's rule: overlap / union of the two intervals > OVERLAP_RATIO."""
    overlap = min(truth[1], end) - max(truth[0], start)
    if overlap <= 0:
        return False
    union = max(truth[1], end) - min(truth[0], start)
    return overlap / union > OVERLAP_RATIO


def project_stable_path(path: str, path_length: int, path_start: int, path_end: int, reference: str) -> tuple[int, int, bool] | None:
    """Project a stable-coordinate GAF path onto the reference.

    Only the two end segments matter: both must be reference segments traversed in
    the same orientation (inner nodes, e.g. 1 bp inversion bubbles, are ignored).
    Returns (start, end, reversed), or None when that does not hold.
    """
    segments = STABLE_SEGMENT.findall(path)
    if not segments or "".join(f"{o}{n}:{s}-{e}" for o, n, s, e in segments) != path:
        return None
    first, last = segments[0], segments[-1]
    if first[1] != reference or last[1] != reference or first[0] != last[0]:
        return None
    tail = path_length - path_end
    if first[0] == ">":
        return int(first[2]) + path_start, int(last[3]) - tail, False
    return int(last[2]) + tail, int(first[3]) - path_start, True


def load_alignments(gaf: Path) -> dict[str, list[list[str]]]:
    by_read: dict[str, list[list[str]]] = defaultdict(list)
    with gaf.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[4] not in ("+", "-"):
                continue
            tags = fields[12:]
            # mapeval skips minimap2-style secondaries: an s1 tag without an s2 tag.
            if any(tag.startswith("s1:i:") for tag in tags) and not any(tag.startswith("s2:i:") for tag in tags):
                continue
            by_read[fields[0]].append(fields)
    return by_read


def classify_reads(read_names: list[str], gaf: Path, reference: str) -> list[dict]:
    """Per-read outcome for the longest alignment (mapeval mode 0: most matching bases, first wins)."""
    alignments = load_alignments(gaf)
    outcomes = []
    for name in read_names:
        _, target, true_start, true_end, true_strand = parse_truth_name(name)
        candidates = alignments.get(name, [])
        if not candidates:
            outcomes.append({"read": name, "status": "unaligned", "mapq": None, "mapeval_correct": False,
                             "projected_correct": False, "off_reference_path": False, "unprojectable": False})
            continue
        best = candidates[0]
        for fields in candidates[1:]:
            if int(fields[9]) > int(best[9]):
                best = fields
        strand, path = best[4], best[5]
        path_length, path_start, path_end = int(best[6]), int(best[7]), int(best[8])
        mapq = min(int(best[11]), 60)
        truth = (true_start, true_end)
        on_reference = path == target
        mapeval_correct = on_reference and strand == true_strand and overlaps(truth, path_start, path_end)
        projected_correct = mapeval_correct
        unprojectable = False
        if not on_reference:
            projection = project_stable_path(path, path_length, path_start, path_end, reference)
            unprojectable = projection is None
            if projection is not None:
                start, end, flipped = projection
                projected_strand = {"+": "-", "-": "+"}[strand] if flipped else strand
                projected_correct = projected_strand == true_strand and overlaps(truth, start, end)
        outcomes.append({
            "read": name, "status": "aligned", "mapq": mapq,
            "mapeval_correct": mapeval_correct, "projected_correct": projected_correct,
            "off_reference_path": not on_reference, "unprojectable": unprojectable,
        })
    return outcomes


def summarize(outcomes: list[dict]) -> dict:
    total = len(outcomes)
    aligned = [o for o in outcomes if o["status"] == "aligned"]
    summary = {
        "reads": total,
        "aligned": len(aligned),
        "unaligned": total - len(aligned),
        "mapeval_correct": sum(o["mapeval_correct"] for o in aligned),
        "mapeval_wrong": sum(not o["mapeval_correct"] for o in aligned),
        "projected_correct": sum(o["projected_correct"] for o in aligned),
        "projected_wrong": sum(not o["projected_correct"] and not o["unprojectable"] for o in aligned),
        "unprojectable": sum(o["unprojectable"] for o in aligned),
        "off_reference_path": sum(o["off_reference_path"] for o in aligned),
    }
    for key in ("aligned", "unaligned", "mapeval_correct", "mapeval_wrong", "projected_correct", "projected_wrong"):
        summary[f"{key}_pct"] = 100.0 * summary[key] / total if total else 0.0
    return summary


def mapq_curve(outcomes: list[dict]) -> list[dict]:
    """Fig. 6-style points: reads aligned at MAPQ >= Q and the share of them that are wrong."""
    total = len(outcomes)
    rows = []
    for threshold in MAPQ_THRESHOLDS:
        kept = [o for o in outcomes if o["status"] == "aligned" and o["mapq"] >= threshold]
        wrong = sum(not o["mapeval_correct"] for o in kept)
        rows.append({
            "mapq_threshold": threshold, "aligned": len(kept), "wrong": wrong,
            "aligned_fraction": len(kept) / total if total else 0.0,
            "wrong_fraction_of_aligned": wrong / len(kept) if kept else 0.0,
        })
    return rows


def read_fasta_names(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line[1:].strip() for line in handle if line.startswith(">")]


def tool_version(binary: Path) -> str:
    completed = subprocess.run([str(binary), "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    return completed.stdout.strip()


def git_state() -> dict:
    def git(*args: str) -> str:
        completed = subprocess.run(["git", *args], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=False)
        return completed.stdout.strip()
    return {"commit": git("rev-parse", "HEAD"), "dirty_paths": len(git("status", "--porcelain").splitlines())}


def plot_curves(curves: dict[tuple, list[dict]], output_png: Path, title: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    for (aligner, graph), rows in sorted(curves.items()):
        axis.plot(
            [row["wrong_fraction_of_aligned"] for row in rows],
            [row["aligned_fraction"] for row in rows],
            marker="o", markersize=3, label=f"{aligner} / {graph}",
        )
    axis.set_xlabel("Fraction of aligned reads that are wrong (mapeval)")
    axis.set_ylabel("Fraction of reads aligned at MAPQ >= Q")
    axis.set_title(title)
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output_png, dpi=150)
    plt.close(figure)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Score PanAligner and minigraph on truth-labelled simulated reads.")
    parser.add_argument("--source", choices=("reference", "heldout"), default="reference")
    parser.add_argument("--graphs", nargs="+", choices=("combined", "linear"), default=["combined", "linear"])
    parser.add_argument("--aligners", nargs="+", choices=tuple(ALIGNERS), default=list(ALIGNERS))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, default=TRUTH_DIR)
    args = parser.parse_args(argv)
    output_dir = args.output_dir.resolve()

    simulation = read_json(output_dir / "reads" / args.source / "simulation_manifest.json")
    result_dir = ensure_dir(output_dir / "results" / args.source)
    rows: list[dict] = []
    curves: dict[str, dict[tuple, list[list[dict]]]] = defaultdict(lambda: defaultdict(list))
    inputs: dict[str, str] = {}

    for run in simulation["runs"]:
        gene, seed = run["gene"], run["seed"]
        reads = ROOT / run["reads_fasta"]
        read_names = read_fasta_names(reads)
        graph_paths = {
            "combined": ROOT / "graphs" / "train" / f"{gene.lower()}.combined.gfa",
            "linear": ROOT / run["reference_fasta"],
        }
        for graph in args.graphs:
            inputs[project_relative(graph_paths[graph])] = sha256_file(graph_paths[graph])
            for aligner in args.aligners:
                run_dir = ensure_dir(result_dir / gene / f"seed{seed}")
                gaf = run_dir / f"{aligner}.{graph}.gaf"
                seconds = run_aligner(ALIGNERS[aligner], graph_paths[graph], reads, gaf, args.threads)
                mapeval_text = run_mapeval(gaf)
                (run_dir / f"{aligner}.{graph}.mapeval.txt").write_text(mapeval_text, encoding="utf-8")
                mapeval_rows = parse_mapeval(mapeval_text)
                outcomes = classify_reads(read_names, gaf, run["reference_name"])
                summary = summarize(outcomes)
                mapeval_wrong = sum(row["wrong"] for row in mapeval_rows)
                mapeval_reads = sum(row["reads"] for row in mapeval_rows)
                if (mapeval_reads, mapeval_wrong) != (summary["aligned"], summary["mapeval_wrong"]):
                    raise RuntimeError(
                        f"{gene} seed{seed} {aligner}/{graph}: own strict count {summary['aligned']}/{summary['mapeval_wrong']} "
                        f"disagrees with mapeval {mapeval_reads}/{mapeval_wrong}"
                    )
                curve = mapq_curve(outcomes)
                curves[gene][(aligner, graph)].append(curve)
                rows.append({"source": args.source, "gene": gene, "seed": seed, "aligner": aligner, "graph": graph,
                             **summary, "wall_seconds": round(seconds, 3)})
                print(f"{args.source} {gene} seed{seed} {aligner}/{graph}: {summary['aligned']}/{summary['reads']} aligned, "
                      f"mapeval wrong {summary['mapeval_wrong']}, projected wrong {summary['projected_wrong']}, "
                      f"unprojectable {summary['unprojectable']}, off-reference paths {summary['off_reference_path']}", flush=True)

    with (result_dir / "per_run.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["gene"], row["aligner"], row["graph"])].append(row)
    condition_rows = []
    for (gene, aligner, graph), members in sorted(grouped.items()):
        condition = {"source": args.source, "gene": gene, "aligner": aligner, "graph": graph, "seeds": len(members),
                     "reads_total": sum(m["reads"] for m in members)}
        for key in ("aligned_pct", "unaligned_pct", "mapeval_wrong_pct", "projected_wrong_pct"):
            values = [m[key] for m in members]
            condition[f"{key}_mean"] = round(statistics.mean(values), 4)
            condition[f"{key}_sd"] = round(statistics.stdev(values), 4) if len(values) > 1 else 0.0
        condition["off_reference_path_total"] = sum(m["off_reference_path"] for m in members)
        condition["unprojectable_total"] = sum(m["unprojectable"] for m in members)
        condition_rows.append(condition)
    with (result_dir / "by_condition.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(condition_rows[0]))
        writer.writeheader()
        writer.writerows(condition_rows)

    curve_rows = []
    for gene, by_condition in curves.items():
        pooled: dict[tuple, list[dict]] = {}
        for condition, per_seed in by_condition.items():
            pooled[condition] = []
            for index, threshold in enumerate(MAPQ_THRESHOLDS):
                aligned = sum(curve[index]["aligned"] for curve in per_seed)
                wrong = sum(curve[index]["wrong"] for curve in per_seed)
                reads_total = sum(r["reads"] for r in rows if r["gene"] == gene and (r["aligner"], r["graph"]) == condition)
                point = {"gene": gene, "aligner": condition[0], "graph": condition[1], "mapq_threshold": threshold,
                         "aligned": aligned, "wrong": wrong,
                         "aligned_fraction": aligned / reads_total if reads_total else 0.0,
                         "wrong_fraction_of_aligned": wrong / aligned if aligned else 0.0}
                pooled[condition].append(point)
                curve_rows.append(point)
        plot_curves(pooled, result_dir / f"{gene.lower()}.mapq_curve.png", f"{gene} ({args.source} reads, seeds pooled)")
    with (result_dir / "mapq_curve.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(curve_rows[0]))
        writer.writeheader()
        writer.writerows(curve_rows)

    write_json(result_dir / "run_manifest.json", {
        "source": args.source,
        "git": git_state(),
        "host": platform.platform(),
        "threads": args.threads,
        "aligner_command": "<aligner> -t <threads> -cx lr <graph> <reads>",
        "mapeval_command": "k8 paftools.js mapeval <gaf>  (defaults: -r 0.1, -m 0 longest alignment)",
        "tool_versions": {
            **{name: tool_version(binary) for name, binary in ALIGNERS.items() if name in args.aligners},
            "paftools": subprocess.run([str(K8_BINARY), str(PAFTOOLS_JS), "version"], stdout=subprocess.PIPE, text=True, check=False).stdout.strip(),
        },
        "simulation_manifest": project_relative(output_dir / "reads" / args.source / "simulation_manifest.json"),
        "input_sha256": inputs,
        "notes": [
            "Locus-scale benchmark on APP/PSEN1/PSEN2 regions; not the paper's whole-genome 10H-95H benchmark.",
            "wall_seconds is laptop/WSL wall time for one aligner call; it is not a runtime benchmark.",
            "mapeval counts alignments through non-reference nodes as wrong; projected_* re-scores those paths.",
        ],
    })


if __name__ == "__main__":
    main()
