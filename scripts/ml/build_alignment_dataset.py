"""Model A dataset: per-alignment features labelled by simulation truth.

Reads are simulated with PBSIM2 from the held-out haplotypes (never in the graphs)
over a grid of read accuracies and lengths, then mapped by PanAligner and minigraph
to every gene's combined train graph:

  own gene graph   label from the true interval: `correct` if the longest alignment
                   overlaps it (paftools mapeval rule, overlap/union > 0.1, same
                   strand), with paths through non-reference nodes first projected
                   onto the reference; `incorrect` otherwise; `unprojectable` if the
                   path cannot be projected (kept out of training)
  other gene graph a decoy: the read comes from another chromosome, so any
                   alignment is `incorrect`

Labels come only from the simulator, never from the aligners being judged; strict
mapeval "wrong" (any path through a non-reference node) is not used, because those
are scoring artefacts. Reads with no alignment are kept as `unaligned` rows without
features. Linux only (PBSIM2, k8, PanAligner, minigraph); on the laptop run in WSL.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
import sys
import time
from argparse import Namespace
from collections import Counter
from pathlib import Path

# Runnable as a script: put scripts/ on the import path for common, truth_eval, ml.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, project_relative, write_json
from simulate_truth_reads import DEFAULT_MODEL, DEFAULT_PBSIM, GENES, parse_truth_name, read_fasta, sha256_file, simulate_gene
from truth_eval import ALIGNERS, STABLE_SEGMENT, git_state, load_alignments, overlaps, project_stable_path, read_fasta_names, run_aligner, tool_version


FEATURE_COLUMNS = [
    "mapq", "query_length", "query_coverage", "residue_matches", "block_length", "identity", "nm_rate",
    "chain_minimizers", "chain_score", "second_chain_score", "second_to_best_chain_ratio", "divergence",
    "path_segments", "off_reference_path", "alignment_count",
]
ROW_COLUMNS = ["setting", "accuracy", "read_length", "seed", "aligner", "own_gene", "graph_gene", "decoy", "group", "read", "label", *FEATURE_COLUMNS]
DEFAULT_OUTPUT = OUTPUTS_DIR / "ml_alignment"


def setting_name(accuracy: float, length: int) -> str:
    return f"acc{round(accuracy * 100)}_len{length}"


def read_group(read_name: str) -> str:
    """Source haplotype of a held-out read (`<haplotype>.<pbsim id>!...`), used for group-disjoint splits."""
    read_id = read_name.split("!", 1)[0]
    return read_id.rsplit(".", 1)[0]


def _tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def alignment_features(fields: list[str], alignment_count: int) -> dict:
    tags = _tags(fields)
    query_length, query_start, query_end = int(fields[1]), int(fields[2]), int(fields[3])
    matches, block = int(fields[9]), int(fields[10])
    chain_score, second_score = int(tags.get("s1", 0)), int(tags.get("s2", 0))
    segments = STABLE_SEGMENT.findall(fields[5])
    return {
        "mapq": min(int(fields[11]), 60),
        "query_length": query_length,
        "query_coverage": max(0, query_end - query_start) / query_length if query_length else 0.0,
        "residue_matches": matches,
        "block_length": block,
        "identity": matches / block if block else 0.0,
        "nm_rate": int(tags.get("NM", 0)) / block if block else 0.0,
        "chain_minimizers": int(tags.get("cm", 0)),
        "chain_score": chain_score,
        "second_chain_score": second_score,
        "second_to_best_chain_ratio": second_score / chain_score if chain_score else 0.0,
        "divergence": float(tags.get("dv", -1.0)),
        "path_segments": len(segments) if segments else 1,
        "off_reference_path": int(bool(segments)),
        "alignment_count": alignment_count,
    }


def label_alignment(fields: list[str], read_name: str, own_gene: str, graph_gene: str, reference_name: str) -> str:
    """`correct`, `incorrect` or `unprojectable` for one alignment of a simulated read."""
    if graph_gene != own_gene:
        return "incorrect"
    _, target, true_start, true_end, true_strand = parse_truth_name(read_name)
    strand, path = fields[4], fields[5]
    path_length, path_start, path_end = int(fields[6]), int(fields[7]), int(fields[8])
    truth = (true_start, true_end)
    if path == target:
        return "correct" if strand == true_strand and overlaps(truth, path_start, path_end) else "incorrect"
    projection = project_stable_path(path, path_length, path_start, path_end, reference_name)
    if projection is None:
        return "unprojectable"
    start, end, flipped = projection
    projected_strand = {"+": "-", "-": "+"}[strand] if flipped else strand
    return "correct" if projected_strand == true_strand and overlaps(truth, start, end) else "incorrect"


def rows_for_gaf(read_names: list[str], gaf: Path, *, setting: str, seed: int, aligner: str, own_gene: str,
                 graph_gene: str, reference_name: str, accuracy: float = 0.0, read_length: int = 0) -> list[dict]:
    """One row per simulated read: its longest alignment (mapeval mode 0, first wins), or `unaligned`."""
    alignments = load_alignments(gaf)
    rows = []
    for name in read_names:
        base = {"setting": setting, "accuracy": accuracy, "read_length": read_length, "seed": seed, "aligner": aligner,
                "own_gene": own_gene, "graph_gene": graph_gene, "decoy": int(graph_gene != own_gene),
                "group": read_group(name), "read": name}
        candidates = alignments.get(name, [])
        if not candidates:
            rows.append({**base, "label": "unaligned", **{column: "" for column in FEATURE_COLUMNS}})
            continue
        best = candidates[0]
        for fields in candidates[1:]:
            if int(fields[9]) > int(best[9]):
                best = fields
        rows.append({**base, "label": label_alignment(best, name, own_gene, graph_gene, reference_name),
                     **alignment_features(best, len(candidates))})
    return rows


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build the Model A alignment-reliability dataset from simulated reads.")
    parser.add_argument("--accuracies", nargs="+", type=float, default=[0.95, 0.85, 0.75])
    parser.add_argument("--lengths", nargs="+", type=int, default=[10000, 2000])
    parser.add_argument("--depth", type=float, default=2.0, help="Coverage per held-out haplotype.")
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--genes", nargs="+", choices=GENES, default=list(GENES))
    parser.add_argument("--aligners", nargs="+", choices=tuple(ALIGNERS), default=list(ALIGNERS))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--pbsim", type=Path, default=DEFAULT_PBSIM)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if sys.platform.startswith("win"):
        raise SystemExit("build_alignment_dataset.py needs Linux tools (PBSIM2, k8, PanAligner, minigraph); run it inside WSL.")

    output_dir = ensure_dir(args.output_dir.resolve())
    started = time.perf_counter()
    graphs = {gene: ROOT / "graphs" / "train" / f"{gene.lower()}.combined.gfa" for gene in args.genes}
    references = {gene: read_fasta(ROOT / "data" / "raw" / f"{gene.lower()}_reference.fa")[0][0] for gene in args.genes}
    rows: list[dict] = []
    simulations: list[dict] = []
    for accuracy in args.accuracies:
        for length in args.lengths:
            setting = setting_name(accuracy, length)
            sim_dir = output_dir / "sim" / setting
            sim_args = Namespace(depth=args.depth, length=length, accuracy=accuracy, pbsim=args.pbsim.resolve(),
                                 model=args.model.resolve(), source="heldout", test_manifest=args.test_manifest.resolve(),
                                 output_dir=sim_dir)
            for gene in args.genes:
                for seed in args.seeds:
                    record = simulate_gene(gene, seed, sim_args)
                    simulations.append({"setting": setting, "accuracy": accuracy, "read_length": length, **record})
                    reads = ROOT / record["reads_fasta"]
                    names = read_fasta_names(reads)
                    for graph_gene, graph in graphs.items():
                        for aligner in args.aligners:
                            gaf = sim_dir / "results" / gene / f"seed{seed}" / f"{aligner}.{graph_gene.lower()}.gaf"
                            run_aligner(ALIGNERS[aligner], graph, reads, gaf, args.threads)
                            rows.extend(rows_for_gaf(names, gaf, setting=setting, seed=seed, aligner=aligner, own_gene=gene,
                                                     graph_gene=graph_gene, reference_name=references[graph_gene],
                                                     accuracy=accuracy, read_length=length))
                    print(f"{setting} {gene} seed{seed}: {record['reads']} reads", flush=True)

    features_path = output_dir / "alignment_features.csv.gz"
    with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter((row["setting"], row["aligner"], "decoy" if row["decoy"] else "own", row["label"]) for row in rows)
    class_counts = [{"setting": s, "aligner": a, "graph": g, "label": l, "rows": n} for (s, a, g, l), n in sorted(counts.items())]
    write_json(output_dir / "class_counts.json", {"totals": dict(Counter(row["label"] for row in rows)), "by_condition": class_counts})
    write_json(output_dir / "dataset_manifest.json", {
        "command": "python3 " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "git": git_state(),
        "parameters": {"accuracies": args.accuracies, "lengths": args.lengths, "depth": args.depth, "seeds": args.seeds,
                       "genes": args.genes, "aligners": args.aligners, "source": "heldout haplotypes (test_manifest)"},
        "tool_versions": {name: tool_version(path) for name, path in ALIGNERS.items() if name in args.aligners},
        "graphs": {gene: {"path": project_relative(path), "sha256": sha256_file(path)} for gene, path in graphs.items()},
        "features": project_relative(features_path),
        "rows": len(rows),
        "simulations": simulations,
        "host": platform.platform(),
        "python": platform.python_version(),
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "labels": "simulation truth; own-gene alignments projected onto the reference; every decoy alignment is incorrect",
    })
    print(json.dumps({"rows": len(rows), "totals": dict(Counter(row["label"] for row in rows))}, indent=2))


if __name__ == "__main__":
    main()
