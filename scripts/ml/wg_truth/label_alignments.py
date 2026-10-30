"""Label one GAF against simulator truth and write one row per read: its longest alignment, or `unaligned`.

Truth states, in the order they are decided:
  unaligned            no alignment
  truth_unprojectable  the source interval has no clean reference interval
  truth_ambiguous      at least AMBIGUOUS_FRACTION of the source lies in low-mappability or segdup regions
  path_unprojectable   the alignment has no reference interval, so it cannot be judged
  correct / incorrect  overlap/union above the threshold, same chromosome and strand

truth_state uses overlap 0.1 (paftools mapeval); truth_state_strict repeats it at 0.5.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from ml.build_alignment_dataset import FEATURE_COLUMNS, alignment_features  # noqa: E402
from ml.wg_truth.project_paths import PathProjection, project_path  # noqa: E402
from ml.wg_truth.project_truth import RefInterval, load_alignments, project_interval  # noqa: E402
from simulate_truth_reads import parse_truth_name  # noqa: E402
from truth_eval import load_alignments as load_gaf  # noqa: E402

AMBIGUOUS_FRACTION = 0.9
OVERLAP_DEFAULT = 0.1
OVERLAP_STRICT = 0.5
ROW_COLUMNS = [
    "read", "setting", "accuracy", "read_length", "seed", "aligner", "haplotype",
    "source_contig", "source_start", "source_end", "source_strand",
    "truth_ref_chrom", "truth_ref_start", "truth_ref_end", "truth_ref_strand", "truth_reason",
    "ambiguous_fraction", "difficult_fraction",
    "align_chrom", "align_start", "align_end", "align_strand", "path_projected",
    "truth_state", "truth_state_strict", "error_type", *FEATURE_COLUMNS,
]


class BedCover:
    """Merged intervals of one BED, by contig, answering what fraction of an interval they cover."""

    def __init__(self, paths: Path | list[Path] | None):
        self.contigs: dict[str, tuple[list[int], list[int], list[int]]] = {}
        if paths is None:
            return
        raw: dict[str, list[tuple[int, int]]] = {}
        for path in [paths] if isinstance(paths, Path) else paths:
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(("#", "track", "browser")) or not line.strip():
                        continue
                    fields = line.split("\t")
                    raw.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
        for contig, intervals in raw.items():
            merged: list[list[int]] = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            starts, ends, cumulative, total = [], [], [], 0
            for start, end in merged:
                starts.append(start)
                ends.append(end)
                cumulative.append(total)
                total += end - start
            self.contigs[contig] = (starts, ends, cumulative)

    def covered_to(self, contig: str, position: int) -> int:
        starts, ends, cumulative = self.contigs[contig]
        index = bisect.bisect_right(starts, position) - 1
        if index < 0:
            return 0
        return cumulative[index] + min(position, ends[index]) - starts[index]

    def fraction(self, contig: str, start: int, end: int) -> float:
        if contig not in self.contigs or end <= start:
            return 0.0
        return (self.covered_to(contig, end) - self.covered_to(contig, start)) / (end - start)


def overlap_ratio(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = min(a_end, b_end) - max(a_start, b_start)
    if overlap <= 0:
        return 0.0
    return overlap / (max(a_end, b_end) - min(a_start, b_start))


def judge(truth: RefInterval, path: PathProjection, threshold: float) -> tuple[str, str]:
    """(state, error type) for an alignment with both intervals known."""
    if path.chrom != truth.chrom:
        return "incorrect", "wrong_chromosome"
    if path.strand != truth.strand:
        return "incorrect", "wrong_strand"
    if overlap_ratio(truth.start, truth.end, path.start, path.end) > threshold:
        return "correct", ""
    return "incorrect", "same_chromosome_wrong_place"


def label_read(name: str, candidates: list[list[str]], alignments: dict, ambiguous: BedCover, difficult: BedCover, reference_names: set[str]) -> dict:
    _, contig, source_start, source_end, source_strand = parse_truth_name(name)
    row = {column: "" for column in ROW_COLUMNS}
    row.update(read=name, source_contig=contig, source_start=source_start, source_end=source_end, source_strand=source_strand)
    truth, reason = project_interval(alignments.get(contig, []), source_start, source_end, source_strand)
    row["truth_reason"] = reason
    if truth is not None:
        row.update(truth_ref_chrom=truth.chrom, truth_ref_start=truth.start, truth_ref_end=truth.end, truth_ref_strand=truth.strand)
    ambiguous_fraction = ambiguous.fraction(contig, source_start, source_end)
    row["ambiguous_fraction"] = round(ambiguous_fraction, 4)
    row["difficult_fraction"] = round(difficult.fraction(contig, source_start, source_end), 4)
    if not candidates:
        row["truth_state"] = row["truth_state_strict"] = "unaligned"
        return row
    best = candidates[0]
    for fields in candidates[1:]:
        if int(fields[9]) > int(best[9]):
            best = fields
    row.update(alignment_features(best, len(candidates)))
    path = project_path(best[4], best[5], int(best[7]), int(best[8]), reference_names)
    if path is not None:
        row.update(align_chrom=path.chrom, align_start=path.start, align_end=path.end, align_strand=path.strand, path_projected=path.kind)
    if truth is None:
        state = strict = "truth_unprojectable"
    elif ambiguous_fraction >= AMBIGUOUS_FRACTION:
        state = strict = "truth_ambiguous"
    elif path is None:
        state = strict = "path_unprojectable"
    else:
        state, error = judge(truth, path, OVERLAP_DEFAULT)
        strict, _ = judge(truth, path, OVERLAP_STRICT)
        row["error_type"] = error
    row["truth_state"], row["truth_state_strict"] = state, strict
    return row


def read_names_from_fasta(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line[1:].strip() for line in handle if line.startswith(">")]


def require_contigs(cover: BedCover, contigs: set[str], label: str) -> None:
    """A BED that uses other contig names would silently report no ambiguity, so that is an error."""
    missing = contigs - set(cover.contigs)
    if cover.contigs and missing:
        raise SystemExit(f"{label} BED lacks source contigs {sorted(missing)}; it has e.g. {sorted(cover.contigs)[:5]}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--gaf", type=Path, required=True)
    parser.add_argument("--alignments", type=Path, nargs="+", required=True, help="Haplotype-to-reference PAF or SAM files.")
    parser.add_argument("--reference-names", type=Path, required=True, help=".fai of the graph's reference, or a list with names first.")
    parser.add_argument("--ambiguous-bed", type=Path, nargs="*", default=[])
    parser.add_argument("--difficult-bed", type=Path, nargs="*", default=[])
    parser.add_argument("--setting", required=True)
    parser.add_argument("--accuracy", type=float, required=True)
    parser.add_argument("--read-length", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; refusing to overwrite it")

    names = read_names_from_fasta(args.reads)
    contigs = {parse_truth_name(name)[1] for name in names}
    reference_names = {line.split("\t")[0].strip() for line in args.reference_names.read_text(encoding="utf-8").splitlines() if line.strip()}
    alignments = load_alignments(args.alignments)
    missing = contigs - set(alignments)
    if missing:
        raise SystemExit(f"no haplotype-to-reference alignment for {sorted(missing)}; alignments cover {sorted(alignments)[:5]}")
    ambiguous = BedCover(args.ambiguous_bed)
    difficult = BedCover(args.difficult_bed)
    require_contigs(ambiguous, contigs, "ambiguous")
    require_contigs(difficult, contigs, "difficult")
    gaf = load_gaf(args.gaf)
    fixed = {"setting": args.setting, "accuracy": args.accuracy, "read_length": args.read_length, "seed": args.seed, "aligner": args.aligner}
    counts: dict[str, int] = {}
    with gzip.open(args.output, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for name in names:
            row = label_read(name, gaf.get(name, []), alignments, ambiguous, difficult, reference_names)
            row.update(fixed, haplotype=name.split(".", 1)[0])
            writer.writerow(row)
            counts[row["truth_state"]] = counts.get(row["truth_state"], 0) + 1
    print(args.output.name, dict(sorted(counts.items())), flush=True)


if __name__ == "__main__":
    main()
