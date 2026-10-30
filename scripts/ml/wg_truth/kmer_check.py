"""Independence check on truth projection, using no aligner code.

For a random sample of reads, take the error-free source interval from the haplotype and find where its
31-mers occur in the reference. K-mers occurring exactly once in the reference vote for a location, which is
compared with the projected truth. A read with too few unique k-mers is `not_assessable`.

Only the sampled k-mers are kept, so memory is set by the sample size, not the genome.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
from collections import Counter
from pathlib import Path

import numpy as np

K = 31
STEP = 50
CHUNK = 20_000_000
MIN_VOTES = 3
SLACK = 1000
CODE = np.full(256, 4, dtype=np.uint8)
for index, base in enumerate(b"ACGT"):
    CODE[base] = index


def canonical_kmers(sequence: bytes, k: int = K) -> tuple[np.ndarray, np.ndarray]:
    """(canonical k-mer of each window, whether the window has no ambiguous base). k <= 31."""
    codes = CODE[np.frombuffer(sequence, dtype=np.uint8)]
    n = len(codes) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    ambiguous = np.concatenate(([0], np.cumsum(codes == 4)))
    valid = ambiguous[k:] == ambiguous[:-k]
    clipped = np.minimum(codes, 3).astype(np.uint64)
    forward = np.zeros(n, dtype=np.uint64)
    reverse = np.zeros(n, dtype=np.uint64)
    # One pass per k-mer position, so memory stays a few arrays rather than a k-by-n matrix.
    for offset in range(k):
        column = clipped[offset:offset + n]
        forward += column << np.uint64(2 * (k - 1 - offset))
        reverse += (np.uint64(3) - column) << np.uint64(2 * offset)
    return np.minimum(forward, reverse), valid


def read_fasta_bytes(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        yield name, b"".join(chunks)


def scan_reference(path: Path, targets: np.ndarray) -> tuple[np.ndarray, dict[int, tuple[str, int]]]:
    """Occurrence count of each sorted target k-mer, and the last place each was seen."""
    counts = np.zeros(len(targets), dtype=np.int64)
    places: dict[int, tuple[str, int]] = {}
    for chrom, sequence in read_fasta_bytes(path):
        for start in range(0, max(1, len(sequence) - K + 1), CHUNK):
            kmers, valid = canonical_kmers(sequence[start:start + CHUNK + K - 1])
            positions = np.flatnonzero(valid)
            kmers = kmers[positions]
            slots = np.searchsorted(targets, kmers)
            slots[slots == len(targets)] = 0
            hits = np.flatnonzero(targets[slots] == kmers)
            np.add.at(counts, slots[hits], 1)
            for slot, position in zip(slots[hits].tolist(), positions[hits].tolist()):
                places[slot] = (chrom, start + position)
    return counts, places


def verdict(votes: list[tuple[str, int]], chrom: str, start: int, end: int) -> str:
    if len(votes) < MIN_VOTES:
        return "not_assessable"
    agreeing = sum(1 for vote_chrom, position in votes if vote_chrom == chrom and start - SLACK <= position < end + SLACK)
    return "agree" if agreeing * 2 > len(votes) else "disagree"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--features", type=Path, nargs="+", required=True, help="Labelled rows from label_alignments.py.")
    parser.add_argument("--haplotypes", nargs="+", required=True, help="Haplotype FASTA files.")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--sample", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; refusing to overwrite it")

    reads: dict[str, dict] = {}
    for path in args.features:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["truth_ref_chrom"]:
                    reads.setdefault(row["read"], row)
    names = sorted(reads)
    sample = random.Random(args.seed).sample(names, min(args.sample, len(names)))
    sources = {name: sequence for path in args.haplotypes for name, sequence in read_fasta_bytes(Path(path))}

    read_kmers: dict[str, np.ndarray] = {}
    for name in sample:
        row = reads[name]
        kmers, valid = canonical_kmers(sources[row["source_contig"]][int(row["source_start"]):int(row["source_end"])])
        read_kmers[name] = np.unique(kmers[::STEP][valid[::STEP]])
    targets = np.unique(np.concatenate(list(read_kmers.values()))) if read_kmers else np.empty(0, dtype=np.uint64)
    counts, places = scan_reference(args.reference, targets)

    outcomes: Counter[str] = Counter()
    with gzip.open(args.output, "xt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["read", "unique_kmers", "verdict", "truth_chrom", "truth_start", "truth_end"])
        for name in sample:
            row = reads[name]
            slots = np.searchsorted(targets, read_kmers[name])
            votes = [places[int(slot)] for slot in slots if counts[slot] == 1]
            result = verdict(votes, row["truth_ref_chrom"], int(row["truth_ref_start"]), int(row["truth_ref_end"]))
            outcomes[result] += 1
            writer.writerow([name, len(votes), result, row["truth_ref_chrom"], row["truth_ref_start"], row["truth_ref_end"]])
    print(json.dumps({"sampled": len(sample), **outcomes}, indent=2))


if __name__ == "__main__":
    main()
