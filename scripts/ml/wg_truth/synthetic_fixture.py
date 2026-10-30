"""Write a small synthetic genome for testing the pipeline end to end with the real tools. Not biological data.

Reference: chr1 and chr21, sharing a near-identical 20 kb paralog, so misplaced reads are possible.
Graph haplotype: chr21 with an insertion, added to the graph, which creates non-reference nodes.
Source haplotype: chr21 with SNPs, small indels and a large insertion, kept out of the graph.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from simulate_truth_reads import write_fasta_with_fai  # noqa: E402


def random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choices("ACGT", k=length))


def mutate(rng: random.Random, sequence: str, rate: float, indel_rate: float) -> str:
    out = []
    for base in sequence:
        draw = rng.random()
        if draw < rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        elif draw < rate + indel_rate / 2:
            continue
        elif draw < rate + indel_rate:
            out.append(base + rng.choice("ACGT"))
        else:
            out.append(base)
    return "".join(out)


def main() -> None:
    output = Path(sys.argv[1])
    rng = random.Random(20261030)
    paralog = random_sequence(rng, 20_000)
    left, right = random_sequence(rng, 150_000), random_sequence(rng, 230_000)
    chr21 = left + paralog + right
    chr1 = random_sequence(rng, 90_000) + mutate(rng, paralog, 0.001, 0.0) + random_sequence(rng, 190_000)
    shared_insertion = random_sequence(rng, 3_000)
    insert_at = 130_000  # within `right`, so at 300,000 on chr21
    graph_hap = left + paralog + right[:insert_at] + shared_insertion + right[insert_at:]
    # The source also carries the graph's insertion, so its reads must pass through a non-reference node.
    source_left = mutate(rng, left[:60_000], 0.001, 0.0002) + random_sequence(rng, 8_000) + mutate(rng, left[60_000:], 0.001, 0.0002)
    source_paralog = mutate(rng, paralog, 0.0005, 0.0)
    source = source_left + source_paralog + mutate(rng, right[:insert_at], 0.001, 0.0002) + shared_insertion + mutate(rng, right[insert_at:], 0.001, 0.0002)
    write_fasta_with_fai([("chr1", chr1), ("chr21", chr21)], output / "reference.fa")
    write_fasta_with_fai([("graphhap#1#chr21", graph_hap)], output / "graph_hap.fa")
    write_fasta_with_fai([("chr21_MATERNAL", source)], output / "source.fa")
    start = len(source_left)
    (output / "ambiguous.bed").write_text(f"chr21_MATERNAL\t{start}\t{start + len(source_paralog)}\n", encoding="ascii", newline="\n")


if __name__ == "__main__":
    main()
