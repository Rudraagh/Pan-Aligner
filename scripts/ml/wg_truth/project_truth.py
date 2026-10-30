"""Project a read's source interval from a haplotype onto the reference through whole-assembly alignments.

The interval itself comes from the simulator. Alignments only translate its coordinates, and any read
that cannot be translated cleanly is reported as unprojectable rather than guessed.

Input alignments are PAF with a cg tag (minimap2 -c) or SAM (samtools view of the dipcall BAMs); secondary
records are ignored. Coordinates are 0-based half-open throughout.
"""

from __future__ import annotations

import bisect
import gzip
import re
from dataclasses import dataclass
from pathlib import Path

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
MAX_END_SNAP = 1000
MIN_SPAN_RATIO = 0.5
MAX_SPAN_RATIO = 2.0


@dataclass(frozen=True)
class RefInterval:
    chrom: str
    start: int
    end: int
    strand: str


@dataclass
class ContigAlignment:
    """One alignment. Matched blocks are sorted by forward-strand query position."""

    query: str
    query_start: int
    query_end: int
    target: str
    strand: str
    block_query_starts: list[int]
    block_query_ends: list[int]
    block_target_starts: list[int]

    def _lift(self, index: int, position: int) -> int:
        offset = position - self.block_query_starts[index]
        if self.strand == "+":
            return self.block_target_starts[index] + offset
        length = self.block_query_ends[index] - self.block_query_starts[index]
        return self.block_target_starts[index] + length - 1 - offset

    def snap_inward(self, position: int, direction: int) -> tuple[int, int] | None:
        """(target base, bases moved): lift `position`, moving towards the interval's inside (direction +1 or -1) if it sits in a gap."""
        index = bisect.bisect_right(self.block_query_starts, position) - 1
        if index >= 0 and position < self.block_query_ends[index]:
            return self._lift(index, position), 0
        if direction > 0:
            index += 1
            if index >= len(self.block_query_starts):
                return None
            moved = self.block_query_starts[index] - position
            return self._lift(index, self.block_query_starts[index]), moved
        if index < 0:
            return None
        last = self.block_query_ends[index] - 1
        return self._lift(index, last), position - last


def _build(query: str, query_length: int, query_start: int, target: str, target_start: int, strand: str, cigar: str) -> ContigAlignment:
    """`query_start` is the start of the aligned part on the aligned strand (clips are already skipped)."""
    blocks: list[tuple[int, int, int]] = []
    query_pos, target_pos = query_start, target_start
    for length_text, op in CIGAR_OP.findall(cigar):
        length = int(length_text)
        if op in "M=X":
            if blocks and blocks[-1][0] + blocks[-1][2] == query_pos and blocks[-1][1] + blocks[-1][2] == target_pos:
                blocks[-1] = (blocks[-1][0], blocks[-1][1], blocks[-1][2] + length)
            else:
                blocks.append((query_pos, target_pos, length))
            query_pos += length
            target_pos += length
        elif op == "I":
            query_pos += length
        elif op in "DN":
            target_pos += length
    if strand == "-":
        blocks = sorted((query_length - q - n, t, n) for q, t, n in blocks)
        first, last = query_length - query_pos, query_length - query_start
    else:
        first, last = query_start, query_pos
    return ContigAlignment(query, first, last, target, strand, [b[0] for b in blocks], [b[0] + b[2] for b in blocks], [b[1] for b in blocks])


def parse_paf_line(line: str) -> ContigAlignment | None:
    fields = line.rstrip("\n").split("\t")
    tags = dict((tag[:2], tag[5:]) for tag in fields[12:])
    if tags.get("tp", "P") not in ("P", "I"):
        return None
    if "cg" not in tags:
        raise ValueError(f"PAF record without a cg tag for {fields[0]}; run minimap2 with -c")
    query_length, query_start, query_end = int(fields[1]), int(fields[2]), int(fields[3])
    strand = fields[4]
    aligned_start = query_start if strand == "+" else query_length - query_end
    return _build(fields[0], query_length, aligned_start, fields[5], int(fields[7]), strand, tags["cg"])


def parse_sam_line(line: str) -> ContigAlignment | None:
    if line.startswith("@"):
        return None
    fields = line.split("\t", 6)
    flag = int(fields[1])
    if flag & 0x4 or flag & 0x100:
        return None
    cigar = fields[5]
    ops = CIGAR_OP.findall(cigar)
    if not ops:
        return None
    query_length = sum(int(n) for n, op in ops if op in "MIS=XH")
    clip = int(ops[0][0]) if ops[0][1] in "SH" else 0
    return _build(fields[0], query_length, clip, fields[2], int(fields[3]) - 1, "-" if flag & 0x10 else "+", cigar)


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def load_alignments(paths: list[Path]) -> dict[str, list[ContigAlignment]]:
    """Alignments by query contig, from .paf or .sam files, optionally gzipped."""
    by_query: dict[str, list[ContigAlignment]] = {}
    for path in paths:
        parse = parse_sam_line if ".sam" in Path(path).name else parse_paf_line
        with _open_text(Path(path)) as handle:
            for line in handle:
                record = parse(line)
                if record is not None and record.block_query_starts:
                    by_query.setdefault(record.query, []).append(record)
    return by_query


def project_interval(alignments: list[ContigAlignment], start: int, end: int, read_strand: str) -> tuple[RefInterval | None, str]:
    """Reference interval for source interval [start, end), or (None, reason)."""
    containing = [a for a in alignments if a.query_start <= start and end <= a.query_end]
    if not containing:
        return None, "no_single_alignment"
    if len(containing) > 1:
        return None, "multiple_alignments"
    alignment = containing[0]
    first = alignment.snap_inward(start, +1)
    last = alignment.snap_inward(end - 1, -1)
    if first is None or last is None:
        return None, "end_outside_blocks"
    if first[1] > MAX_END_SNAP or last[1] > MAX_END_SNAP:
        return None, "end_in_insertion"
    if alignment.strand == "+":
        target_start, target_end = first[0], last[0] + 1
    else:
        target_start, target_end = last[0], first[0] + 1
    span = target_end - target_start
    if span <= 0:
        return None, "end_in_insertion"
    ratio = span / (end - start)
    if not MIN_SPAN_RATIO <= ratio <= MAX_SPAN_RATIO:
        return None, "large_structural_difference"
    strand = read_strand if alignment.strand == "+" else {"+": "-", "-": "+"}[read_strand]
    return RefInterval(alignment.target, target_start, target_end, strand), "ok"
