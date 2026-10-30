"""Project a GAF path from a whole-genome rGFA onto the reference chromosomes.

`clean`: the alignment starts and ends on reference segments.
`flank`: it starts or ends on non-reference nodes, so the interval runs from the nearest reference segments
         and is widened by the bases skipped. Reported apart from `clean`, because its ends are approximate.
Paths with no reference segment, with ends on different chromosomes or in different orientations, or whose
ends are out of order on the reference cannot be projected, and return None.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from truth_eval import STABLE_SEGMENT  # noqa: E402

FLIP = {"+": "-", "-": "+"}


@dataclass(frozen=True)
class PathProjection:
    chrom: str
    start: int
    end: int
    strand: str
    kind: str


def project_path(strand: str, path: str, path_start: int, path_end: int, reference_names: set[str]) -> PathProjection | None:
    """`strand`, `path`, `path_start` and `path_end` are GAF columns 5, 6, 8 and 9."""
    if path[0] not in "<>":
        if path not in reference_names:
            return None
        return PathProjection(path, path_start, path_end, strand, "clean")
    segments = STABLE_SEGMENT.findall(path)
    if not segments or "".join(f"{o}{n}:{s}-{e}" for o, n, s, e in segments) != path:
        return None
    offsets, position = [], 0
    for _, _, seg_start, seg_end in segments:
        offsets.append(position)
        position += int(seg_end) - int(seg_start)
    on_reference = [i for i, segment in enumerate(segments) if segment[1] in reference_names]
    if not on_reference:
        return None
    first_index, last_index = on_reference[0], on_reference[-1]
    first, last = segments[first_index], segments[last_index]
    if first[1] != last[1] or first[0] != last[0]:
        return None
    after_last = offsets[last_index] + int(last[3]) - int(last[2])
    before = path_start - offsets[first_index]
    beyond = after_last - path_end
    if first[0] == ">":
        start, end, out_strand = int(first[2]) + before, int(last[3]) - beyond, strand
    else:
        start, end, out_strand = int(last[2]) + beyond, int(first[3]) - before, FLIP[strand]
    start = max(0, start)
    if start >= end:
        return None
    kind = "clean" if first_index == 0 and last_index == len(segments) - 1 else "flank"
    return PathProjection(first[1], start, end, out_strand, kind)
