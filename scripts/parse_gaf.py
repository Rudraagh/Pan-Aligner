from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
from pathlib import Path

from common import write_json


@dataclass
class AlignmentResult:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    path: str
    path_length: int
    path_start: int
    path_end: int
    residue_matches: int
    alignment_block_length: int
    mapping_quality: int
    tags: dict[str, str] = field(default_factory=dict)

    @property
    def coverage(self) -> float:
        return self.alignment_block_length / self.query_length if self.query_length else 0.0

    @property
    def identity(self) -> float:
        return self.residue_matches / self.alignment_block_length if self.alignment_block_length else 0.0

    @property
    def path_span(self) -> int:
        return max(0, self.path_end - self.path_start)

    @property
    def alignment_score(self) -> float:
        if "AS" in self.tags:
            try:
                return float(self.tags["AS"])
            except ValueError:
                return 0.0
        return float(self.residue_matches)

    @property
    def normalized_score(self) -> float:
        if self.query_length == 0:
            return 0.0
        return self.alignment_score / self.query_length

    @property
    def traversed_nodes(self) -> list[str]:
        nodes: list[str] = []
        current = []
        for char in self.path:
            if char in "<>":
                if current:
                    nodes.append("".join(current))
                    current = []
            else:
                current.append(char)
        if current:
            nodes.append("".join(current))
        return nodes


def parse_gaf_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for tag_field in fields:
        parts = tag_field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_gaf_line(line: str) -> AlignmentResult:
    fields = line.rstrip().split("\t")
    if len(fields) < 12:
        raise ValueError(f"Malformed GAF line with {len(fields)} fields: {line}")
    return AlignmentResult(
        query_name=fields[0],
        query_length=int(fields[1]),
        query_start=int(fields[2]),
        query_end=int(fields[3]),
        strand=fields[4],
        path=fields[5],
        path_length=int(fields[6]),
        path_start=int(fields[7]),
        path_end=int(fields[8]),
        residue_matches=int(fields[9]),
        alignment_block_length=int(fields[10]),
        mapping_quality=int(fields[11]),
        tags=parse_gaf_tags(fields[12:]),
    )


def parse_gaf_file(gaf_path: Path) -> list[AlignmentResult]:
    results: list[AlignmentResult] = []
    with gaf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                results.append(parse_gaf_line(line))
    return results


def best_alignment(results: list[AlignmentResult]) -> AlignmentResult | None:
    if not results:
        return None
    return max(results, key=lambda item: (item.normalized_score, item.identity, item.coverage, item.mapping_quality))


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse a GAF file into JSON records.")
    parser.add_argument("gaf", type=Path, help="Input GAF file.")
    parser.add_argument("--output-json", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    results = parse_gaf_file(args.gaf.resolve())
    output = [asdict(item) | {"coverage": item.coverage, "identity": item.identity, "normalized_score": item.normalized_score, "traversed_nodes": item.traversed_nodes} for item in results]
    if args.output_json:
        write_json(args.output_json.resolve(), output)
    else:
        print(output)


if __name__ == "__main__":
    main()

