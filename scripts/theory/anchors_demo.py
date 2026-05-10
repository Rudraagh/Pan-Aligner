from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.common import THEORY_OUTPUT_DIR, DemoGraphContext, build_demo_context, draw_graph, write_text


@dataclass
class Anchor:
    anchor_id: str
    vertex: str
    x: int
    y: int
    c: int
    d: int
    weight: int
    note: str


def _create_noise_anchors(context: DemoGraphContext, query_length: int) -> list[Anchor]:
    noise_anchors: list[Anchor] = []
    excluded = set(context.walk_nodes)
    other_nodes = [node for node in context.cyclic_demo_graph.nodes() if node not in excluded]
    for index, node in enumerate(other_nodes[:2], start=1):
        sequence = context.cyclic_demo_graph.nodes[node].get("sequence", "A")
        span = min(6, len(sequence))
        noise_anchors.append(
            Anchor(
                anchor_id=f"N{index}",
                vertex=node,
                x=1,
                y=span,
                c=max(1, query_length - (index * 4)),
                d=max(1, query_length - (index * 4) + span - 1),
                weight=max(3, span - 2),
                note="synthetic distractor anchor",
            )
        )
    return noise_anchors


def generate_demo_anchors(context: DemoGraphContext) -> list[Anchor]:
    anchors: list[Anchor] = []
    query_index = 1
    for index, node in enumerate(context.walk_nodes, start=1):
        sequence = context.cyclic_demo_graph.nodes[node].get("sequence", "A")
        span = min(8, len(sequence))
        anchors.append(
            Anchor(
                anchor_id=f"A{index}",
                vertex=node,
                x=1,
                y=span,
                c=query_index,
                d=query_index + span - 1,
                weight=span,
                note="walk anchor",
            )
        )
        query_index += span
        if index < len(context.walk_nodes) and index % 2 == 1:
            query_index += 1
    anchors.extend(_create_noise_anchors(context, max(1, len(context.query_sequence))))
    return sorted(anchors, key=lambda item: (item.c, item.d, item.anchor_id))


def write_anchor_csv(anchors: list[Anchor], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["anchor_id", "vertex", "x", "y", "c", "d", "weight", "note"],
        )
        writer.writeheader()
        for anchor in anchors:
            writer.writerow(anchor.__dict__)


def run_anchors_demo(graph_path: Path | None = None, output_dir: Path | None = None) -> dict:
    theory_output_dir = output_dir or THEORY_OUTPUT_DIR
    context = build_demo_context(graph_path)
    anchors = generate_demo_anchors(context)
    write_anchor_csv(anchors, theory_output_dir / "anchors.csv")

    node_anchor_counts = {node: 0 for node in context.cyclic_demo_graph.nodes()}
    for anchor in anchors:
        node_anchor_counts[anchor.vertex] = node_anchor_counts.get(anchor.vertex, 0) + 1
    labels = {
        node: f"{node}\nanchors={node_anchor_counts[node]}"
        for node in context.cyclic_demo_graph.nodes()
    }
    draw_graph(
        context.cyclic_demo_graph,
        theory_output_dir / "anchors_graph.png",
        "Anchor placement on educational graph",
        labels=labels,
    )

    lines = [
        "Anchor Representation Demo",
        "==========================",
        "Each anchor uses the simplified PanAligner-style tuple (vertex, [x..y], [c..d], weight).",
        f"Query sequence used for the demo: {context.query_sequence}",
        "",
    ]
    for anchor in anchors:
        lines.append(
            f"{anchor.anchor_id}: ({anchor.vertex}, [{anchor.x}..{anchor.y}], [{anchor.c}..{anchor.d}], {anchor.weight}) [{anchor.note}]"
        )
    write_text(theory_output_dir / "anchor_examples.txt", "\n".join(lines) + "\n")
    return {"context": context, "anchors": anchors}


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate simplified anchors for the PanAligner theory demo.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    result = run_anchors_demo(args.graph.resolve() if args.graph else None)
    print(f"Anchor count: {len(result['anchors'])}")


if __name__ == "__main__":
    main()
