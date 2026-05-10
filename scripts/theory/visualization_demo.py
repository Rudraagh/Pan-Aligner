from __future__ import annotations

import argparse
from pathlib import Path

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from theory.anchors_demo import run_anchors_demo
from theory.chaining_demo import run_chaining_demo
from theory.dag_converter_demo import run_dag_converter_demo
from theory.graph_analysis_demo import run_graph_analysis
from theory.path_cover_demo import run_path_cover_demo
from theory.precedence_demo import run_precedence_demo
from theory.scc_demo import run_scc_demo


def run_visualization_demo(graph_path: Path | None = None) -> None:
    run_graph_analysis(graph_path)
    run_scc_demo(graph_path)
    run_dag_converter_demo(graph_path)
    run_path_cover_demo(graph_path)
    run_anchors_demo(graph_path)
    run_precedence_demo(graph_path)
    run_chaining_demo(graph_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate all educational theory visualizations.")
    parser.add_argument("--graph", type=Path, default=None, help="Optional graph GFA.")
    args = parser.parse_args()
    run_visualization_demo(args.graph.resolve() if args.graph else None)
    print("All theory visualizations generated.")


if __name__ == "__main__":
    main()
