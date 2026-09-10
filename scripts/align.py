from __future__ import annotations

import argparse
from pathlib import Path

from common import OUTPUTS_DIR, ensure_dir, resolve_binary, run_command
from parse_gaf import AlignmentResult, best_alignment, parse_gaf_file


def run_panaligner(graph_gfa: Path, query_fasta: Path, output_gaf: Path, panaligner_bin: Path, threads: int) -> list[AlignmentResult]:
    ensure_dir(output_gaf.parent)
    command = [str(panaligner_bin), "-t", str(threads), "-cx", "lr", str(graph_gfa), str(query_fasta)]
    run_command(command, stdout_path=output_gaf)
    return parse_gaf_file(output_gaf)


def run_minigraph(graph_gfa: Path, query_fasta: Path, output_gaf: Path, minigraph_bin: Path, threads: int) -> list[AlignmentResult]:
    ensure_dir(output_gaf.parent)
    command = [str(minigraph_bin), "-t", str(threads), "-cx", "lr", str(graph_gfa), str(query_fasta)]
    run_command(command, stdout_path=output_gaf)
    return parse_gaf_file(output_gaf)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run PanAligner on a graph/query pair and summarize the best alignment.")
    parser.add_argument("--graph", required=True, type=Path, help="Input graph in GFA/rGFA format.")
    parser.add_argument("--query", required=True, type=Path, help="Query FASTA file.")
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to PanAligner.")
    parser.add_argument("--minigraph-bin", type=str, default=None, help="Optional minigraph binary for heuristic alignment.")
    parser.add_argument(
        "--panaligner-bin",
        type=str,
        default=None,
        help="Path to the PanAligner binary. Defaults to ./PanAligner/PanAligner if present.",
    )
    parser.add_argument(
        "--output-gaf",
        type=Path,
        default=OUTPUTS_DIR / "alignments" / "alignment.gaf",
        help="Path for the generated GAF output.",
    )
    args = parser.parse_args()

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    results = run_panaligner(args.graph.resolve(), args.query.resolve(), args.output_gaf.resolve(), panaligner_bin, args.threads)
    best = best_alignment(results)
    print(f"Alignments written to {args.output_gaf.resolve()}")
    if best:
        print(
            f"Best alignment: query={best.query_name} path={best.path} "
            f"identity={best.identity:.4f} coverage={best.coverage:.4f} score={best.normalized_score:.4f}"
        )
    else:
        print("No alignments reported.")


if __name__ == "__main__":
    main()
