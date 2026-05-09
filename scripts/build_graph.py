from __future__ import annotations

import argparse
from pathlib import Path

from common import GRAPHS_DIR, METADATA_DIR, ensure_dir, read_json, resolve_binary, run_command, write_json


def build_gene_graph(reference_fasta: Path, sample_fastas: list[Path], output_gfa: Path, minigraph_bin: Path, threads: int) -> None:
    if not sample_fastas:
        raise ValueError(f"No sample FASTA files provided for graph output {output_gfa}")
    command = [
        str(minigraph_bin),
        "-t",
        str(threads),
        "-cxggs",
        "-l1k",
        "-L1",
        str(reference_fasta),
        *[str(path) for path in sample_fastas],
    ]
    run_command(command, stdout_path=output_gfa)


def build_graphs(preprocess_manifest_path: Path, minigraph_bin: Path, threads: int) -> dict:
    preprocess_manifest = read_json(preprocess_manifest_path)
    graphs_manifest = {"graphs": {}}

    ensure_dir(GRAPHS_DIR)

    for gene, payload in preprocess_manifest["genes"].items():
        reference_fasta = Path(payload["reference"]["path"])
        healthy_fastas = [Path(entry["path"]) for entry in payload["healthy"]]
        unhealthy_fastas = [Path(entry["path"]) for entry in payload["unhealthy"]]
        combined_fastas = healthy_fastas + unhealthy_fastas

        gene_graphs = {
            "reference_fasta": str(reference_fasta.resolve()),
            "healthy_graph": str((GRAPHS_DIR / f"{gene.lower()}.healthy.gfa").resolve()),
            "unhealthy_graph": str((GRAPHS_DIR / f"{gene.lower()}.unhealthy.gfa").resolve()),
            "combined_graph": str((GRAPHS_DIR / f"{gene.lower()}.combined.gfa").resolve()),
        }

        build_gene_graph(reference_fasta, healthy_fastas, Path(gene_graphs["healthy_graph"]), minigraph_bin, threads)
        build_gene_graph(reference_fasta, unhealthy_fastas, Path(gene_graphs["unhealthy_graph"]), minigraph_bin, threads)
        build_gene_graph(reference_fasta, combined_fastas, Path(gene_graphs["combined_graph"]), minigraph_bin, threads)
        graphs_manifest["graphs"][gene] = gene_graphs

    manifest_path = METADATA_DIR / "graph_manifest.json"
    write_json(manifest_path, graphs_manifest)
    return graphs_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build healthy, unhealthy, and combined pangenome graphs with minigraph.")
    parser.add_argument(
        "--preprocess-manifest",
        type=Path,
        default=METADATA_DIR / "preprocess_manifest.json",
        help="Path to the preprocess manifest JSON.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Thread count passed to minigraph.")
    parser.add_argument(
        "--minigraph-bin",
        type=str,
        default=None,
        help="Path to the minigraph binary. Defaults to ./minigraph/minigraph if present.",
    )
    args = parser.parse_args()

    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph"), Path("minigraph") / "minigraph"])
    manifest = build_graphs(args.preprocess_manifest.resolve(), minigraph_bin, args.threads)
    print(f"Graph construction complete. Manifest written to {(METADATA_DIR / 'graph_manifest.json').resolve()}")
    print(f"Graphs built for: {', '.join(sorted(manifest['graphs']))}")


if __name__ == "__main__":
    main()
