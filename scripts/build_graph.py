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


def resolve_manifest_fasta_path(candidate: Path, gene: str, label: str | None = None, sequence_id: str | None = None) -> Path:
    if candidate.exists():
        return candidate

    if label is None:
        fallback = Path("data/raw") / f"{gene.lower()}_reference.fa"
    else:
        fallback = Path("data") / label.lower() / gene.lower() / f"{sequence_id}.fa"

    resolved_fallback = fallback.resolve()
    if resolved_fallback.exists():
        return resolved_fallback
    raise FileNotFoundError(f"Could not locate FASTA file. Checked {candidate} and {resolved_fallback}.")


def build_graphs(
    dataset_manifest_path: Path,
    minigraph_bin: Path,
    threads: int,
    output_dir: Path | None = None,
    output_manifest_path: Path | None = None,
) -> dict:
    preprocess_manifest = read_json(dataset_manifest_path)
    graphs_manifest = {"graphs": {}}

    target_graph_dir = ensure_dir(output_dir or GRAPHS_DIR)

    for gene, payload in preprocess_manifest["genes"].items():
        reference_fasta = resolve_manifest_fasta_path(Path(payload["reference"]["path"]), gene)
        healthy_fastas = [resolve_manifest_fasta_path(Path(entry["path"]), gene, "healthy", entry["id"]) for entry in payload["healthy"]]
        unhealthy_fastas = [resolve_manifest_fasta_path(Path(entry["path"]), gene, "unhealthy", entry["id"]) for entry in payload["unhealthy"]]
        combined_fastas = healthy_fastas + unhealthy_fastas

        gene_graphs = {
            "reference_fasta": str(reference_fasta.resolve()),
            "healthy_graph": str((target_graph_dir / f"{gene.lower()}.healthy.gfa").resolve()),
            "unhealthy_graph": str((target_graph_dir / f"{gene.lower()}.unhealthy.gfa").resolve()),
            "combined_graph": str((target_graph_dir / f"{gene.lower()}.combined.gfa").resolve()),
        }

        build_gene_graph(reference_fasta, healthy_fastas, Path(gene_graphs["healthy_graph"]), minigraph_bin, threads)
        build_gene_graph(reference_fasta, unhealthy_fastas, Path(gene_graphs["unhealthy_graph"]), minigraph_bin, threads)
        build_gene_graph(reference_fasta, combined_fastas, Path(gene_graphs["combined_graph"]), minigraph_bin, threads)
        graphs_manifest["graphs"][gene] = gene_graphs

    manifest_path = output_manifest_path or (METADATA_DIR / "graph_manifest.json")
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
    parser.add_argument("--output-dir", type=Path, default=GRAPHS_DIR, help="Directory where graph GFA files will be written.")
    parser.add_argument(
        "--output-manifest",
        type=Path,
        default=METADATA_DIR / "graph_manifest.json",
        help="Path where the graph manifest JSON will be written.",
    )
    parser.add_argument(
        "--minigraph-bin",
        type=str,
        default=None,
        help="Path to the minigraph binary. Defaults to ./minigraph/minigraph if present.",
    )
    args = parser.parse_args()

    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph"), Path("minigraph") / "minigraph"])
    manifest = build_graphs(
        args.preprocess_manifest.resolve(),
        minigraph_bin,
        args.threads,
        args.output_dir.resolve(),
        args.output_manifest.resolve(),
    )
    print(f"Graph construction complete. Manifest written to {args.output_manifest.resolve()}")
    print(f"Graphs built for: {', '.join(sorted(manifest['graphs']))}")


if __name__ == "__main__":
    main()
