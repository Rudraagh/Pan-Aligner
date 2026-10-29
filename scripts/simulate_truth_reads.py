"""Simulate truth-labelled long reads with PBSIM2 for the local gene regions.

Follows PanAligner/data/simulated_reads.sh (the paper's script): fixed 10 kb reads,
95% accuracy, P6C4 HMM model, and `paftools.js pbsim2fq` to put the true interval in
each read name (`<id>!<ref>!<start>!<end>!<strand>`). Coverage is a parameter because
the paper's 0.5x over a whole genome gives only a handful of reads on a single locus.

Sources:
  reference  reads come from the gene's reference backbone (as the paper simulated from CHM13)
  heldout    reads come from the held-out test haplotypes, which are not in the train graphs;
             their true intervals are lifted to reference coordinates (exact, see build_liftover)

Linux only (runs pbsim and k8); on the laptop run it inside WSL.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
from pathlib import Path
from typing import Callable

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, project_relative, resolve_project_path, write_json


GENES = ("APP", "PSEN1", "PSEN2")
DEFAULT_PBSIM = ROOT / ".local-build" / "pbsim2" / "bin" / "pbsim"
DEFAULT_MODEL = ROOT / ".local-build" / "pbsim2" / "data" / "P6C4.model"
K8_BINARY = ROOT / "tools" / "third_party" / "k8-bin" / "k8-1.2" / "k8-x86_64-Linux"
PAFTOOLS_JS = ROOT / "tools" / "third_party" / "minimap2" / "misc" / "paftools.js"
TRUTH_DIR = OUTPUTS_DIR / "truth_eval"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Records as (first header token, upper-case sequence)."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def write_fasta_with_fai(records: list[tuple[str, str]], path: Path) -> None:
    """Unwrapped FASTA plus a samtools-compatible .fai (samtools is not needed)."""
    ensure_dir(path.parent)
    offset = 0
    fasta_lines: list[str] = []
    fai_lines: list[str] = []
    for name, sequence in records:
        header = f">{name}\n"
        offset += len(header.encode("ascii"))
        fai_lines.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n")
        fasta_lines.append(header + sequence + "\n")
        offset += len(sequence) + 1
    path.write_text("".join(fasta_lines), encoding="ascii", newline="\n")
    Path(f"{path}.fai").write_text("".join(fai_lines), encoding="ascii", newline="\n")


def build_liftover(reference: str, haplotype: str, max_window: int = 50) -> Callable[[int], int]:
    """Map haplotype coordinates onto the reference.

    Equal-length haplotypes (substitutions only) map one-to-one. Otherwise the
    haplotype must differ from the reference inside one short window bounded by a
    shared prefix and suffix; coordinates before the window are unchanged, after
    it they shift by the length difference, and inside it they clamp to the window
    start. Anything else is refused rather than approximated.
    """
    if len(reference) == len(haplotype):
        return lambda position: position
    limit = min(len(reference), len(haplotype))
    prefix = 0
    while prefix < limit and reference[prefix] == haplotype[prefix]:
        prefix += 1
    suffix = 0
    while suffix < limit - prefix and reference[-1 - suffix] == haplotype[-1 - suffix]:
        suffix += 1
    hap_window_end = len(haplotype) - suffix
    ref_window_end = len(reference) - suffix
    if hap_window_end - prefix > max_window or ref_window_end - prefix > max_window:
        raise ValueError(
            f"haplotype differs from the reference over more than {max_window} bp "
            f"(prefix {prefix}, suffix {suffix}); exact liftover is not possible"
        )
    shift = len(reference) - len(haplotype)

    def lift(position: int) -> int:
        if position <= prefix:
            return position
        if position >= hap_window_end:
            return position + shift
        return prefix

    return lift


def parse_truth_name(name: str) -> tuple[str, str, int, int, str]:
    """Split a pbsim2fq read name `<id>!<target>!<start>!<end>!<strand>`."""
    read_id, target, start, end, strand = name.split("!")
    if strand not in "+-" or len(strand) != 1:
        raise ValueError(f"bad strand in read name: {name}")
    return read_id, target, int(start), int(end), strand


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_pbsim(pbsim: Path, model: Path, source_fasta: Path, work_dir: Path, seed: int, args: argparse.Namespace) -> list[Path]:
    command = [
        str(pbsim),
        "--depth", str(args.depth),
        "--length-min", str(args.length), "--length-mean", str(args.length), "--length-max", str(args.length),
        "--accuracy-mean", str(args.accuracy), "--accuracy-max", str(args.accuracy), "--accuracy-min", str(args.accuracy),
        "--seed", str(seed),
        "--hmm_model", str(model),
        "--prefix", "sim",
        str(source_fasta),
    ]
    with (work_dir / "pbsim.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=work_dir, stdout=log, stderr=subprocess.STDOUT, check=True)
    return sorted(work_dir.glob("sim_*.maf"))


def run_pbsim2fq(source_fasta: Path, maf_files: list[Path]) -> list[tuple[str, str]]:
    completed = subprocess.run(
        [str(K8_BINARY), str(PAFTOOLS_JS), "pbsim2fq", f"{source_fasta}.fai", *map(str, maf_files)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True,
    )
    records: list[tuple[str, str]] = []
    lines = completed.stdout.splitlines()
    for header, sequence in zip(lines[0::2], lines[1::2]):
        records.append((header[1:], sequence))
    return records


def heldout_haplotypes(gene: str, test_manifest: Path) -> list[tuple[str, str]]:
    manifest = json.loads(test_manifest.read_text(encoding="utf-8"))
    payload = manifest["genes"][gene]
    records: list[tuple[str, str]] = []
    for bucket in ("healthy", "unhealthy"):
        for entry in payload[bucket]:
            records.extend(read_fasta(resolve_project_path(entry["path"])))
    return records


def simulate_gene(gene: str, seed: int, args: argparse.Namespace) -> dict:
    reference_name, reference_seq = read_fasta(ROOT / "data" / "raw" / f"{gene.lower()}_reference.fa")[0]
    gene_dir = ensure_dir(args.output_dir / "reads" / args.source / gene)
    reference_fasta = gene_dir / "reference.fa"
    write_fasta_with_fai([(reference_name, reference_seq)], reference_fasta)

    if args.source == "reference":
        sources = [(reference_name, reference_seq)]
    else:
        sources = heldout_haplotypes(gene, args.test_manifest)
    liftovers = {name: build_liftover(reference_seq, sequence) for name, sequence in sources}

    seed_dir = ensure_dir(gene_dir / f"seed{seed}")
    with tempfile.TemporaryDirectory() as tmp:
        work_dir = Path(tmp)
        source_fasta = work_dir / "source.fa"
        write_fasta_with_fai(sources, source_fasta)
        maf_files = run_pbsim(args.pbsim, args.model, source_fasta, work_dir, seed, args)
        simulated = run_pbsim2fq(source_fasta, maf_files)
        (seed_dir / "pbsim.log").write_text((work_dir / "pbsim.log").read_text(encoding="utf-8"), encoding="utf-8")

    reads: list[tuple[str, str]] = []
    for name, sequence in simulated:
        read_id, source_name, start, end, strand = parse_truth_name(name)
        lift = liftovers[source_name]
        prefix = "" if args.source == "reference" else f"{source_name}."
        reads.append((f"{prefix}{read_id}!{reference_name}!{lift(start)}!{lift(end)}!{strand}", sequence))
    reads_fasta = seed_dir / "reads.fa"
    reads_fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in reads), encoding="ascii", newline="\n")

    return {
        "gene": gene,
        "seed": seed,
        "source": args.source,
        "source_sequences": [name for name, _ in sources],
        "source_bases": sum(len(sequence) for _, sequence in sources),
        "reference_name": reference_name,
        "reference_fasta": project_relative(reference_fasta),
        "reads_fasta": project_relative(reads_fasta),
        "reads": len(reads),
        "reads_sha256": sha256_file(reads_fasta),
        "pbsim_log": project_relative(seed_dir / "pbsim.log"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Simulate truth-labelled PBSIM2 reads for the local gene regions.")
    parser.add_argument("--source", choices=("reference", "heldout"), default="reference")
    parser.add_argument("--genes", nargs="+", default=list(GENES), choices=GENES)
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3])
    parser.add_argument("--depth", type=float, default=30.0, help="Coverage per source sequence (paper: 0.5 over a whole genome).")
    parser.add_argument("--length", type=int, default=10000, help="Fixed read length (paper: 10000).")
    parser.add_argument("--accuracy", type=float, default=0.95, help="Read accuracy (paper: 0.95).")
    parser.add_argument("--pbsim", type=Path, default=DEFAULT_PBSIM)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json")
    parser.add_argument("--output-dir", type=Path, default=TRUTH_DIR)
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()

    runs = [simulate_gene(gene, seed, args) for gene in args.genes for seed in args.seeds]
    manifest = {
        "source": args.source,
        "parameters": {
            "depth": args.depth, "length": args.length, "accuracy": args.accuracy,
            "model": args.model.name, "seeds": args.seeds,
        },
        "paper_parameters": {"depth": 0.5, "length": 10000, "accuracy": 0.95, "model": "P6C4.model"},
        "tools": {
            "pbsim": project_relative(args.pbsim),
            "pbsim_sha256": sha256_file(args.pbsim),
            "paftools_js": project_relative(PAFTOOLS_JS),
            "k8": project_relative(K8_BINARY),
        },
        "host": platform.platform(),
        "runs": runs,
    }
    write_json(args.output_dir / "reads" / args.source / "simulation_manifest.json", manifest)
    for run in runs:
        print(f"{run['source']} {run['gene']} seed{run['seed']}: {run['reads']} reads -> {run['reads_fasta']}")


if __name__ == "__main__":
    main()
