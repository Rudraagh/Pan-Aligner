"""Simulate reads from each haplotype for each setting, with PBSIM2.

The read name carries the origin on the haplotype: `<hap>.<id>!<contig>!<start>!<end>!<strand>`, written by
`paftools.js pbsim2fq` from PBSIM2's own MAF. No aligner is involved. Reads are streamed to disk and the MAF
is kept compressed beside them. Refuses to write into a setting that already has reads. Linux only.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from ml.build_alignment_dataset import setting_name  # noqa: E402
from simulate_truth_reads import K8_BINARY, PAFTOOLS_JS, sha256_file  # noqa: E402


def setting_seed(base_seed: int, accuracy: float, length: int, haplotype_index: int) -> int:
    """A separate, reproducible seed for each haplotype and setting."""
    return base_seed + 1000 * round(accuracy * 100) + length + 7 * haplotype_index


def simulate_one(source: Path, haplotype: str, accuracy: float, length: int, depth: float, seed: int, pbsim: Path, model: Path, output_dir: Path, k8: Path, paftools: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    reads_path = output_dir / f"{haplotype}.reads.fa"
    if reads_path.exists():
        raise SystemExit(f"{reads_path} already exists; refusing to overwrite it")
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(dir=output_dir) as work:
        command = [str(pbsim.resolve()), "--depth", str(depth),
                   "--length-min", str(length), "--length-mean", str(length), "--length-max", str(length),
                   "--accuracy-mean", str(accuracy), "--accuracy-max", str(accuracy), "--accuracy-min", str(accuracy),
                   "--seed", str(seed), "--hmm_model", str(model.resolve()), "--prefix", "sim", str(source.resolve())]
        with (output_dir / f"{haplotype}.pbsim.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
        mafs = sorted(Path(work).glob("sim_*.maf"))
        if not mafs:
            raise SystemExit(f"PBSIM2 wrote no MAF for {source}")
        reads = bases = 0
        partial = reads_path.with_suffix(".fa.partial")
        converter = subprocess.Popen([str(k8), str(paftools), "pbsim2fq", f"{source.resolve()}.fai", *map(str, mafs)], stdout=subprocess.PIPE, text=True)
        with partial.open("w", encoding="ascii", newline="\n") as out:
            for line in converter.stdout:
                if line.startswith(">"):
                    out.write(f">{haplotype}.{line[1:]}")
                    reads += 1
                else:
                    out.write(line)
                    bases += len(line) - 1
        if converter.wait() != 0:
            raise SystemExit(f"pbsim2fq failed for {source}")
        partial.rename(reads_path)
        for maf in mafs:
            with maf.open("rb") as raw, gzip.open(output_dir / f"{haplotype}.{maf.name}.gz", "xb") as packed:
                shutil.copyfileobj(raw, packed)
    return {"haplotype": haplotype, "source": str(source), "source_sha256": sha256_file(source), "accuracy": accuracy, "read_length": length,
            "depth": depth, "seed": seed, "reads": reads, "bases": bases, "reads_file": reads_path.name, "reads_sha256": sha256_file(reads_path),
            "seconds": round(time.perf_counter() - started, 3), "command": command}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--haplotypes", nargs="+", required=True, help="name=path.fa; each FASTA needs a .fai and is treated separately.")
    parser.add_argument("--accuracies", nargs="+", type=float, default=[0.95, 0.85, 0.75])
    parser.add_argument("--lengths", nargs="+", type=int, default=[10000, 2000])
    parser.add_argument("--depth", type=float, default=1.0, help="Coverage per haplotype.")
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--pbsim", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--k8", type=Path, default=K8_BINARY)
    parser.add_argument("--paftools", type=Path, default=PAFTOOLS_JS)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    haplotypes = [item.split("=", 1) for item in args.haplotypes]
    records = []
    for accuracy in args.accuracies:
        for length in args.lengths:
            setting = setting_name(accuracy, length)
            for index, (name, path) in enumerate(haplotypes):
                record = simulate_one(Path(path), name, accuracy, length, args.depth, setting_seed(args.seed, accuracy, length, index), args.pbsim, args.model, args.output_dir / setting, args.k8, args.paftools)
                records.append({"setting": setting, **record})
                print(setting, name, record["reads"], "reads", record["seconds"], "s", flush=True)
    manifest = args.output_dir / "simulation_manifest.json"
    with manifest.open("x", encoding="utf-8") as handle:
        json.dump({"base_seed": args.seed, "depth": args.depth, "model": str(args.model), "model_sha256": sha256_file(args.model), "runs": records}, handle, indent=2)


if __name__ == "__main__":
    main()
