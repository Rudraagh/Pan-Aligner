"""CPU/GPU minimizer parity and minimizer-only benchmark on a local CUDA device.

Measures only minimizer extraction (scripts/proposed/minimizer_engine.py). Graph
construction, consistency scoring, anchors and feedback always run on the CPU and are
not timed here. Results describe the machine the script runs on; a laptop run is not
HPC GPU validation.

  parity      seeded random sequences (IUPAC symbols, k up to the kernel limit 55,
              short and empty sequences): GPU records must equal CPU records
  benchmark   seeded random DNA and the three gene reference FASTAs
              cpu_end_to_end / gpu_end_to_end   compute_minimizers() wall time, median
                                                of repetitions after warm-up
              gpu_phases_ms                     CUDA-event times per device phase
                                                (median), from _compute_gpu_minimizers
              gpu_host_total_ms                 wall time of the same call, including
                                                Python record construction
              gpu_first_call_s                  first call of this process for the case
              gpu_pool_bytes                    CuPy memory-pool bytes held after one call
              python_heap_peak_bytes            tracemalloc peak for one call (Python
                                                heap only, not process RSS)

Needs CuPy and a usable CUDA device (run with the GPU environment's Python):
  python scripts/gpu_minimizer_bench.py --output-dir outputs/gpu_validation/laptop
"""

from __future__ import annotations

import argparse
import platform
import random
import statistics
import subprocess
import sys
import time
import tracemalloc
from pathlib import Path

from common import OUTPUTS_DIR, ROOT, ensure_dir, project_relative, write_json
from proposed.minimizer_engine import _compute_gpu_minimizers, _load_usable_cupy, compute_minimizers
from simulate_truth_reads import read_fasta, sha256_file


SEED = 20261029
PARITY_ALPHABET = "ACGT" * 8 + "NRYKMSWBDHV"
RANDOM_LENGTHS = (10_000, 100_000, 1_000_000)
RANDOM_PARAMETERS = ((15, 10), (31, 20), (55, 10))
GENE_PARAMETERS = (11, 6)  # proposed final-evaluation configuration
GENES = ("APP", "PSEN1", "PSEN2")


def git(*args: str) -> str:
    completed = subprocess.run(["git", *args], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    return completed.stdout.strip()


def timed(function, repetitions: int, warmups: int) -> dict:
    for _ in range(warmups):
        function()
    samples = []
    for _ in range(repetitions):
        start = time.perf_counter()
        function()
        samples.append(time.perf_counter() - start)
    return {"median_s": statistics.median(samples), "min_s": min(samples), "repetitions": repetitions, "warmups": warmups}


def heap_peak(function) -> int:
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def run_parity(cases: int, rng: random.Random) -> dict:
    mismatches = []
    for _ in range(cases):
        sequence = "".join(rng.choice(PARITY_ALPHABET) for _ in range(rng.randint(0, 400)))
        k = rng.choice([1, 2, 3, 5, 11, 15, 21, 31, 54, 55, rng.randint(1, 55)])
        window = rng.randint(1, 60)
        cpu = compute_minimizers(sequence, k, window, "hap", "cpu")["records"]
        gpu = compute_minimizers(sequence, k, window, "hap", "gpu")
        if gpu["backend_used"] != "gpu" or gpu["records"] != cpu:
            mismatches.append({"length": len(sequence), "k": k, "window_size": window})
    return {"cases": cases, "mismatches": len(mismatches), "examples": mismatches[:5]}


def benchmark_case(cp, name: str, sequence: str, k: int, window: int, repetitions: int) -> dict:
    pool = cp.get_default_memory_pool()
    pool.free_all_blocks()
    first_start = time.perf_counter()
    gpu_payload = compute_minimizers(sequence, k, window, "hap", "gpu")
    first_call = time.perf_counter() - first_start
    cpu_payload = compute_minimizers(sequence, k, window, "hap", "cpu")
    identical = gpu_payload["backend_used"] == "gpu" and gpu_payload["records"] == cpu_payload["records"]
    gpu_pool_bytes = pool.total_bytes()

    cpu_repetitions = repetitions if len(sequence) <= 100_000 else max(2, repetitions // 2)
    cpu_time = timed(lambda: compute_minimizers(sequence, k, window, "hap", "cpu"), cpu_repetitions, 1)
    gpu_time = timed(lambda: compute_minimizers(sequence, k, window, "hap", "gpu"), repetitions, 2)

    phases: dict[str, list[float]] = {}
    host_totals = []
    for _ in range(repetitions):
        timings: dict = {}
        start = time.perf_counter()
        _compute_gpu_minimizers(sequence, k, window, "hap", cp, timings=timings)
        host_totals.append((time.perf_counter() - start) * 1000.0)
        for phase, value in timings.items():
            phases.setdefault(phase, []).append(value)

    return {
        "case": name, "length": len(sequence), "k": k, "window_size": window,
        "minimizers": len(cpu_payload["records"]), "identical": identical,
        "cpu_end_to_end": cpu_time, "gpu_end_to_end": gpu_time,
        "speedup_end_to_end": cpu_time["median_s"] / gpu_time["median_s"],
        "gpu_first_call_s": first_call,
        "gpu_phases_ms": {phase: statistics.median(values) for phase, values in phases.items()},
        "gpu_device_total_ms": sum(statistics.median(values) for values in phases.values()),
        "gpu_host_total_ms": statistics.median(host_totals),
        "gpu_pool_bytes": gpu_pool_bytes,
        "python_heap_peak_bytes": {
            "cpu": heap_peak(lambda: compute_minimizers(sequence, k, window, "hap", "cpu")),
            "gpu": heap_peak(lambda: compute_minimizers(sequence, k, window, "hap", "gpu")),
        },
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="CPU/GPU minimizer parity and minimizer-only benchmark.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "gpu_validation" / "laptop")
    parser.add_argument("--parity-cases", type=int, default=400)
    parser.add_argument("--repetitions", type=int, default=5)
    args = parser.parse_args(argv)

    cp, reason = _load_usable_cupy()
    if cp is None:
        raise SystemExit(f"No usable CUDA device; nothing was measured. {reason}")
    output_dir = ensure_dir(args.output_dir.resolve())
    rng = random.Random(SEED)
    started = time.perf_counter()

    parity = run_parity(args.parity_cases, rng)
    print(f"parity: {parity['cases']} cases, {parity['mismatches']} mismatches", flush=True)

    cases = []
    for length in RANDOM_LENGTHS:
        sequence = "".join(rng.choice("ACGT") for _ in range(length))
        for k, window in RANDOM_PARAMETERS:
            cases.append((f"random_{length}", sequence, k, window))
    reference_inputs = {}
    for gene in GENES:
        path = ROOT / "data" / "raw" / f"{gene.lower()}_reference.fa"
        reference_inputs[project_relative(path)] = sha256_file(path)
        cases.append((f"{gene}_reference", read_fasta(path)[0][1].upper(), *GENE_PARAMETERS))

    rows = []
    for name, sequence, k, window in cases:
        row = benchmark_case(cp, name, sequence, k, window, args.repetitions)
        rows.append(row)
        print(f"{name:>16} len={row['length']:>9} k={k:2} w={window:2} identical={row['identical']} "
              f"cpu={row['cpu_end_to_end']['median_s']:.4f}s gpu={row['gpu_end_to_end']['median_s']:.4f}s "
              f"device={row['gpu_device_total_ms']:.2f}ms speedup={row['speedup_end_to_end']:.2f}x", flush=True)

    properties = cp.cuda.runtime.getDeviceProperties(0)
    write_json(output_dir / "results.json", {"parity": parity, "benchmark": rows})
    write_json(output_dir / "run_manifest.json", {
        "command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "git_commit": git("rev-parse", "HEAD"),
        "git_uncommitted_code": git("status", "--porcelain", "--", "scripts", "tests", "proposed_methodology.py", "main.py").splitlines(),
        "seed": SEED,
        "parameters": {"parity_cases": args.parity_cases, "repetitions": args.repetitions,
                       "random_lengths": RANDOM_LENGTHS, "random_parameters": RANDOM_PARAMETERS, "gene_parameters": GENE_PARAMETERS},
        "reference_input_sha256": reference_inputs,
        "device": {
            "name": properties["name"].decode(),
            "compute_capability": f"{properties['major']}.{properties['minor']}",
            "total_memory_bytes": properties["totalGlobalMem"],
            "cuda_runtime": cp.cuda.runtime.runtimeGetVersion(),
            "cuda_driver": cp.cuda.runtime.driverGetVersion(),
        },
        "cupy": cp.__version__,
        "python": platform.python_version(),
        "host": platform.platform(),
        "processor": platform.processor(),
        "pip_freeze": subprocess.run([sys.executable, "-m", "pip", "freeze", "--disable-pip-version-check"],
                                     stdout=subprocess.PIPE, text=True, check=False).stdout.split(),
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "scope": [
            "Minimizer extraction only; graph construction and every later stage run on the CPU and are not measured.",
            "End-to-end times include Python record construction and host-device transfer.",
            "gpu_phases_ms are CUDA-event times on the default stream; python_heap_peak_bytes is tracemalloc, not RSS.",
            "Single local device; not HPC GPU validation and not a whole-pipeline speed-up.",
        ],
    })
    print(f"wrote {project_relative(output_dir)}")


if __name__ == "__main__":
    main()
