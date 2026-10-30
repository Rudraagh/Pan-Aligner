# GPU minimizer engine

`scripts/proposed/minimizer_engine.py` retains the CPU implementation as the reference and optionally uses CuPy with CUDA for minimizer extraction.

## What runs on the GPU

For `backend="gpu"`, CUDA kernels compute the SHA-1-prefix key for every k-mer and the minimum key in every minimizer window. Consecutive repeated selected positions are compacted on the device. Only the compact minimizer positions and keys are copied to Python to create the existing `MinimizerRecord` objects. Graph construction remains CPU-side.

The GPU key is exactly the first 64 bits of the SHA-1 digest used by the CPU reference. Equal keys retain the earliest k-mer position, matching the CPU `(hash, position, kmer)` ordering. The current CUDA SHA-1 kernel supports ASCII sequences, `k <= 55` and sequences shorter than 2^31 bases; this is sufficient for the prototype defaults. Kernel scalar arguments are passed as `numpy.int32` to match the kernels' `int` parameters. `auto` uses CPU if CUDA is unavailable or this GPU kernel cannot execute. Explicit `gpu` mode raises a clear error instead of silently using CPU. The graph records the requested and used backend per stage in `metadata["backend"]`; every stage after minimizer extraction is CPU.

## Installation and use

CPU-only installation remains `pip install -r requirements.txt`. For CUDA, install the CuPy package matching the installed CUDA runtime, for example the package documented at <https://docs.cupy.dev/en/stable/install.html>. CuPy deliberately remains optional because its package name depends on the user’s CUDA version.

Run with `python proposed_methodology.py --backend cpu` or `python proposed_methodology.py --backend gpu`. `--backend auto` attempts a usable CUDA device and otherwise records an explicit CPU fallback reason in the minimizer payload.

Run the focused tests with:

```powershell
python -m unittest tests/test_proposed_minimizer_engine.py -v
```

The CUDA tests are skipped with their device/CuPy reason when CUDA is unavailable. They compare ordered records, including hashes, positions, and window starts, on fixed and 120 seeded random cases; check that kernel scalars are `int32`; check that phase timings are recorded; and (in `tests/test_proposed_backend_selection.py`) compare GPU-built and CPU-built graphs. No GPU speed claim should be made unless these tests pass on the target machine and a separate benchmark is run. GPU minimizer execution does not accelerate graph construction, scoring, anchor selection, or the baseline PanAligner workflow.

## Parity and benchmark

`scripts/gpu_minimizer_bench.py` (run with a Python that has CuPy) checks parity on 400 seeded random cases and benchmarks seeded random DNA (10 kb, 100 kb, 1 Mb) and the APP/PSEN1/PSEN2 references. It writes `results.json` and `run_manifest.json` (commit, uncommitted code, seed, input checksums, device, CUDA versions, `pip freeze`). End-to-end times include Python record construction and transfers; `gpu_phases_ms` are CUDA-event times for upload, hash kernel, window kernel, compaction and download; memory is reported as CuPy pool bytes and `tracemalloc` Python-heap peaks, not process RSS.

## Measured results: one laptop GPU

From `outputs/gpu_validation/laptop/` (commit `f2d3530`, no uncommitted code): NVIDIA GeForce RTX 3070 Ti Laptop GPU, CuPy 14.2.0, CUDA 13.2 runtime / 13.3 driver, Python 3.13.7, Windows 11. The environment was an isolated virtual environment (`cupy-cuda13x==14.2.0`), not part of `requirements.txt`.

| Case | CPU | GPU end to end | Speed-up |
|---|---:|---:|---:|
| Random 10 kb (k/w 15/10, 31/20, 55/10) | 31–37 ms | 3–6 ms | 4.8–11.6× |
| Random 100 kb | 331–370 ms | 30–64 ms | 5.2–12.4× |
| Random 1 Mb | 3.1–3.7 s | 0.32–0.64 s | 4.9–11.4× |
| APP reference (290 kb, k=11, w=6) | 917 ms | 328 ms | 2.8× |
| PSEN1 reference (87 kb) | 288 ms | 89 ms | 3.2× |
| PSEN2 reference (58 kb) | 181 ms | 58 ms | 3.1× |

- Parity: all 400 random cases and all 12 benchmark cases identical.
- For inputs of 100 kb or more, GPU device time is about 3–7% of the GPU path; Python record construction dominates and caps the speed-up. Copying results back is the largest device phase for every input of 100 kb or more; compaction is the largest at 10 kb.
- GPU memory pool: at most 27 MB.
- Individual speed-ups varied by up to about 20% between two runs; treat them as approximate.

Not covered: HPC GPUs, process RSS, graph construction (CPU), and whole-pipeline time. The final proposed evaluation in `outputs/proposed_final_evaluation/` was run with `--backend cpu`.
