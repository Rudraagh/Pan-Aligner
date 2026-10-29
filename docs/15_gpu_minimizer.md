# GPU minimizer engine

`scripts/proposed/minimizer_engine.py` retains the CPU implementation as the reference and optionally uses CuPy with CUDA for minimizer extraction.

## What runs on the GPU

For `backend="gpu"`, CUDA kernels compute the SHA-1-prefix key for every k-mer and the minimum key in every minimizer window. Consecutive repeated selected positions are compacted on the device. Only the compact minimizer positions and keys are copied to Python to create the existing `MinimizerRecord` objects. Graph construction remains CPU-side.

The GPU key is exactly the first 64 bits of the SHA-1 digest used by the CPU reference. Equal keys retain the earliest k-mer position, matching the CPU `(hash, position, kmer)` ordering. The current CUDA SHA-1 kernel supports ASCII sequences and `k <= 55`; this is sufficient for the prototype defaults. `auto` uses CPU if CUDA is unavailable or this GPU kernel cannot execute. Explicit `gpu` mode raises a clear error instead of silently using CPU.

## Installation and use

CPU-only installation remains `pip install -r requirements.txt`. For CUDA, install the CuPy package matching the installed CUDA runtime, for example the package documented at <https://docs.cupy.dev/en/stable/install.html>. CuPy deliberately remains optional because its package name depends on the user’s CUDA version.

Run with `python proposed_methodology.py --backend cpu` or `python proposed_methodology.py --backend gpu`. `--backend auto` attempts a usable CUDA device and otherwise records an explicit CPU fallback reason in the minimizer payload.

Run the focused tests with:

```powershell
python -m unittest tests/test_proposed_minimizer_engine.py -v
```

The CUDA parity test is skipped with its device/CuPy reason when CUDA is unavailable. It compares ordered records, including hashes, positions, and window starts. No GPU speed claim should be made unless this test passes on the target machine and a separate benchmark is run. GPU minimizer execution does not accelerate graph construction, scoring, anchor selection, or the baseline PanAligner workflow.
