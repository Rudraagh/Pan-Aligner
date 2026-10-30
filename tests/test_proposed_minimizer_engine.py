from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed.minimizer_engine import _compute_cpu_minimizers, _compute_gpu_minimizers, _load_usable_cupy, compute_minimizers


class MinimizerEngineTests(unittest.TestCase):
    def test_cpu_matches_known_reference_output(self) -> None:
        payload = compute_minimizers("ACGTACGTACGT", k=3, window_size=4, assembly_id="hap1", backend="cpu")
        self.assertEqual(payload["backend_used"], "cpu")
        self.assertEqual(payload["records"], [
            {"minimizer": "TAC", "hash_value": 6524611422502594526, "position": 3, "window_start": 0, "source_assembly": "hap1"},
            {"minimizer": "TAC", "hash_value": 6524611422502594526, "position": 7, "window_start": 4, "source_assembly": "hap1"},
        ])

    def test_cpu_boundaries_ambiguity_and_determinism(self) -> None:
        cases = [("", 3, 4), ("AC", 3, 4), ("ACG", 3, 4), ("ACGT", 3, 4), ("NNRYACGTN", 3, 8), ("AAAAAAAAAA", 3, 4)]
        for sequence, k, window in cases:
            with self.subTest(sequence=sequence, k=k, window=window):
                first = compute_minimizers(sequence, k, window, "hap", "cpu")
                second = compute_minimizers(sequence, k, window, "hap", "cpu")
                self.assertEqual(first["records"], second["records"])
        self.assertEqual(compute_minimizers("", 3, 4, "hap", "cpu")["records"], [])
        self.assertEqual(len(compute_minimizers("ACGT", 3, 4, "hap", "cpu")["records"]), 1)

    def test_invalid_parameters_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            compute_minimizers("ACGT", 0, 4, "hap", "cpu")
        with self.assertRaises(ValueError):
            compute_minimizers("ACGT", 3, 0, "hap", "cpu")

    def test_auto_falls_back_or_executes_gpu_truthfully(self) -> None:
        payload = compute_minimizers("ACGTACGT", 3, 4, "hap", "auto")
        self.assertIn(payload["backend_used"], {"cpu", "gpu"})
        self.assertIn("CPU" if payload["backend_used"] == "cpu" else "executed on the GPU", payload["backend_note"])

    def test_explicit_gpu_requires_a_usable_device(self) -> None:
        cp, _ = _load_usable_cupy()
        if cp is None:
            with self.assertRaisesRegex(RuntimeError, "GPU backend was requested but is unavailable"):
                compute_minimizers("ACGTACGT", 3, 4, "hap", "gpu")
        else:
            self.skipTest("GPU correctness is covered by the CUDA-gated parity test.")

    def test_gpu_matches_cpu_exactly_when_cuda_is_available(self) -> None:
        cp, reason = _load_usable_cupy()
        if cp is None:
            self.skipTest(reason)
        cases = [("ACGTACGTACGT", 3, 4), ("AAAAAAAAAA", 3, 4), ("NNRYACGTN", 3, 8), ("ACGT", 3, 4), ("AC", 3, 4)]
        for sequence, k, window in cases:
            with self.subTest(sequence=sequence, k=k, window=window):
                cpu = compute_minimizers(sequence, k, window, "hap", "cpu")
                gpu = compute_minimizers(sequence, k, window, "hap", "gpu")
                self.assertEqual(gpu["backend_used"], "gpu")
                self.assertEqual(gpu["records"], cpu["records"])

    def test_gpu_matches_cpu_on_seeded_random_cases(self) -> None:
        cp, reason = _load_usable_cupy()
        if cp is None:
            self.skipTest(reason)
        rng = random.Random(20261030)
        alphabet = "ACGT" * 8 + "NRYKMSWBDHV"
        for case in range(120):
            sequence = "".join(rng.choice(alphabet) for _ in range(rng.randint(0, 300)))
            k = rng.choice([1, 2, 3, 11, 15, 31, 54, 55, rng.randint(1, 55)])
            window = rng.randint(1, 60)
            with self.subTest(case=case, length=len(sequence), k=k, window=window):
                cpu = compute_minimizers(sequence, k, window, "hap", "cpu")
                gpu = compute_minimizers(sequence, k, window, "hap", "gpu")
                self.assertEqual(gpu["backend_used"], "gpu")
                self.assertEqual(gpu["records"], cpu["records"])

    def test_gpu_kernel_scalar_arguments_are_int32(self) -> None:
        # CuPy passes Python ints as long long; the kernels declare `int` parameters.
        cp, reason = _load_usable_cupy()
        if cp is None:
            self.skipTest(reason)
        real_kernel, seen = cp.RawKernel, []

        class RecordingKernel:
            def __init__(self, code, name, *args, **kwargs):
                self.name, self.kernel = name, real_kernel(code, name, *args, **kwargs)

            def __call__(self, grid, block, args, *rest, **kwargs):
                seen.append((self.name, [type(value) for value in args if not isinstance(value, cp.ndarray)]))
                return self.kernel(grid, block, args, *rest, **kwargs)

        with mock.patch.object(cp, "RawKernel", RecordingKernel):
            records = _compute_gpu_minimizers("ACGTACGTACGTTGCA", 3, 4, "hap", cp)
        self.assertEqual(records, _compute_cpu_minimizers("ACGTACGTACGTTGCA", 3, 4, "hap"))
        self.assertEqual([name for name, _ in seen], ["sha1_prefix64", "window_minima"])
        for name, scalar_types in seen:
            with self.subTest(kernel=name):
                self.assertTrue(scalar_types)
                self.assertTrue(all(kind is numpy.int32 for kind in scalar_types), scalar_types)

    def test_gpu_phase_timings_are_recorded_when_requested(self) -> None:
        cp, reason = _load_usable_cupy()
        if cp is None:
            self.skipTest(reason)
        timings: dict = {}
        sequence = "ACGTTGCATGACCTAGGATCCAGTTACG" * 20
        records = _compute_gpu_minimizers(sequence, 11, 6, "hap", cp, timings=timings)
        self.assertEqual(records, _compute_cpu_minimizers(sequence, 11, 6, "hap"))
        self.assertEqual(set(timings), {"host_to_device_ms", "hash_kernel_ms", "window_kernel_ms", "compaction_ms", "device_to_host_ms"})
        self.assertTrue(all(value >= 0.0 for value in timings.values()))


if __name__ == "__main__":
    unittest.main()
