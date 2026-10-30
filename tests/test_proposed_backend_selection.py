from __future__ import annotations

import sys
import unittest
from dataclasses import asdict
from unittest import mock

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from proposed import graph_constructor, minimizer_engine
from proposed.final_evaluation import _backend_statement, run_final_proposed_evaluation
from proposed.graph_constructor import build_proposed_pangenome_graph
from proposed.models import AssemblyInput


ASSEMBLIES = [
    AssemblyInput("h1", "PSEN2", "HEALTHY", "ACGTACGTTTGACCATGCAAGTCCGATAGCTAGGCT"),
    AssemblyInput("h2", "PSEN2", "HEALTHY", "ACGTACGTTTGACCATGCTAGTCCGATAGCTAGGCT"),
    AssemblyInput("u1", "PSEN2", "UNHEALTHY", "ACGTACGTTTGACGATGCAAGTCCGATAGCTAGGCA"),
]


def _graph_signature(graph):
    return (
        [asdict(node) for node in graph.nodes],
        [asdict(edge) for edge in graph.edges],
        [asdict(path) for path in graph.paths],
        graph.quality,
    )


class BackendSelectionTests(unittest.TestCase):
    def test_cpu_backend_is_recorded_per_stage(self) -> None:
        graph = build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="cpu")
        self.assertEqual(graph.backend_used, "cpu")
        backend = graph.metadata["backend"]
        self.assertEqual((backend["requested"], backend["minimizer_extraction"], backend["graph_construction"]), ("cpu", "cpu", "cpu"))

    def test_explicit_gpu_fails_clearly_when_cuda_is_unusable(self) -> None:
        with mock.patch.object(minimizer_engine, "_load_usable_cupy", return_value=(None, "No CUDA device (simulated).")):
            with self.assertRaisesRegex(RuntimeError, "GPU backend was requested but is unavailable"):
                build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="gpu")

    def test_auto_fallback_is_recorded_not_hidden(self) -> None:
        with mock.patch.object(minimizer_engine, "_load_usable_cupy", return_value=(None, "No CUDA device (simulated).")):
            graph = build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="auto")
        self.assertEqual(graph.backend_used, "cpu")
        self.assertEqual(graph.metadata["backend"]["requested"], "auto")
        self.assertTrue(any("Auto mode selected CPU" in note for note in graph.metadata["backend"]["notes"]))

    def test_partial_gpu_fallback_is_reported_as_mixed(self) -> None:
        real = minimizer_engine.compute_minimizers

        def fake(sequence, k, window_size, assembly_id, backend="auto"):
            payload = real(sequence, k, window_size, assembly_id, "cpu")
            if assembly_id == "h1":
                payload["backend_used"], payload["backend_note"] = "gpu", "simulated GPU run"
            return payload

        with mock.patch.object(graph_constructor, "compute_minimizers", side_effect=fake):
            graph = build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="auto")
        self.assertEqual(graph.backend_used, "mixed")
        self.assertEqual(graph.metadata["backend"]["graph_construction"], "cpu")

    def test_invalid_backend_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="tpu")
        with self.assertRaises(ValueError):
            run_final_proposed_evaluation(backend="tpu")

    def test_readme_statement_never_claims_gpu_for_cpu_runs(self) -> None:
        self.assertIn("no GPU acceleration", _backend_statement("auto", [{"backend_used": "cpu"}]))
        statement = _backend_statement("gpu", [{"backend_used": "gpu"}])
        self.assertIn("Only minimizer extraction can run on the GPU", statement)

    def test_gpu_graph_matches_cpu_graph_when_cuda_is_available(self) -> None:
        cp, reason = minimizer_engine._load_usable_cupy()
        if cp is None:
            self.skipTest(reason)
        cpu = build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="cpu")
        gpu = build_proposed_pangenome_graph(ASSEMBLIES, k=5, window_size=3, backend="gpu")
        self.assertEqual(gpu.backend_used, "gpu")
        self.assertEqual(_graph_signature(gpu), _graph_signature(cpu))


if __name__ == "__main__":
    unittest.main()
