"""Entry-point wiring: main.py forwarding modes, CPU backend defaults, CLI flag handling."""

from __future__ import annotations

import inspect
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "scripts"
for path in (ROOT, SCRIPTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import main as project_main
import proposed_methodology
import shared_query_eval
import simulate_truth_reads
import truth_eval
from proposed import final_evaluation, graph_constructor, minimizer_engine, quality_feedback, workflow


class StopRun(Exception):
    """Raised by mocks to stop a run once its arguments have been captured."""


class MainForwardingTests(unittest.TestCase):
    def test_each_mode_forwards_the_remaining_arguments(self) -> None:
        cases = {
            "--proposed": (proposed_methodology, ["--final-proposed-evaluation", "--gene", "APP", "--output-dir", "x"]),
            "--shared-query-eval": (shared_query_eval, ["--genes", "PSEN2", "--k", "11"]),
            "--simulate-truth-reads": (simulate_truth_reads, ["--source", "heldout", "--seeds", "1"]),
            "--truth-eval": (truth_eval, ["--source", "heldout", "--threads", "2"]),
        }
        for mode, (module, rest) in cases.items():
            with self.subTest(mode=mode), mock.patch.object(project_main.sys, "platform", "linux"), \
                    mock.patch.object(module, "main") as target:
                project_main.main([mode, *rest])
                target.assert_called_once_with(rest)

    def test_forwarding_mode_must_come_first(self) -> None:
        with mock.patch.object(proposed_methodology, "main") as target, self.assertRaises(SystemExit):
            project_main.main(["--threads", "2", "--proposed"])
        target.assert_not_called()

    def test_truth_modes_stop_on_windows_with_a_wsl_message(self) -> None:
        for mode, module in (("--simulate-truth-reads", simulate_truth_reads), ("--truth-eval", truth_eval)):
            with self.subTest(mode=mode), mock.patch.object(project_main.sys, "platform", "win32"), \
                    mock.patch.object(module, "main") as target:
                with self.assertRaisesRegex(SystemExit, "WSL"):
                    project_main.main([mode])
                target.assert_not_called()

    def test_ml_modes_forward_to_the_ml_tools(self) -> None:
        import importlib
        import importlib.util

        cases = {
            "--ml-alignment-data": ("ml.build_alignment_dataset", ["--seeds", "1"], False),
            "--ml-alignment": ("ml.alignment_reliability", ["--folds", "3"], True),
            "--ml-variant": ("ml.variant_feasibility", ["--seeds", "0"], True),
        }
        for mode, (module_name, rest, needs_sklearn) in cases.items():
            with self.subTest(mode=mode):
                if needs_sklearn and importlib.util.find_spec("sklearn") is None:
                    self.skipTest("scikit-learn not installed in this environment")
                module = importlib.import_module(module_name)
                with mock.patch.object(project_main.sys, "platform", "linux"), mock.patch.object(module, "main") as target:
                    project_main.main([mode, *rest])
                target.assert_called_once_with(rest)

    def test_ml_dataset_mode_stops_on_windows(self) -> None:
        with mock.patch.object(project_main.sys, "platform", "win32"), self.assertRaisesRegex(SystemExit, "WSL"):
            project_main.main(["--ml-alignment-data"])

    def test_linux_modes_run_without_the_baseline_dependencies(self) -> None:
        # WSL's Python has no networkx/Biopython; forwarding modes must not import the baseline stack first.
        import subprocess

        script = (
            "import importlib.abc, runpy, sys\n"
            "class Block(importlib.abc.MetaPathFinder):\n"
            "    def find_spec(self, name, path, target=None):\n"
            "        if name.split('.')[0] in ('networkx', 'Bio'):\n"
            "            raise ModuleNotFoundError(f'No module named {name!r} (blocked)', name=name)\n"
            "sys.meta_path.insert(0, Block())\n"
            "sys.platform = 'linux'\n"
            "sys.argv = ['main.py', sys.argv[1], '--help']\n"
            "runpy.run_path('main.py', run_name='__main__')\n"
        )
        for mode in ("--truth-eval", "--simulate-truth-reads", "--ml-alignment-data"):
            with self.subTest(mode=mode):
                completed = subprocess.run([sys.executable, "-B", "-c", script, mode], cwd=ROOT, capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr[-800:])
                self.assertIn("usage:", completed.stdout)

    def test_tool_mains_accept_argv(self) -> None:
        for module in (proposed_methodology, shared_query_eval, simulate_truth_reads, truth_eval):
            with self.subTest(module=module.__name__):
                self.assertIn("argv", inspect.signature(module.main).parameters)


class BackendDefaultTests(unittest.TestCase):
    def test_library_backend_defaults_are_cpu(self) -> None:
        for function in (
            minimizer_engine.compute_minimizers,
            graph_constructor.build_proposed_pangenome_graph,
            quality_feedback.run_pg_scunk_feedback_loop,
            workflow.run_proposed_methodology,
            final_evaluation.run_final_proposed_evaluation,
        ):
            with self.subTest(function=function.__name__):
                self.assertEqual(inspect.signature(function).parameters["backend"].default, "cpu")


class ProposedCliTests(unittest.TestCase):
    def test_final_evaluation_receives_every_passed_flag(self) -> None:
        with mock.patch.object(proposed_methodology, "run_final_proposed_evaluation", side_effect=StopRun) as target:
            with self.assertRaises(StopRun):
                proposed_methodology.main([
                    "--final-proposed-evaluation", "--backend", "gpu", "--gene", "app",
                    "--dataset-manifest", "tests/data/proposed/toy_manifest.json",
                    "--quality-threshold", "0.5", "--max-rounds", "4", "--improvement-tolerance", "0.01",
                    "--stagnation-patience", "3", "--anchor-alpha", "0.5", "--anchor-beta", "0.3", "--anchor-gamma", "0.2",
                ])
        kwargs = target.call_args.kwargs
        self.assertEqual(kwargs["backend"], "gpu")
        self.assertEqual(kwargs["genes"], ("APP",))
        self.assertEqual(kwargs["dataset_manifest"], (ROOT / "tests/data/proposed/toy_manifest.json").resolve())
        self.assertEqual(
            (kwargs["quality_threshold"], kwargs["max_rounds"], kwargs["improvement_tolerance"], kwargs["stagnation_patience"]),
            (0.5, 4, 0.01, 3),
        )
        self.assertEqual((kwargs["anchor_alpha"], kwargs["anchor_beta"], kwargs["anchor_gamma"]), (0.5, 0.3, 0.2))

    def test_final_evaluation_defaults_cover_all_three_genes(self) -> None:
        with mock.patch.object(proposed_methodology, "run_final_proposed_evaluation", side_effect=StopRun) as target:
            with self.assertRaises(StopRun):
                proposed_methodology.main(["--final-proposed-evaluation"])
        kwargs = target.call_args.kwargs
        self.assertEqual(kwargs["genes"], ("APP", "PSEN1", "PSEN2"))
        self.assertEqual(kwargs["backend"], "cpu")

    def test_final_evaluation_rejects_an_assembly_cap(self) -> None:
        with mock.patch.object(proposed_methodology, "run_final_proposed_evaluation") as target, \
                mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            proposed_methodology.main(["--final-proposed-evaluation", "--max-assemblies-per-bucket", "3"])
        target.assert_not_called()

    def test_normal_run_keeps_the_default_cap_of_eight(self) -> None:
        with mock.patch.object(proposed_methodology, "run_proposed_methodology", side_effect=StopRun) as target:
            with self.assertRaises(StopRun):
                proposed_methodology.main([])
        self.assertEqual(target.call_args.kwargs["max_assemblies_per_bucket"], 8)

    def test_stress_demo_honours_backend(self) -> None:
        with mock.patch.object(proposed_methodology, "run_pg_scunk_stress_experiment", side_effect=StopRun) as target:
            with self.assertRaises(StopRun):
                proposed_methodology.main(["--pg-scunk-stress-demo", "--backend", "auto"])
        self.assertEqual(target.call_args.kwargs["backend"], "auto")


class FinalEvaluationParameterTests(unittest.TestCase):
    def test_flags_reach_the_feedback_loop_with_the_chosen_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(final_evaluation, "run_pg_scunk_feedback_loop", side_effect=StopRun) as loop:
            with self.assertRaises(StopRun):
                final_evaluation.run_final_proposed_evaluation(
                    output_dir=Path(tmp), backend="cpu", genes=("TOY",),
                    dataset_manifest=ROOT / "tests" / "data" / "proposed" / "toy_manifest.json",
                    quality_threshold=0.5, max_rounds=4, improvement_tolerance=0.01, stagnation_patience=3,
                    anchor_alpha=0.5, anchor_beta=0.3, anchor_gamma=0.2,
                )
        assemblies = loop.call_args.args[0]
        kwargs = loop.call_args.kwargs
        self.assertEqual({item.gene for item in assemblies}, {"TOY"})
        self.assertEqual(
            (kwargs["backend"], kwargs["quality_threshold"], kwargs["max_rounds"], kwargs["improvement_tolerance"], kwargs["stagnation_patience"]),
            ("cpu", 0.5, 4, 0.01, 3),
        )
        self.assertEqual(kwargs["base_parameters"], {"anchor_alpha": 0.5, "anchor_beta": 0.3, "anchor_gamma": 0.2})


if __name__ == "__main__":
    unittest.main()
