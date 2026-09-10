import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import networkx as nx


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from predictor import generate_custom_alignment_overlay, has_valid_alignment, plot_custom_query_scores
from parse_gaf import AlignmentResult
from visualize import alignment_highlighted_nodes


class FakeAxis:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def record(*args, **kwargs):
            self.calls.append((name, args, kwargs))
        return record


class FakeFigure:
    def __init__(self):
        self.calls = []

    def suptitle(self, *args, **kwargs):
        self.calls.append(("suptitle", args, kwargs))

    def tight_layout(self):
        self.calls.append(("tight_layout", (), {}))

    def savefig(self, path, **kwargs):
        self.calls.append(("savefig", (path,), kwargs))


def unavailable_prediction() -> dict:
    unavailable = {"available": False, "identity": 0.0, "coverage": 0.0}
    zero_score = {"composite_score": 0.0}
    return {
        "selected_gene": "PSEN2",
        "alignment_detected": False,
        "combined_alignment": unavailable,
        "class_alignments": {"HEALTHY": unavailable, "UNHEALTHY": unavailable},
        "score_breakdown": {"HEALTHY": zero_score, "UNHEALTHY": zero_score, "COMBINED": zero_score},
        "gene_detection": {},
    }


def valid_prediction() -> dict:
    available = {"available": True, "identity": 1.0, "coverage": 0.99}
    return {
        "selected_gene": "PSEN2",
        "alignment_detected": True,
        "combined_alignment": available,
        "class_alignments": {"HEALTHY": available, "UNHEALTHY": available},
        "score_breakdown": {
            "HEALTHY": {"composite_score": 0.9},
            "UNHEALTHY": {"composite_score": 0.91},
            "COMBINED": {"composite_score": 0.91},
        },
        "gene_detection": {},
    }


class CustomQueryPresentationTests(unittest.TestCase):
    def _plot_with_fakes(self, prediction: dict):
        figure = FakeFigure()
        axes = [FakeAxis(), FakeAxis()]
        with patch("predictor.plt.subplots", return_value=(figure, axes)), patch("predictor.plt.close"):
            path = plot_custom_query_scores(prediction, ROOT)
        self.assertEqual(path, ROOT / "custom_query_scores.png")
        return figure, axes

    def test_no_alignment_renders_unavailable_status_not_zero_measurements(self):
        figure, axes = self._plot_with_fakes(unavailable_prediction())
        self.assertFalse(has_valid_alignment(unavailable_prediction()))
        self.assertFalse(any(name in {"bar", "plot"} for axis in axes for name, _, _ in axis.calls))
        text = " ".join(str(args[2]) for axis in axes for name, args, _ in axis.calls if name == "text")
        self.assertIn("No valid alignment detected", text)
        self.assertIn("Scores not applicable", text)

    def test_valid_alignment_keeps_score_plotting(self):
        _, axes = self._plot_with_fakes(valid_prediction())
        self.assertTrue(has_valid_alignment(valid_prediction()))
        bar_call = next(call for call in axes[0].calls if call[0] == "bar")
        self.assertEqual(bar_call[1][1], [0.9, 0.91, 0.91])
        score_labels = [args[2] for name, args, _ in axes[0].calls if name == "text"]
        self.assertEqual(score_labels, ["0.900000", "0.910000", "0.910000"])
        self.assertEqual(sum(call[0] == "plot" for call in axes[1].calls), 2)

    def test_existing_real_psen2_output_remains_compatible(self):
        prediction_path = ROOT / "outputs" / "alignments" / "custom_query_real_psen2_demo" / "prediction.json"
        if not prediction_path.exists():
            self.skipTest("Real PSEN2 demonstration output is not present.")
        import json

        prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
        self.assertTrue(has_valid_alignment(prediction))
        self.assertEqual(prediction["prediction"], "UNHEALTHY")
        self.assertIn("healthy_score", prediction)
        self.assertIn("unhealthy_score", prediction)
        self.assertIn("score_breakdown", prediction)

    def test_no_combined_alignment_skips_overlay(self):
        prediction = unavailable_prediction()
        result = generate_custom_alignment_overlay(prediction, ROOT / "data" / "metadata" / "graph_manifest.json", ROOT / "tests" / "_tmp")
        self.assertIsNone(result)
        self.assertEqual(prediction["alignment_overlay_status"], "SKIPPED_NO_COMBINED_ALIGNMENT")

    def test_overlay_failure_is_non_fatal(self):
        output_dir = ROOT / "outputs" / "alignments" / "custom_query_real_psen2_presentation_test"
        if not (output_dir / "psen2.combined.gaf").exists():
            self.skipTest("Real PSEN2 custom-query GAF is not present.")
        prediction = valid_prediction()
        with patch("predictor.visualize_alignment", side_effect=RuntimeError("render error")):
            result = generate_custom_alignment_overlay(prediction, ROOT / "data" / "metadata" / "graph_manifest.json", output_dir)
        self.assertIsNone(result)
        self.assertEqual(prediction["alignment_overlay_status"], "FAILED")
        self.assertEqual(prediction["alignment_overlay_error"], "render error")

    def test_gaf_reference_ranges_highlight_minigraph_segment_ids(self):
        graph = nx.DiGraph()
        graph.add_node("s1", length=10, reference_name="ref|NC_000001.11|:region", reference_offset=0)
        graph.add_node("s2", length=10, reference_name="ref|NC_000001.11|:region", reference_offset=10)
        result = AlignmentResult("query", 20, 0, 20, "+", ">ref|NC_000001.11|:region:0-20", 20, 0, 20, 20, 20, 60)
        self.assertEqual(alignment_highlighted_nodes(graph, [result]), {"s1", "s2", "ref|NC_000001.11|:region:0-20"})
