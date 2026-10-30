"""Mechanics of Model A training: grouping, calibration metrics, abstention.

The arrays here are small constructed inputs that exercise the code paths; they are
not data and nothing is concluded from them.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import importlib.util

from ml.build_alignment_dataset import FEATURE_COLUMNS

# scikit-learn is only installed where requirements-ml.txt is; a missing module of our own must still fail.
HAVE_ML_LIBRARIES = all(importlib.util.find_spec(name) for name in ("numpy", "sklearn"))
if HAVE_ML_LIBRARIES:
    import numpy as np
    from ml import alignment_reliability as ar
else:
    ar = None


def _rows(n_groups=10, per_group=30, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for group in range(n_groups):
        for index in range(per_group):
            incorrect = index % 5 == 0
            mapq = int(rng.integers(0, 10)) if incorrect else int(rng.integers(40, 61))
            row = {column: 1.0 for column in FEATURE_COLUMNS}
            row.update({"mapq": mapq, "identity": 0.7 if incorrect else 0.95, "group": f"g{group}", "aligner": "panaligner",
                        "label": "incorrect" if incorrect else "correct", "setting": "s", "decoy": 0})
            rows.append(row)
    rows.append({**rows[0], "label": "unaligned"})
    rows.append({**rows[1], "label": "unprojectable"})
    return rows


@unittest.skipIf(ar is None, "scikit-learn/numpy not installed in this environment")
class MetricTests(unittest.TestCase):
    def test_expected_calibration_error_known_answers(self) -> None:
        y = np.array([1] * 8 + [0] * 2)
        self.assertAlmostEqual(ar.expected_calibration_error(y, np.full(10, 0.8)), 0.0)
        self.assertAlmostEqual(ar.expected_calibration_error(np.array([1, 0] * 5), np.ones(10)), 0.5)

    def test_thresholds_meet_targets_on_the_calibration_set(self) -> None:
        prob = np.array([0.01, 0.02, 0.03, 0.05, 0.4, 0.5, 0.92, 0.95, 0.97, 0.99])
        y = np.array([0, 0, 0, 0, 1, 0, 1, 1, 1, 1])
        low, high = ar.choose_thresholds(prob, y, max_accept_error=0.0, min_reject_precision=1.0)
        self.assertEqual((low, high), (0.05, 0.92))
        decisions = ar.three_way(prob, low, high)
        self.assertEqual(list(decisions), ["correct"] * 4 + ["insufficient_evidence"] * 2 + ["incorrect"] * 4)

    def test_unreachable_targets_abstain_instead_of_guessing(self) -> None:
        prob, y = np.array([0.1, 0.2, 0.8, 0.9]), np.array([1, 0, 1, 0])
        low, high = ar.choose_thresholds(prob, y, max_accept_error=0.0, min_reject_precision=1.0)
        self.assertEqual(set(ar.three_way(prob, low, high)), {"insufficient_evidence"})


@unittest.skipIf(ar is None, "scikit-learn/numpy not installed in this environment")
class CrossValidationTests(unittest.TestCase):
    def test_outer_folds_never_share_a_group(self) -> None:
        groups = np.array([f"g{i // 7}" for i in range(70)])
        for train, test in ar.group_folds(groups, n_splits=5):
            self.assertFalse(set(groups[train]) & set(groups[test]))

    def test_only_labelled_rows_are_used_and_every_one_gets_a_decision(self) -> None:
        rows = _rows()
        result = ar.cross_validate(rows, "logistic", n_splits=5, seed=0)
        labelled = [row for row in rows if row["label"] in ("correct", "incorrect")]
        self.assertEqual(len(result["prob_incorrect"]), len(labelled))
        self.assertLessEqual(set(result["decision"]), {"correct", "incorrect", "insufficient_evidence"})
        self.assertTrue(np.all((result["prob_incorrect"] >= 0) & (result["prob_incorrect"] <= 1)))

    def test_feature_matrix_excludes_simulation_metadata(self) -> None:
        names = ar.model_feature_names()
        for leaked in ("accuracy", "read_length", "setting", "decoy", "graph_gene", "own_gene", "seed", "group"):
            self.assertNotIn(leaked, names)
        self.assertIn("mapq", names)
        self.assertIn("aligner_panaligner", names)

    def test_mapq_baseline_uses_mapq_only(self) -> None:
        self.assertEqual(ar.MODEL_FEATURES["mapq_only"], ["mapq"])


@unittest.skipIf(ar is None, "scikit-learn/numpy not installed in this environment")
class ReportTests(unittest.TestCase):
    def test_report_states_numbers_and_caveats_from_the_results(self) -> None:
        import tempfile

        rows = _rows()
        for row in rows[:12]:
            row["decoy"] = 1 if row["label"] == "incorrect" else 0
        summaries = [ar.summarize(ar.cross_validate(rows, name, n_splits=5, seed=0)) for name in ("mapq_only", "logistic")]
        summaries[1]["thresholds"] = [{"fold": 0, "low": 0.4, "high": 0.4}, {"fold": 1, "low": 0.1, "high": 0.9}]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "REPORT.md"
            ar.write_report(path, rows, summaries, {"max_accept_error": 0.01, "min_reject_precision": 0.9})
            text = path.read_text(encoding="utf-8")
        for expected in ("mapq_only", "logistic", f"{summaries[1]['overall']['pr_auc_incorrect']:.3f}",
                         "simulation truth", "abstention band collapsed in 1 of 2 folds", "off-target"):
            with self.subTest(expected=expected):
                self.assertIn(expected, text)
        labelled = [row for row in rows if row["label"] == "incorrect"]
        self.assertIn(f"{sum(int(row['decoy']) == 1 for row in labelled)} of {len(labelled)} incorrect", text)


if __name__ == "__main__":
    unittest.main()
