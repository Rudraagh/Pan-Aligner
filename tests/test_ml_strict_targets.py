"""Strict-target evaluation of Model A: validation-only thresholds, abstention, equal-coverage comparison.

Small constructed inputs exercise the code paths; nothing is concluded from them.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from ml.build_alignment_dataset import FEATURE_COLUMNS

HAVE_ML_LIBRARIES = all(importlib.util.find_spec(name) for name in ("numpy", "sklearn"))
if HAVE_ML_LIBRARIES:
    import numpy as np
    from ml import alignment_reliability as ar
    from ml import strict_targets as st

GENES = ("APP", "PSEN1", "PSEN2")


def _rows(n_groups=18, per_group=40):
    rng = np.random.default_rng(0)
    rows = []
    for group in range(n_groups):
        for index in range(per_group):
            incorrect = index % 4 == 0
            overlap = rng.random() < 0.15  # some rows are ambiguous, so an abstention band can exist
            good = not incorrect
            row = {column: 1.0 for column in FEATURE_COLUMNS}
            row.update({
                "mapq": int(rng.integers(0, 61)) if overlap else (int(rng.integers(0, 15)) if incorrect else 60),
                "identity": 0.8 if overlap else (0.6 if incorrect else 0.95),
                "group": f"g{group}", "aligner": "minigraph", "label": "correct" if good else "incorrect", "setting": "s",
                "decoy": int(incorrect and index % 8 == 0), "own_gene": GENES[group % 3], "graph_gene": GENES[group % 3],
                "read": f"g{group}.S1_{index}!ref!{group * 25000}!{group * 25000 + 1000}!+",
            })
            rows.append(row)
    return rows


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class ValidationSplitTests(unittest.TestCase):
    def test_separate_parts_are_group_disjoint_and_exclude_the_test_fold(self) -> None:
        groups = np.array([f"w{i // 10}" for i in range(400)])
        y = np.array([i % 5 == 0 for i in range(400)], dtype=int)
        train, test = np.arange(0, 300), np.arange(300, 400)
        fit, calibration, thresholds = ar.fold_parts(train, y, groups, seed=0, separate_threshold_set=True)
        parts = [set(groups[part]) for part in (fit, calibration, thresholds, test)]
        for a in range(4):
            for b in range(a + 1, 4):
                self.assertFalse(parts[a] & parts[b], (a, b))
        self.assertEqual(len(fit) + len(calibration) + len(thresholds), len(train))

    def test_default_keeps_one_shared_part(self) -> None:
        groups = np.array([f"w{i // 10}" for i in range(200)])
        y = np.array([i % 5 == 0 for i in range(200)], dtype=int)
        fit, calibration, thresholds = ar.fold_parts(np.arange(200), y, groups, seed=0, separate_threshold_set=False)
        self.assertEqual(calibration.tolist(), thresholds.tolist())
        self.assertFalse(set(fit.tolist()) & set(calibration.tolist()))

    def test_recorded_threshold_groups_never_meet_the_test_fold(self) -> None:
        result = ar.cross_validate(_rows(), "logistic", n_splits=3, seed=0, grouping="window", separate_threshold_set=True)
        for fold in result["thresholds"]:
            self.assertFalse(set(fold["threshold_groups"]) & set(fold["test_groups"]))
            self.assertFalse(set(fold["threshold_groups"]) & set(fold["calibration_groups"]))


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class TargetTests(unittest.TestCase):
    def test_coverage_threshold_accepts_the_requested_share(self) -> None:
        prob = np.arange(1, 101) / 100.0
        self.assertEqual(ar.coverage_threshold(prob, 0.9), 0.9)

    def test_apply_targets_matches_a_direct_run(self) -> None:
        rows = _rows()
        direct = ar.cross_validate(rows, "logistic", n_splits=3, seed=0, grouping="window", separate_threshold_set=True,
                                   max_accept_error=0.05, min_reject_precision=0.95)
        other = ar.cross_validate(rows, "logistic", n_splits=3, seed=0, grouping="window", separate_threshold_set=True)
        decision, thresholds = ar.apply_targets(other, max_accept_error=0.05, min_reject_precision=0.95)
        self.assertEqual(decision.tolist(), direct["decision"].tolist())
        self.assertEqual([(t["low"], t["high"]) for t in thresholds], [(t["low"], t["high"]) for t in direct["thresholds"]])

    def test_coverage_target_replaces_the_error_target(self) -> None:
        result = ar.cross_validate(_rows(), "logistic", n_splits=3, seed=0, grouping="window", separate_threshold_set=True)
        decision, _ = ar.apply_targets(result, max_accept_error=0.0, min_reject_precision=0.99, target_coverage=0.5)
        self.assertGreater((decision == "correct").mean(), 0.2)


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class AbstentionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.Y = np.array([0, 0, 0, 0, 1, 0, 1, 1])
        self.DECISION = np.array(["correct"] * 4 + ["insufficient_evidence"] * 2 + ["incorrect"] * 2, dtype=object)
        self.OPEN = [{"low": 0.1, "high": 0.9}, {"low": 0.2, "high": 0.8}]

    def test_meaningful_when_abstained_error_lies_between_accepted_and_rejected(self) -> None:
        report = st.assess_abstention(self.Y, self.DECISION, self.OPEN)
        self.assertEqual((report["abstentions"], report["folds_with_band"]), (2, 2))
        self.assertEqual(report["abstained_error_rate"], 0.5)
        self.assertTrue(report["meaningful"])

    def test_not_meaningful_without_abstentions(self) -> None:
        decision = np.array(["correct"] * 6 + ["incorrect"] * 2, dtype=object)
        self.assertFalse(st.assess_abstention(self.Y, decision, self.OPEN)["meaningful"])

    def test_not_meaningful_when_a_fold_has_no_band(self) -> None:
        report = st.assess_abstention(self.Y, self.DECISION, [{"low": 0.5, "high": 0.5}, {"low": 0.1, "high": 0.9}])
        self.assertEqual(report["folds_with_band"], 1)
        self.assertFalse(report["meaningful"])

    def test_not_meaningful_when_abstained_rows_are_as_clean_as_accepted(self) -> None:
        y = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        self.assertFalse(st.assess_abstention(y, self.DECISION, self.OPEN)["meaningful"])


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class EqualCoverageTests(unittest.TestCase):
    def test_weights_share_a_tied_tier_evenly(self) -> None:
        weights = st.acceptance_weights(np.array([0.1, 0.5, 0.5, 0.9]), 0.5)
        self.assertEqual(weights.tolist(), [1.0, 0.5, 0.5, 0.0])

    def test_row_order_does_not_change_a_tied_comparison(self) -> None:
        prob, y = np.array([0.5, 0.5, 0.5, 0.5]), np.array([1, 0, 0, 0])
        forward = st.weighted_error(y, st.acceptance_weights(prob, 0.5))
        backward = st.weighted_error(y[::-1], st.acceptance_weights(prob, 0.5))
        self.assertEqual((forward, backward), (0.25, 0.25))

    def test_comparison_splits_within_locus_and_off_target(self) -> None:
        y = np.array([0, 0, 1, 1, 0, 0, 1, 0])
        prob = np.array([0.1, 0.2, 0.9, 0.8, 0.1, 0.2, 0.3, 0.2])
        mapq = np.array([60, 60, 60, 0, 60, 60, 0, 60])
        decoy = np.array([0, 0, 1, 0, 0, 0, 0, 0])
        groups = np.array(["a", "a", "b", "b", "c", "c", "d", "d"])
        row = st.equal_coverage_comparison(y, prob, mapq, decoy, groups, thresholds=(60,), resamples=50, seed=0)[0]
        self.assertEqual(row["coverage"], 0.75)
        self.assertAlmostEqual(row["mapq_error"], 1 / 6)
        self.assertEqual(row["mapq_accepted_off_target"], 1)
        self.assertEqual(row["mapq_accepted_within_locus"], 0)
        self.assertAlmostEqual(row["model_error"], 1 / 6)  # only five rows are correct, so a sixth must be an error
        self.assertEqual(row["model_accepted_within_locus"], 1.0)
        self.assertEqual(len(row["difference_ci95"]), 2)


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class StudyTests(unittest.TestCase):
    def test_targets_are_fixed_before_results_and_strictly_tighter_than_defaults(self) -> None:
        names = {target["name"] for target in st.TARGETS}
        self.assertTrue(any(target.get("target_coverage") for target in st.TARGETS))
        self.assertEqual(len(names), len(st.TARGETS))
        for target in st.TARGETS:
            self.assertLessEqual(target["max_accept_error"], 0.01)
            self.assertGreaterEqual(target["min_reject_precision"], 0.9)

    def test_run_writes_a_report_without_claiming_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "window"
            results = st.run_grouping(_rows(), "window", output, n_splits=3, resamples=20, seed=0)
            self.assertEqual({row["target"] for row in results["targets"]}, {target["name"] for target in st.TARGETS})
            for name in ("REPORT.md", "strict_results.json", "risk_coverage_split.png", "reliability.png"):
                self.assertTrue((output / name).exists(), name)
            text = (output / "REPORT.md").read_text(encoding="utf-8")
        for expected in ("Within-locus", "Off-target", "Abstentions", "validation", "MAPQ"):
            self.assertIn(expected, text)
        self.assertNotIn("replaces MAPQ", text)
        self.assertIn(results["mapq_comparison_verdict"], text)

    def test_verdict_requires_the_interval_to_favour_the_model_everywhere(self) -> None:
        win = [{"threshold": 30, "difference_ci95": (-0.01, -0.001)}, {"threshold": 60, "difference_ci95": (-0.02, -0.01)}]
        mixed = [{"threshold": 30, "difference_ci95": (-0.01, -0.001)}, {"threshold": 60, "difference_ci95": (-0.001, 0.001)}]
        self.assertIn("at every", st.mapq_verdict(win))
        self.assertIn("not at every", st.mapq_verdict(mixed))

    def test_existing_results_are_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "window"
            output.mkdir()
            (output / "strict_results.json").write_text("{}", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                st.run_grouping(_rows(), "window", output, n_splits=3, resamples=5, seed=0)

    def test_default_output_is_separate_from_existing_results(self) -> None:
        self.assertEqual(st.DEFAULT_OUTPUT.name, "strict_targets")
        self.assertEqual(st.DEFAULT_OUTPUT.parent.name, "ml_alignment_v2")


if __name__ == "__main__":
    unittest.main()
