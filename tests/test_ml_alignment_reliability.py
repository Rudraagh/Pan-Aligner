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


GENES = ("APP", "PSEN1", "PSEN2")


def _rows(n_groups=10, per_group=30, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for group in range(n_groups):
        for index in range(per_group):
            incorrect = index % 5 == 0
            mapq = int(rng.integers(0, 10)) if incorrect else int(rng.integers(40, 61))
            row = {column: 1.0 for column in FEATURE_COLUMNS}
            row.update({"mapq": mapq, "identity": 0.7 if incorrect else 0.95, "group": f"g{group}", "aligner": "panaligner",
                        "label": "incorrect" if incorrect else "correct", "setting": "s", "decoy": 0,
                        "own_gene": GENES[group % 3], "graph_gene": GENES[group % 3],
                        "read": f"g{group}.S1_{index}!ref!{group * 25000 + index * 100}!{group * 25000 + index * 100 + 1000}!+"})
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
class GroupingTests(unittest.TestCase):
    def _labelled(self):
        rows = [row for row in _rows() if row["label"] in ("correct", "incorrect")]
        rows[0] = {**rows[0], "decoy": 1, "graph_gene": "PSEN2", "label": "incorrect"}  # APP read on the PSEN2 graph
        return rows

    def test_window_folds_keep_windows_apart_and_purge_neighbours(self) -> None:
        rows = self._labelled()
        folds, _ = ar.evaluation_folds(rows, "window", n_splits=3)
        key = [(row["own_gene"], ar.window_index(row)) for row in rows]
        tested = []
        for train, test in folds:
            tested.extend(test.tolist())
            test_keys = {key[i] for i in test}
            for i in train:
                gene, index = key[i]
                self.assertFalse({(gene, index - 1), (gene, index), (gene, index + 1)} & test_keys)
        self.assertEqual(sorted(tested), list(range(len(rows))))

    def test_leave_one_gene_out_removes_the_gene_and_its_graph_from_training(self) -> None:
        rows = self._labelled()
        folds, _ = ar.evaluation_folds(rows, "leave_one_gene_out", n_splits=3)
        self.assertEqual(len(folds), 3)
        tested = []
        for train, test in folds:
            tested.extend(test.tolist())
            genes = {rows[i]["own_gene"] for i in test}
            self.assertEqual(len(genes), 1)
            gene = genes.pop()
            self.assertFalse(any(rows[i]["own_gene"] == gene or rows[i]["graph_gene"] == gene for i in train))
        self.assertEqual(sorted(tested), list(range(len(rows))))

    def test_cross_validate_accepts_each_grouping(self) -> None:
        for grouping in ("haplotype", "window", "leave_one_gene_out"):
            with self.subTest(grouping=grouping):
                result = ar.cross_validate(_rows(n_groups=12), "logistic", n_splits=3, seed=0, grouping=grouping)
                self.assertEqual(result["grouping"], grouping)
                self.assertEqual(len(result["prob_incorrect"]), len(result["y"]))

    def test_default_output_never_overwrites_the_committed_results(self) -> None:
        path = ar.output_dir_for(None, "window")
        self.assertEqual(path, ar.DEFAULT_DIR / "regrouped" / "window")
        self.assertEqual(ar.output_dir_for(Path("x"), "window"), Path("x"))


@unittest.skipIf(ar is None, "scikit-learn/numpy not installed in this environment")
class MapqThresholdTests(unittest.TestCase):
    def test_threshold_sweep_known_answers(self) -> None:
        y, mapq = np.array([0, 0, 1, 1]), np.array([60, 30, 60, 0])
        sweep = {row["threshold"]: row for row in ar.mapq_threshold_sweep(y, mapq, thresholds=(0, 30, 60))}
        self.assertEqual((sweep[0]["coverage"], sweep[0]["accepted_error_rate"], sweep[0]["incorrect_caught"]), (1.0, 0.5, 0.0))
        self.assertAlmostEqual(sweep[30]["accepted_error_rate"], 1 / 3)
        self.assertEqual((sweep[30]["coverage"], sweep[30]["incorrect_caught"]), (0.75, 0.5))
        self.assertEqual((sweep[60]["coverage"], sweep[60]["accepted_error_rate"]), (0.5, 0.5))

    def test_model_error_at_matching_coverage(self) -> None:
        y, prob = np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.9, 0.3])
        self.assertEqual(ar.model_error_at_coverage(y, prob, 0.5), 0.0)
        self.assertAlmostEqual(ar.model_error_at_coverage(y, prob, 0.75), 1 / 3)

    def test_summary_separates_within_locus_from_off_target(self) -> None:
        rows = _rows()
        rows[0] = {**rows[0], "decoy": 1}
        overall = ar.summarize(ar.cross_validate(rows, "logistic", n_splits=5, seed=0))["overall"]
        self.assertEqual(overall["off_target"]["incorrect_rows"], 1)
        self.assertEqual(overall["within_locus"]["incorrect_rows"], sum(r["label"] == "incorrect" for r in rows) - 1)
        self.assertIn("pr_auc_incorrect", overall["within_locus"])


@unittest.skipIf(ar is None, "scikit-learn/numpy not installed in this environment")
class ReportLimitsTests(unittest.TestCase):
    def _limits(self, decoys: int, within: int, within_aligner: str = "minigraph") -> str:
        base = {"group": "g", "setting": "s", "own_gene": "APP", "graph_gene": "APP", "read": "r!x!1!2!+"}
        rows = [{**base, "label": "incorrect", "decoy": 1, "aligner": "minigraph"}] * decoys
        rows += [{**base, "label": "incorrect", "decoy": 0, "aligner": within_aligner}] * within
        return ar.report_limits(rows)

    def test_dominated_wording_only_when_off_target_really_dominates(self) -> None:
        self.assertIn("dominated by off-target", self._limits(99, 1)[0])
        text = self._limits(66, 34)[0]
        self.assertNotIn("dominated", text)
        self.assertIn("66%", text)

    def test_names_the_aligner_share_of_within_locus_errors(self) -> None:
        self.assertIn("100% from minigraph", " ".join(self._limits(66, 34)))


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
