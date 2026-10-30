"""Model B feasibility study: variant extraction, grouping, baselines, abstention."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

HAVE_ML_LIBRARIES = all(importlib.util.find_spec(name) for name in ("numpy", "sklearn"))
if HAVE_ML_LIBRARIES:
    import numpy as np
    from ml import variant_feasibility as vf


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class VariantExtractionTests(unittest.TestCase):
    REF = "ACGTACGTACCGTTAGCAGT"

    def test_snv_insertion_deletion_and_multi_snv(self) -> None:
        snv = vf.extract_variant(self.REF, self.REF[:5] + "T" + self.REF[6:])
        self.assertEqual((snv["position"], snv["ref"], snv["alt"], snv["kind"], snv["length_change"]), (5, "C", "T", "SNV", 0))
        insertion = vf.extract_variant(self.REF, self.REF[:8] + "GGG" + self.REF[8:])
        self.assertEqual((insertion["kind"], insertion["length_change"]), ("insertion", 3))
        deletion = vf.extract_variant(self.REF, self.REF[:4] + self.REF[7:])
        self.assertEqual((deletion["kind"], deletion["length_change"]), ("deletion", -3))
        double = vf.extract_variant(self.REF, "T" + self.REF[1:10] + "A" + self.REF[11:])
        self.assertEqual((double["kind"], double["position"], double["snv_count"]), ("multi_SNV", 0, 2))

    def test_transition_and_cpg_context(self) -> None:
        self.assertEqual((vf.is_transition("C", "T"), vf.is_transition("A", "G"), vf.is_transition("C", "A")), (1, 1, 0))
        self.assertEqual(vf.is_transition("CT", "C"), 0)
        reference = "AACGTTTTTT"
        self.assertEqual((vf.in_cpg(reference, 2), vf.in_cpg(reference, 3), vf.in_cpg(reference, 6)), (1, 1, 0))
        self.assertAlmostEqual(vf.gc_fraction("AACGTTTTTT", 4, flank=2), 2 / 5)

    def test_window_groups(self) -> None:
        self.assertEqual(vf.window_group("APP", 11228), "APP:11")
        self.assertEqual(vf.window_group("APP", 11999), "APP:11")
        self.assertNotEqual(vf.window_group("APP", 12000), vf.window_group("APP", 11999))


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class StudyMechanicsTests(unittest.TestCase):
    def _table(self):
        rows = []
        for index in range(60):
            gene = "APP" if index < 40 else "PSEN2"
            label = "UNHEALTHY" if index % 2 else "HEALTHY"
            position = index * 1500 + (200 if label == "UNHEALTHY" else 0)
            rows.append({"id": f"v{index}", "gene": gene, "label": label, "position": position, "position_fraction": position / 100000,
                         "kind": "SNV", "length_change": 0, "transition": index % 3 == 0, "gc_fraction": 0.5, "cpg": 0,
                         "group": vf.window_group(gene, position), "split": "train"})
        return rows

    def test_grouped_folds_keep_windows_together(self) -> None:
        rows = self._table()
        y = np.array([row["label"] == "UNHEALTHY" for row in rows], dtype=int)
        groups = np.array([row["group"] for row in rows])
        for train, test in vf.grouped_folds(y, groups, n_splits=5, seed=0):
            self.assertFalse(set(groups[train]) & set(groups[test]))

    def test_nearest_neighbour_uses_closest_same_gene_training_variant(self) -> None:
        train = [{"gene": "APP", "position": 100, "label": "HEALTHY"}, {"gene": "APP", "position": 900, "label": "UNHEALTHY"},
                 {"gene": "PSEN2", "position": 120, "label": "UNHEALTHY"}]
        test = [{"gene": "APP", "position": 110}, {"gene": "APP", "position": 800}, {"gene": "PSEN2", "position": 5000}]
        self.assertEqual(vf.nearest_neighbour_predict(train, test), [0, 1, 1])

    def test_feature_sets_separate_position_from_context(self) -> None:
        self.assertIn("position_fraction", vf.FEATURE_SETS["position_only"])
        self.assertNotIn("position_fraction", vf.FEATURE_SETS["context_only"])
        self.assertLessEqual(set(vf.FEATURE_SETS["context_only"]) | {"position_fraction"}, set(vf.FEATURE_SETS["context_and_position"]))

    def test_abstention_margin_meets_target_or_abstains(self) -> None:
        prob = np.array([0.05, 0.1, 0.45, 0.55, 0.9, 0.95])
        y = np.array([0, 0, 1, 0, 1, 1])
        self.assertEqual(vf.choose_margin(prob, y, target_accuracy=1.0), 0.1)  # smallest margin that excludes 0.45 and 0.55
        self.assertEqual(vf.choose_margin(np.array([0.9, 0.1]), np.array([0, 1]), target_accuracy=0.8), 0.5)
        self.assertEqual(list(vf.decide(np.array([0.1, 0.5, 0.9]), 0.3)), ["healthy_label", "uncertain", "unhealthy_label"])

    def test_study_runs_and_reports_every_model(self) -> None:
        results = vf.run_study(self._table(), seeds=[0, 1], n_splits=4)
        self.assertEqual(set(results["models"]), {"majority", "nearest_neighbour", "position_knn", "position_only", "context_only", "context_and_position"})
        for name, summary in results["models"].items():
            with self.subTest(model=name):
                self.assertIn("balanced_accuracy_mean", summary)
                if name in vf.FEATURE_SETS:  # accuracy-when-confident is meaningless without how many were confident
                    self.assertIn("confident_rows_mean", summary)


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class ConfoundAnalysisTests(unittest.TestCase):
    TABLE = [
        {"gene": "APP", "label": "HEALTHY", "position": 100, "kind": "SNV"},
        {"gene": "APP", "label": "HEALTHY", "position": 200, "kind": "SNV"},
        {"gene": "APP", "label": "UNHEALTHY", "position": 900, "kind": "deletion"},
        {"gene": "APP", "label": "UNHEALTHY", "position": 950, "kind": "SNV"},
        {"gene": "PSEN2", "label": "HEALTHY", "position": 500, "kind": "SNV"},
        {"gene": "PSEN2", "label": "UNHEALTHY", "position": 100, "kind": "SNV"},
    ]

    def test_position_knn_votes_within_the_gene(self) -> None:
        train = [{"gene": "APP", "position": p, "label": l} for p, l in
                 ((100, "HEALTHY"), (110, "UNHEALTHY"), (120, "UNHEALTHY"), (5000, "HEALTHY"))]
        self.assertEqual(vf.nearest_neighbour_predict(train, [{"gene": "APP", "position": 105}], k=3), [1])
        self.assertEqual(vf.nearest_neighbour_predict(train, [{"gene": "APP", "position": 105}], k=1), [0])

    def test_class_balance_counts_each_gene(self) -> None:
        balance = vf.class_balance(self.TABLE)
        self.assertEqual(balance["APP"], {"HEALTHY": 2, "UNHEALTHY": 2})
        self.assertEqual(balance["PSEN2"], {"HEALTHY": 1, "UNHEALTHY": 1})

    def test_position_separation_is_one_when_labels_split_by_position(self) -> None:
        separation = vf.position_separation(self.TABLE)
        self.assertEqual(separation["APP"]["position_auc"], 1.0)
        self.assertEqual((separation["APP"]["unhealthy_min"], separation["APP"]["unhealthy_max"]), (900, 950))
        self.assertEqual(separation["PSEN2"]["position_auc"], 1.0)  # direction does not matter

    def test_variant_type_crosstab(self) -> None:
        self.assertEqual(vf.variant_type_by_label(self.TABLE)["deletion"], {"HEALTHY": 0, "UNHEALTHY": 1})

    def test_window_bootstrap_resamples_whole_windows(self) -> None:
        y = np.array([0, 0, 1, 1, 0, 1])
        groups = np.array(["a", "a", "b", "b", "c", "c"])
        low, high = vf.window_bootstrap_ci(y, y.copy(), groups, resamples=200, seed=0)
        self.assertEqual((low, high), (1.0, 1.0))  # perfect predictions stay perfect
        low, high = vf.window_bootstrap_ci(y, 1 - y, groups, resamples=200, seed=0)
        self.assertEqual((low, high), (0.0, 0.0))

    def _study_table(self):
        return StudyMechanicsTests._table(self)

    def test_study_reports_per_gene_results_and_window_is_configurable(self) -> None:
        results = vf.run_study(self._study_table(), seeds=[0], n_splits=3, window=5000)
        self.assertEqual(results["window_bp"], 5000)
        per_gene = results["per_gene"]["nearest_neighbour"]
        self.assertEqual(set(per_gene), {"APP", "PSEN2"})
        self.assertIn("balanced_accuracy_mean", per_gene["APP"])
        self.assertIn("mcc_mean", per_gene["APP"])
        self.assertEqual(len(results["ci95"]["nearest_neighbour"]), 2)

    def test_report_states_limits_and_avoids_clinical_claims(self) -> None:
        import tempfile

        table = self._study_table()
        results = vf.run_study(table, seeds=[0], n_splits=3)
        sensitivity = {1000: results, 5000: vf.run_study(table, seeds=[0], n_splits=3, window=5000)}
        leak = vf.random_split_nearest_neighbour([{**row, "split": "test" if i % 3 == 0 else "train"} for i, row in enumerate(table)])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "REPORT.md"
            vf.write_report(path, table, results, leak, sensitivity)
            text = path.read_text(encoding="utf-8")
        for expected in ("unknown", "Per-gene", "Class balance", "Window-size sensitivity", "Position separation", "Variant type", "Limitations"):
            with self.subTest(expected=expected):
                self.assertIn(expected, text)
        for banned in ("pathogenic", "diagnos", "disease_associated"):
            with self.subTest(banned=banned):
                self.assertNotIn(banned, text.replace("not diagnosis", ""))


@unittest.skipUnless(HAVE_ML_LIBRARIES, "scikit-learn/numpy not installed in this environment")
class RealTableTests(unittest.TestCase):
    def test_variant_table_matches_the_labelled_data(self) -> None:
        table = vf.build_variant_table()
        self.assertEqual(len(table), 220)
        self.assertEqual(sum(row["label"] == "UNHEALTHY" for row in table), 110)
        self.assertEqual(sum(row["split"] == "test" for row in table), 46)
        self.assertTrue(all(row["kind"] in ("SNV", "multi_SNV", "insertion", "deletion", "complex") for row in table))


if __name__ == "__main__":
    unittest.main()
