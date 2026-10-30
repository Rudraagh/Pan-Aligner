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
        self.assertEqual(list(vf.decide(np.array([0.1, 0.5, 0.9]), 0.3)), ["healthy_associated", "uncertain", "disease_associated"])

    def test_study_runs_and_reports_every_model(self) -> None:
        results = vf.run_study(self._table(), seeds=[0, 1], n_splits=4)
        self.assertEqual(set(results["models"]), {"majority", "nearest_neighbour", "position_only", "context_only", "context_and_position"})
        for name, summary in results["models"].items():
            with self.subTest(model=name):
                self.assertIn("balanced_accuracy_mean", summary)


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
