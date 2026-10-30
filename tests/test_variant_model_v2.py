"""Replacement Model B: thresholds, baselines, separation of roles, and interval logic. Synthetic data only."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from ml import variant_dataset_validation as vd
from ml import variant_model_v2 as vm


class Threshold(unittest.TestCase):
    def test_separates_when_possible(self):
        threshold = vm.choose_threshold(np.array([0.1, 0.2, 0.8, 0.9]), np.array([0, 0, 1, 1]))
        self.assertEqual(((np.array([0.1, 0.2, 0.8, 0.9]) >= threshold).astype(int)).tolist(), [0, 0, 1, 1])

    def test_tie_goes_nearest_one_half(self):
        self.assertEqual(vm.choose_threshold(np.array([0.1, 0.9]), np.array([0, 1])), 0.5)

    def test_tie_between_real_candidates_prefers_nearest_one_half(self):
        threshold = vm.choose_threshold(np.array([0.1, 0.2, 0.3, 0.4]), np.array([0, 0, 1, 1]))
        self.assertAlmostEqual(threshold, 0.25)

    def test_identical_scores_are_still_handled(self):
        threshold = vm.choose_threshold(np.array([0.5, 0.5]), np.array([0, 1]))
        self.assertTrue(np.isfinite(threshold))


class Metrics(unittest.TestCase):
    def test_balanced_accuracy(self):
        self.assertEqual(vm.balanced_accuracy(np.array([1, 1, 0, 0]), np.array([1, 0, 0, 0])), 0.75)

    def test_one_class_is_nan(self):
        self.assertTrue(np.isnan(vm.balanced_accuracy(np.array([1, 1]), np.array([1, 0]))))


def table(*rows):
    """(gene, pos, y, consequence)."""
    return pd.DataFrame([{"gene": g, "pos": p, "y": y, "consequence_group": c} for g, p, y, c in rows])


class Baselines(unittest.TestCase):
    def test_gene_score_is_smoothed_train_rate(self):
        train = table(("A", 1, 1, "m"), ("A", 2, 1, "m"), ("B", 3, 0, "m"), ("B", 4, 0, "m"))
        target = table(("A", 9, 0, "m"), ("B", 9, 1, "m"), ("C", 9, 1, "m"))
        scores = vm.group_rate_scores(train, target, "gene")
        self.assertGreater(scores[0], scores[1])
        self.assertAlmostEqual(scores[2], 0.5)

    def test_target_labels_do_not_change_group_scores(self):
        train = table(("A", 1, 1, "m"), ("B", 3, 0, "m"))
        first = vm.group_rate_scores(train, table(("A", 9, 0, "m")), "gene")
        second = vm.group_rate_scores(train, table(("A", 9, 1, "m")), "gene")
        self.assertEqual(first.tolist(), second.tolist())

    def test_nearest_uses_same_gene_only(self):
        train = table(("A", 1000, 1, "m"), ("B", 101, 0, "m"))
        scores = vm.nearest_variant_scores(train, table(("A", 100, 0, "m")))
        self.assertEqual(scores.tolist(), [1.0])

    def test_nearest_picks_closest(self):
        train = table(("A", 100, 1, "m"), ("A", 1000, 0, "m"))
        self.assertEqual(vm.nearest_variant_scores(train, table(("A", 900, 0, "m"))).tolist(), [0.0])

    def test_nearest_falls_back_to_train_rate_for_unseen_gene(self):
        train = table(("A", 100, 1, "m"), ("A", 1000, 0, "m"))
        self.assertEqual(vm.nearest_variant_scores(train, table(("Z", 900, 0, "m"))).tolist(), [0.5])

    def test_position_vote_averages_nearest_five_in_gene(self):
        train = table(*[("A", p, y, "m") for p, y in ((1, 1), (2, 1), (3, 1), (4, 0), (5, 0), (100, 0), (101, 0))], ("B", 3, 0, "m"))
        self.assertEqual(vm.position_vote_scores(train, table(("A", 3, 0, "m"))).tolist(), [0.6])


def frame(n_clusters=20, per=4, signal=True, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for cluster in range(n_clusters):
        for member in range(per):
            y = member % 2
            rows.append({
                "variation_id": f"{cluster}_{member}", "gene": "APP" if cluster % 2 else "PSEN1", "pos": 10_000 * (cluster + 1) + member, "cluster": f"c{cluster}",
                "consequence_group": "missense", "y": y, "gc_50": rng.random(), "gc_200": rng.random(),
                "good": float(y) + rng.normal(0, 0.1) if signal else rng.normal(), "kind": "a" if y else "b",
            })
    return pd.DataFrame(rows)


def roles(data, test_clusters, validation_clusters):
    role = np.full(len(data), "train", dtype=object)
    role[data["cluster"].isin(validation_clusters)] = "validation"
    role[data["cluster"].isin(test_clusters)] = "test"
    return pd.Series(role, index=data.index)


class Arms(unittest.TestCase):
    MAIN = ["variant_type", "gc_50", "gc_200", "in_cpg"]
    AMINO = ["aa_ref", "blosum62"]

    def test_missense_arm_has_both_tables(self):
        self.assertEqual(vm.arm_columns("missense", self.MAIN, self.AMINO), self.MAIN + self.AMINO)

    def test_no_gc_arms_drop_only_gc(self):
        self.assertEqual(vm.arm_columns("missense_no_gc", self.MAIN, self.AMINO), ["variant_type", "in_cpg"] + self.AMINO)
        self.assertEqual(vm.arm_columns("main_no_gc", self.MAIN, self.AMINO), ["variant_type", "in_cpg"])

    def test_amino_only(self):
        self.assertEqual(vm.arm_columns("amino_only", self.MAIN, self.AMINO), self.AMINO)

    def test_every_arm_is_free_of_forbidden_columns(self):
        for arm in vm.ARMS:
            self.assertEqual(vd.forbidden_feature_columns(vm.arm_columns(arm, self.MAIN, self.AMINO)), [])

    def test_forbidden_column_is_refused(self):
        with self.assertRaises(ValueError):
            vm.arm_columns("main", ["gc_50", "gene"], self.AMINO)


class FoldSeparation(unittest.TestCase):
    def run_fold(self, data, role, model="logistic"):
        return vm.fit_fold(data, role, ["good", "kind"], model, seed=0)

    def test_learns_real_signal(self):
        data = frame()
        result = self.run_fold(data, roles(data, ["c0", "c1", "c2", "c3"], ["c4", "c5", "c6", "c7"]))
        test = data[roles(data, ["c0", "c1", "c2", "c3"], ["c4", "c5", "c6", "c7"]) == "test"]
        self.assertEqual(vm.balanced_accuracy(test["y"].to_numpy(), (result["test_scores"] >= result["threshold"]).astype(int)), 1.0)

    def test_gradient_boosting_also_runs(self):
        data = frame()
        result = self.run_fold(data, roles(data, ["c0", "c1"], ["c2", "c3"]), model="gradient_boosting")
        self.assertEqual(len(result["test_scores"]), 8)

    def test_changing_test_labels_changes_nothing(self):
        data = frame()
        role = roles(data, ["c0", "c1", "c2", "c3"], ["c4", "c5", "c6", "c7"])
        before = self.run_fold(data, role)
        changed = data.copy()
        changed.loc[role == "test", "y"] = 1 - changed.loc[role == "test", "y"]
        after = self.run_fold(changed, role)
        self.assertEqual(before["test_scores"].tolist(), after["test_scores"].tolist())
        self.assertEqual((before["threshold"], before["params"]), (after["threshold"], after["params"]))

    def test_changing_test_features_does_not_change_threshold(self):
        data = frame()
        role = roles(data, ["c0", "c1"], ["c2", "c3"])
        before = self.run_fold(data, role)
        changed = data.copy()
        changed.loc[role == "test", "good"] = 100.0
        self.assertEqual(before["threshold"], self.run_fold(changed, role)["threshold"])

    def test_validation_labels_affect_only_threshold_and_settings(self):
        data = frame()
        role = roles(data, ["c0", "c1"], ["c2", "c3", "c4", "c5"])
        changed = data.copy()
        changed.loc[role == "validation", "y"] = 1 - changed.loc[role == "validation", "y"]
        self.assertNotEqual(self.run_fold(data, role)["threshold"], self.run_fold(changed, role)["threshold"])

    def test_category_only_in_test_is_ignored(self):
        data = frame()
        role = roles(data, ["c0", "c1"], ["c2", "c3"])
        changed = data.copy()
        changed.loc[role == "test", "kind"] = "never_seen"
        self.assertEqual(len(self.run_fold(changed, role)["test_scores"]), 8)

    def test_missing_numeric_value_is_imputed(self):
        data = frame()
        data.loc[data.index[:6], "good"] = np.nan
        result = self.run_fold(data, roles(data, ["c0", "c1"], ["c2", "c3"]))
        self.assertTrue(np.isfinite(result["test_scores"]).all())

    def test_row_order_does_not_decide_result(self):
        data = frame()
        role = roles(data, ["c0", "c1"], ["c2", "c3"])
        shuffled = data.sample(frac=1, random_state=3)
        first = self.run_fold(data, role)
        second = self.run_fold(shuffled, role.loc[shuffled.index])
        order = data.index[role == "test"]
        placed = pd.Series(second["test_scores"], index=shuffled.index[role.loc[shuffled.index] == "test"]).loc[order].to_numpy()
        np.testing.assert_allclose(first["test_scores"], placed, atol=1e-6)


class BaselineChoice(unittest.TestCase):
    def test_uses_validation_not_test(self):
        chosen = vm.choose_baseline({"a": 0.6, "b": 0.8})
        self.assertEqual(chosen, "b")

    def test_tie_is_deterministic(self):
        self.assertEqual(vm.choose_baseline({"b": 0.7, "a": 0.7}), "a")


class Bootstrap(unittest.TestCase):
    def data(self):
        y = np.array([1, 0] * 20)
        clusters = np.repeat(np.arange(10), 4)
        return y, clusters

    def test_identical_predictions_have_zero_difference(self):
        y, clusters = self.data()
        result = vm.cluster_bootstrap_difference(y, y, y, clusters, n_resamples=100, seed=1)
        self.assertEqual((result["difference"], result["low"], result["high"]), (0.0, 0.0, 0.0))

    def test_better_model_gives_positive_difference(self):
        y, clusters = self.data()
        result = vm.cluster_bootstrap_difference(y, y, 1 - y, clusters, n_resamples=100, seed=1)
        self.assertEqual(result["difference"], 1.0)
        self.assertGreater(result["low"], 0)

    def test_same_seed_repeats(self):
        y, clusters = self.data()
        noisy = y.copy()
        noisy[::3] = 1 - noisy[::3]
        first = vm.cluster_bootstrap_difference(y, noisy, y, clusters, 50, 7)
        self.assertEqual(first, vm.cluster_bootstrap_difference(y, noisy, y, clusters, 50, 7))

    def test_clumped_errors_widen_the_interval(self):
        y = np.array([1, 0] * 100)
        blocks = np.arange(200) // 25
        a = np.where(blocks % 2 == 0, y, 1 - y)
        b = np.where(blocks % 3 == 0, y, 1 - y)
        rows_as_clusters = vm.cluster_bootstrap_difference(y, a, b, np.arange(200), 300, 1)
        true_blocks = vm.cluster_bootstrap_difference(y, a, b, blocks, 300, 1)
        self.assertGreater(true_blocks["high"] - true_blocks["low"], rows_as_clusters["high"] - rows_as_clusters["low"])


class PermutedControl(unittest.TestCase):
    def test_permutation_moves_whole_clusters(self):
        data = frame()
        permuted = vm.permute_labels_by_cluster(data, seed=1)
        self.assertEqual(permuted.groupby("cluster")["y"].sum().sort_values().tolist(), data.groupby("cluster")["y"].sum().sort_values().tolist())
        self.assertEqual(len(permuted), len(data))

    def test_permutation_keeps_class_share_near(self):
        data = frame()
        self.assertEqual(vm.permute_labels_by_cluster(data, seed=1)["y"].sum(), data["y"].sum())


if __name__ == "__main__":
    unittest.main()
