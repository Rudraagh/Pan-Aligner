"""Corrected Model B snapshot: invalid records, genomic clusters, split manifests, leakage checks and feature audit."""

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


def variants(*items):
    """(id, gene, pos, ref, alt, residue, label)."""
    return pd.DataFrame(
        [{"variation_id": str(i), "gene": g, "pos": str(p), "ref": r, "alt": a, "residue": res, "label": lab} for i, g, p, r, a, res, lab in items]
    )


class InvalidRecords(unittest.TestCase):
    def row(self, **changes):
        base = {"variation_id": "1", "vcv": "VCV000000001", "ref": "A", "alt": "G", "label": "benign", "gene": "APP", "chrom": "21", "pos": "10", "canonical_spdi": "x"}
        return {**base, **changes}

    def test_valid_record_has_no_reasons(self):
        self.assertEqual(vd.invalid_reasons(self.row()), [])

    def test_identical_alleles(self):
        self.assertIn("REF equals ALT: reference-allele record, not a variant", vd.invalid_reasons(self.row(alt="A")))

    def test_non_nucleotide_allele(self):
        self.assertTrue(vd.invalid_reasons(self.row(alt="N")))
        self.assertTrue(vd.invalid_reasons(self.row(ref="")))

    def test_label_outside_the_two_classes(self):
        self.assertTrue(vd.invalid_reasons(self.row(label="uncertain")))
        self.assertTrue(vd.invalid_reasons(self.row(label="")))

    def test_accession_must_match_variation_id(self):
        self.assertTrue(vd.invalid_reasons(self.row(vcv="VCV000000002")))

    def test_unknown_gene(self):
        self.assertTrue(vd.invalid_reasons(self.row(gene="SORL1")))


class Clusters(unittest.TestCase):
    def test_variants_within_gap_share_a_cluster(self):
        frame = variants((1, "APP", 100, "A", "G", "", "benign"), (2, "APP", 1100, "A", "G", "", "benign"))
        clusters = vd.cluster_variants(frame, max_gap=1000)
        self.assertEqual(clusters.nunique(), 1)

    def test_variants_beyond_gap_are_separate(self):
        frame = variants((1, "APP", 100, "A", "G", "", "benign"), (2, "APP", 1101, "A", "G", "", "benign"))
        self.assertEqual(vd.cluster_variants(frame, max_gap=1000).nunique(), 2)

    def test_gap_is_measured_from_the_end_of_a_deletion(self):
        frame = variants((1, "APP", 100, "A" * 10, "A", "", "benign"), (2, "APP", 1109, "A", "G", "", "benign"))
        self.assertEqual(vd.cluster_variants(frame, max_gap=1000).nunique(), 1)

    def test_chain_links_transitively(self):
        frame = variants((1, "APP", 100, "A", "G", "", "benign"), (2, "APP", 1000, "A", "G", "", "benign"), (3, "APP", 1900, "A", "G", "", "benign"))
        self.assertEqual(vd.cluster_variants(frame, max_gap=1000).nunique(), 1)

    def test_genes_never_merge(self):
        frame = variants((1, "APP", 100, "A", "G", "", "benign"), (2, "PSEN1", 100, "A", "G", "", "benign"))
        self.assertEqual(vd.cluster_variants(frame, max_gap=1000).nunique(), 2)

    def test_same_residue_merges_across_a_long_gap(self):
        frame = variants((1, "APP", 100, "A", "G", "Ala5", "benign"), (2, "APP", 9000, "A", "G", "Ala5", "pathogenic"))
        self.assertEqual(vd.cluster_variants(frame, max_gap=1000).nunique(), 1)

    def test_order_of_rows_does_not_change_groups(self):
        frame = variants((1, "APP", 100, "A", "G", "", "benign"), (2, "APP", 500, "A", "G", "", "benign"), (3, "APP", 9000, "A", "G", "", "benign"))
        forward = vd.cluster_variants(frame, max_gap=1000)
        backward = vd.cluster_variants(frame.iloc[::-1].reset_index(drop=True), max_gap=1000)
        self.assertEqual(dict(zip(frame["variation_id"], forward)), dict(zip(frame.iloc[::-1]["variation_id"], backward)))


def dataset(n_groups=10, per_group=4):
    rows = []
    identifier = 1
    for group in range(n_groups):
        for member in range(per_group):
            label = "pathogenic" if member % 2 == 0 else "benign"
            rows.append((identifier, "APP", 10_000 * (group + 1) + member, "A", "G", "", label))
            identifier += 1
    frame = variants(*rows)
    frame["cluster"] = vd.cluster_variants(frame, max_gap=1000)
    return frame


class FoldManifest(unittest.TestCase):
    def test_every_record_has_one_role_per_fold(self):
        frame = dataset()
        manifest = vd.fold_manifest(frame, "s", seeds=[0], n_folds=5)
        self.assertEqual(len(manifest), 5 * len(frame))
        self.assertTrue((manifest.groupby(["seed", "fold", "variation_id"]).size() == 1).all())

    def test_each_record_is_tested_exactly_once_per_seed(self):
        frame = dataset()
        manifest = vd.fold_manifest(frame, "s", seeds=[0, 1], n_folds=5)
        tested = manifest[manifest["role"] == "test"].groupby(["seed", "variation_id"]).size()
        self.assertTrue((tested == 1).all())
        self.assertEqual(len(tested), 2 * len(frame))

    def test_roles_are_train_validation_test(self):
        manifest = vd.fold_manifest(dataset(), "s", seeds=[0], n_folds=5)
        self.assertEqual(set(manifest["role"]), {"train", "validation", "test"})

    def test_seed_fixes_the_result(self):
        frame = dataset()
        pd.testing.assert_frame_equal(vd.fold_manifest(frame, "s", [3], 5), vd.fold_manifest(frame, "s", [3], 5))

    def test_clusters_are_not_split(self):
        frame = dataset()
        manifest = vd.fold_manifest(frame, "s", [0], 5).merge(frame[["variation_id", "cluster"]], on="variation_id")
        self.assertTrue((manifest.groupby(["seed", "fold", "cluster"])["role"].nunique() == 1).all())


class GeneHeldOut(unittest.TestCase):
    def test_held_out_gene_is_only_in_test(self):
        frame = pd.concat([dataset(), dataset().assign(gene="PSEN1", variation_id=lambda d: "9" + d["variation_id"])], ignore_index=True)
        frame["cluster"] = frame["gene"] + frame["cluster"]
        manifest = vd.gene_held_out_manifest(frame, "s", seeds=[0], n_folds=5).merge(frame[["variation_id", "gene"]], on="variation_id")
        for held, part in manifest.groupby("held_out_gene"):
            self.assertEqual(set(part[part["role"] == "test"]["gene"]), {held})
            self.assertNotIn(held, set(part[part["role"] != "test"]["gene"]))


class GeneHeldOutValidation(unittest.TestCase):
    def frame(self):
        rows, identifier = [], 1
        for cluster in range(12):
            rows.append((identifier, "PSEN1", 10_000 * (cluster + 1), "A", "G", "", "pathogenic"))
            identifier += 1
        for cluster in range(6):
            rows.append((identifier, "PSEN2", 10_000 * (cluster + 1), "A", "G", "", "benign"))
            identifier += 1
        rows.append((identifier, "APP", 10_000, "A", "G", "", "benign"))
        rows.append((identifier + 1, "APP", 20_000, "A", "G", "", "pathogenic"))
        frame = variants(*rows)
        frame["cluster"] = vd.cluster_variants(frame, max_gap=1000)
        return frame

    def test_first_fold_with_both_classes_is_chosen(self):
        fold_of = np.array([0, 0, 1, 1, 2, 2])
        labels = np.array(["pathogenic", "pathogenic", "pathogenic", "benign", "pathogenic", "benign"])
        self.assertEqual(vd.choose_validation_fold(fold_of, labels), 1)

    def test_falls_back_to_fold_zero_when_no_fold_has_both(self):
        fold_of = np.array([0, 0, 1, 1])
        labels = np.array(["pathogenic", "pathogenic", "benign", "benign"])
        self.assertEqual(vd.choose_validation_fold(fold_of, labels), 0)

    def test_validation_holds_both_classes_when_possible(self):
        frame = self.frame()
        manifest = vd.gene_held_out_manifest(frame, "s", seeds=range(10), n_folds=5).merge(frame[["variation_id", "label"]], on="variation_id")
        chosen = manifest[(manifest["held_out_gene"] == "APP") & (manifest["role"] == "validation")]
        self.assertTrue((chosen.groupby("seed")["label"].nunique() == 2).all())


class AuditManifest(unittest.TestCase):
    def good(self):
        frame = dataset()
        return frame, vd.fold_manifest(frame, "s", [0], 5)

    def test_clean_manifest_has_no_failures(self):
        frame, manifest = self.good()
        self.assertEqual(vd.audit_manifest(frame, manifest, max_gap=1000), [])

    def test_detects_same_variant_with_two_roles(self):
        frame, manifest = self.good()
        extra = manifest.iloc[[0]].assign(role="test" if manifest.iloc[0]["role"] != "test" else "train")
        self.assertTrue(vd.audit_manifest(frame, pd.concat([manifest, extra]), max_gap=1000))

    def test_detects_nearby_variants_crossing(self):
        frame, manifest = self.good()
        moved = frame.copy()
        target = manifest[(manifest["fold"] == 0) & (manifest["role"] == "test")]["variation_id"].iloc[0]
        other = manifest[(manifest["fold"] == 0) & (manifest["role"] == "train")]["variation_id"].iloc[0]
        moved.loc[moved["variation_id"] == target, "pos"] = str(int(moved.loc[moved["variation_id"] == other, "pos"].iloc[0]) + 5)
        failures = vd.audit_manifest(moved, manifest, max_gap=1000)
        self.assertTrue(any("within" in failure for failure in failures))

    def test_detects_unknown_variant(self):
        frame, manifest = self.good()
        extra = manifest.iloc[[0]].assign(variation_id="999999")
        self.assertTrue(vd.audit_manifest(frame, pd.concat([manifest, extra]), max_gap=1000))

    def test_detects_excluded_variant(self):
        frame, manifest = self.good()
        self.assertTrue(vd.audit_manifest(frame.iloc[1:], manifest, max_gap=1000))

    def test_detects_same_residue_crossing(self):
        frame = dataset()
        frame["residue"] = ""
        manifest = vd.fold_manifest(frame, "s", [0], 5)
        test_id = manifest[(manifest["fold"] == 0) & (manifest["role"] == "test")]["variation_id"].iloc[0]
        train_id = manifest[(manifest["fold"] == 0) & (manifest["role"] == "train")]["variation_id"].iloc[0]
        frame.loc[frame["variation_id"].isin([test_id, train_id]), "residue"] = "Ala5"
        self.assertTrue(any("residue" in failure for failure in vd.audit_manifest(frame, manifest, max_gap=1000)))


class FeatureColumns(unittest.TestCase):
    def test_clean_columns_pass(self):
        self.assertEqual(vd.forbidden_feature_columns(["gc_50", "in_cpg", "blosum62"]), [])

    def test_label_and_evidence_columns_are_caught(self):
        for column in ("label", "label_code", "aggregate_class", "review_status", "stars", "sole_submitter", "condition", "vcv", "rcvs", "last_evaluated", "first_seen_release", "name"):
            self.assertEqual(vd.forbidden_feature_columns([column]), [column])

    def test_gene_position_and_consequence_are_caught(self):
        for column in ("gene", "pos", "chrom", "relative_position", "consequence_group", "vep_consequence", "protein_position", "window_1kb", "cluster"):
            self.assertEqual(vd.forbidden_feature_columns([column]), [column])


class FeatureAudit(unittest.TestCase):
    def frame(self):
        return pd.DataFrame({
            "variation_id": [str(i) for i in range(8)],
            "good": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
            "flat": [1] * 8,
            "gaps": [None, 1.0, 2.0, None, 3.0, 1.0, 2.0, 3.0],
            "kind": ["a", "b", "a", "b", "a", "b", "a", "b"],
        })

    def context(self):
        return pd.DataFrame({
            "variation_id": [str(i) for i in range(8)],
            "label": ["benign"] * 4 + ["pathogenic"] * 4,
            "consequence_group": ["missense"] * 4 + ["intronic"] * 4,
            "gene": ["APP", "PSEN1"] * 4,
        })

    def audit(self):
        return vd.feature_audit(self.frame(), self.context(), table="t").set_index("feature")

    def test_constant_feature_flagged(self):
        self.assertIn("constant", self.audit().loc["flat", "flags"])

    def test_missing_values_counted(self):
        self.assertEqual(self.audit().loc["gaps", "n_missing"], 2)

    def test_perfect_numeric_association_gives_auc_one(self):
        row = self.audit().loc["good"]
        self.assertEqual(row["label_auc"], 1.0)
        self.assertEqual(row["missense_auc"], 0.0)

    def test_categorical_uses_cramers_v(self):
        row = self.audit().loc["kind"]
        self.assertEqual(row["kind"], "categorical")
        self.assertEqual(row["gene_association"], 1.0)
        self.assertEqual(row["label_association"], 0.0)

    def test_proxy_flag(self):
        self.assertIn("possible_label_proxy", self.audit().loc["good", "flags"])

    def test_id_is_not_a_feature_row(self):
        self.assertNotIn("variation_id", self.audit().index)


class Joins(unittest.TestCase):
    def test_one_to_one_passes(self):
        self.assertEqual(vd.join_problems(["1", "2"], ["2", "1"], "t"), [])

    def test_missing_extra_and_duplicate_rows(self):
        self.assertTrue(vd.join_problems(["1", "2"], ["1"], "t"))
        self.assertTrue(vd.join_problems(["1"], ["1", "2"], "t"))
        self.assertTrue(vd.join_problems(["1"], ["1", "1"], "t"))

    def test_subset_mode_allows_missense_table(self):
        self.assertEqual(vd.join_problems(["1", "2"], ["1"], "t", subset=True), [])
        self.assertTrue(vd.join_problems(["1"], ["1", "2"], "t", subset=True))


if __name__ == "__main__":
    unittest.main()
