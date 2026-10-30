"""Labelling and features for the Model A (alignment reliability) dataset."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from ml.build_alignment_dataset import (
    FEATURE_COLUMNS,
    alignment_features,
    label_alignment,
    read_group,
    rows_for_gaf,
    setting_name,
)

REF = "ref|NC_000021.9|:25880550-26171128"
READ = f"APP_H_3.S1_1!{REF}!1000!2000!+"


def gaf_fields(path=REF, strand="+", start=1000, end=2000, matches=950, block=1000, mapq=60, path_length=290579,
               query_length=1000, query_start=0, query_end=1000, tags=("tp:A:P", "NM:i:50", "cm:i:80", "s1:i:900", "s2:i:300", "dv:f:0.05")):
    return [READ, str(query_length), str(query_start), str(query_end), strand, path, str(path_length), str(start), str(end),
            str(matches), str(block), str(mapq), *tags]


class LabelTests(unittest.TestCase):
    def test_on_reference_overlap_is_correct(self) -> None:
        self.assertEqual(label_alignment(gaf_fields(), READ, "APP", "APP", REF), "correct")

    def test_wrong_place_or_strand_is_incorrect(self) -> None:
        self.assertEqual(label_alignment(gaf_fields(start=50000, end=51000), READ, "APP", "APP", REF), "incorrect")
        self.assertEqual(label_alignment(gaf_fields(strand="-"), READ, "APP", "APP", REF), "incorrect")

    def test_any_alignment_to_another_gene_graph_is_incorrect(self) -> None:
        decoy = gaf_fields(path="ref|NC_000014.9|:73136434-73223709", start=1000, end=2000)
        self.assertEqual(label_alignment(decoy, READ, "APP", "PSEN1", "ref|NC_000014.9|:73136434-73223709"), "incorrect")

    def test_off_reference_path_is_projected_not_counted_wrong(self) -> None:
        path = f">{REF}:0-1500>APP_H_9:1500-1501>{REF}:1501-290579"
        self.assertEqual(label_alignment(gaf_fields(path=path, start=1000, end=2000), READ, "APP", "APP", REF), "correct")

    def test_unprojectable_path_gets_its_own_status(self) -> None:
        path = f">{REF}:0-1500>APP_H_9:0-10"
        self.assertEqual(label_alignment(gaf_fields(path=path, start=1000, end=1505), READ, "APP", "APP", REF), "unprojectable")


class FeatureTests(unittest.TestCase):
    def test_features_come_from_columns_and_tags(self) -> None:
        features = alignment_features(gaf_fields(query_start=10, query_end=1000), alignment_count=2)
        self.assertEqual(set(features), set(FEATURE_COLUMNS))
        self.assertEqual(features["mapq"], 60)
        self.assertAlmostEqual(features["identity"], 0.95)
        self.assertAlmostEqual(features["query_coverage"], 0.99)
        self.assertAlmostEqual(features["nm_rate"], 0.05)
        self.assertEqual((features["chain_minimizers"], features["chain_score"], features["second_chain_score"]), (80, 900, 300))
        self.assertAlmostEqual(features["second_to_best_chain_ratio"], 300 / 900)
        self.assertAlmostEqual(features["divergence"], 0.05)
        self.assertEqual((features["path_segments"], features["off_reference_path"], features["alignment_count"]), (1, 0, 2))

    def test_missing_tags_become_explicit_defaults(self) -> None:
        features = alignment_features(gaf_fields(tags=()), alignment_count=1)
        self.assertEqual((features["chain_score"], features["second_chain_score"], features["divergence"]), (0, 0, -1.0))
        self.assertEqual(features["second_to_best_chain_ratio"], 0.0)

    def test_group_is_the_source_haplotype(self) -> None:
        self.assertEqual(read_group(READ), "APP_H_3")
        self.assertEqual(setting_name(0.85, 2000), "acc85_len2000")


class RowsTests(unittest.TestCase):
    def test_best_alignment_per_read_and_unaligned_rows(self) -> None:
        other = f"PSEN1_U_2.S1_9!{REF}!5000!6000!-"
        weaker = gaf_fields(start=50000, end=51000, matches=400)
        with tempfile.TemporaryDirectory() as tmp:
            gaf = Path(tmp) / "x.gaf"
            gaf.write_text("\t".join(gaf_fields()) + "\n" + "\t".join(weaker) + "\n", encoding="utf-8")
            rows = rows_for_gaf([READ, other], gaf, setting="acc95_len10000", seed=1, aligner="panaligner",
                                own_gene="APP", graph_gene="APP", reference_name=REF)
        by_read = {row["read"]: row for row in rows}
        self.assertEqual(by_read[READ]["label"], "correct")  # longest alignment wins, as in mapeval
        self.assertEqual(by_read[READ]["alignment_count"], 2)
        self.assertEqual((by_read[READ]["group"], by_read[READ]["decoy"]), ("APP_H_3", 0))
        self.assertEqual(by_read[other]["label"], "unaligned")
        self.assertEqual(by_read[other]["group"], "PSEN1_U_2")


if __name__ == "__main__":
    unittest.main()
