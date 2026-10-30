"""Feature calculations and label blinding for the replacement Model B."""

from __future__ import annotations

import ast
import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from ml import variant_features_v2 as vf


class Sequence(unittest.TestCase):
    def test_variant_type(self):
        self.assertEqual(vf.variant_type("A", "G"), "snv")
        self.assertEqual(vf.variant_type("AC", "A"), "deletion")
        self.assertEqual(vf.variant_type("A", "AC"), "insertion")
        self.assertEqual(vf.variant_type("AC", "GT"), "mnv")

    def test_substitution_is_strand_neutral(self):
        forward = vf.substitution_and_context("C", "T", "A", "G")
        reverse = vf.substitution_and_context("G", "A", "C", "T")
        self.assertEqual(forward, ("C>T", "ACG", "1"))
        self.assertEqual(reverse, forward)

    def test_transversion_and_indel(self):
        self.assertEqual(vf.substitution_and_context("C", "A", "A", "A")[2], "0")
        self.assertEqual(vf.substitution_and_context("AC", "A", "A", "A"), ("none", "none", ""))

    def test_cpg_on_either_side(self):
        self.assertEqual(vf.in_cpg("AACGA", 2, 3), 1)
        self.assertEqual(vf.in_cpg("AACGA", 3, 4), 1)
        self.assertEqual(vf.in_cpg("AACGA", 1, 2), 0)
        self.assertEqual(vf.in_cpg("AAGCA", 2, 3), 0)

    def test_gc_leaves_out_variant_bases(self):
        self.assertEqual(vf.gc_fraction("AAGAA", 2, 3, 2), 0.0)
        self.assertEqual(vf.gc_fraction("GAAAG", 2, 3, 2), 0.5)

    def test_homopolymer_run(self):
        self.assertEqual(vf.homopolymer_run("CAAAAG", 2, 3), 3)
        self.assertEqual(vf.homopolymer_run("CAGT", 1, 2), 1)

    def test_main_features_refuses_wrong_reference(self):
        row = {"variation_id": "1", "ref": "T", "alt": "G", "pos": "13"}
        with self.assertRaises(ValueError):
            vf.main_features(row, "AAACAAA", 10)

    def test_main_features_uses_offset(self):
        row = {"variation_id": "1", "ref": "C", "alt": "T", "pos": "14"}
        result = vf.main_features(row, "AAACGAA", 10)
        self.assertEqual((result["substitution"], result["trinucleotide"], result["in_cpg"]), ("C>T", "ACG", 1))


class Protein(unittest.TestCase):
    def test_properties(self):
        from Bio.Align import substitution_matrices
        result = vf.missense_features("1", "G/P", substitution_matrices.load("BLOSUM62"))
        self.assertEqual(result["blosum62"], -2)
        self.assertEqual(result["delta_hydropathy"], -1.2)
        self.assertEqual(result["delta_volume"], 52.6)
        self.assertEqual((result["to_proline"], result["from_glycine"], result["delta_charge"]), (1, 1, 0))

    def test_charge(self):
        from Bio.Align import substitution_matrices
        self.assertEqual(vf.missense_features("1", "D/K", substitution_matrices.load("BLOSUM62"))["delta_charge"], 2)

    def test_change_spanning_two_codons(self):
        self.assertEqual(vf.single_amino_change("VI/VV"), ("I", "V"))
        self.assertEqual(vf.single_amino_change("A/C"), ("A", "C"))

    def test_changes_that_are_not_one_residue(self):
        for value in ("VI/LV", "V/VI", "A/*", "", "A/A"):
            self.assertIsNone(vf.single_amino_change(value), value)

    def test_reference_allele_record_is_not_a_variant(self):
        self.assertTrue(vf.is_no_change("T", "T"))
        self.assertFalse(vf.is_no_change("T", "C"))

    def test_tables_cover_all_amino_acids(self):
        self.assertEqual(set(vf.HYDROPATHY), set(vf.VOLUME))
        self.assertEqual(len(vf.VOLUME), 20)

    def test_canonical_transcript_must_match_gene(self):
        item = {"transcript_consequences": [{"canonical": 1, "gene_symbol": "OTHER"}, {"canonical": 1, "gene_symbol": "APP", "amino_acids": "A/V"}, {"gene_symbol": "APP"}]}
        self.assertEqual(vf.canonical_consequence(item, "APP")["amino_acids"], "A/V")
        self.assertIsNone(vf.canonical_consequence(item, "PSEN1"))
        self.assertIsNone(vf.canonical_consequence(None, "APP"))


class Blinding(unittest.TestCase):
    FORBIDDEN = {"label", "label_code", "aggregate_class", "review_status", "stars", "number_of_submitters", "sole_submitter", "last_evaluated", "first_seen_release", "name", "rcvs", "vcv"}

    def test_only_identity_columns_are_read(self):
        self.assertFalse(self.FORBIDDEN & set(vf.IDENTITY_COLUMNS))
        self.assertEqual(set(vf.BASELINE_INPUT_COLUMNS), {"variation_id", "consequence_group"})

    def test_read_columns_returns_nothing_else(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.csv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["variation_id", "label"])
                writer.writerow(["1", "pathogenic"])
            self.assertEqual(vf.read_columns(path, ("variation_id",)), [{"variation_id": "1"}])

    def test_source_does_not_use_forbidden_columns_as_strings(self):
        tree = ast.parse(Path(vf.__file__).read_text(encoding="utf-8"))
        strings = {node.value for node in ast.walk(tree) if isinstance(node, ast.Constant) and isinstance(node.value, str)}
        self.assertFalse(self.FORBIDDEN & strings)

    def test_feature_tables_exclude_gene_position_and_consequence(self):
        banned = {"gene", "pos", "chrom", "relative_position", "consequence_group", "vep_consequence", "protein_position"}
        self.assertFalse(banned & set(vf.MAIN_COLUMNS))
        self.assertFalse(banned & set(vf.MISSENSE_COLUMNS))


if __name__ == "__main__":
    unittest.main()
