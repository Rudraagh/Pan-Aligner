"""Pilot gates, k-mer check and input preparation."""

from __future__ import annotations

import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

try:
    import numpy as np
except ImportError:
    np = None

from ml.wg_truth.pilot_report import compute
from ml.wg_truth.prepare_inputs import read_records, sequence_md5
from ml.wg_truth.simulate_haplotypes import setting_seed


def row(state="correct", aligner="panaligner", setting="a", start=0, kind="clean", error="", reason="ok"):
    return {"truth_state": state, "aligner": aligner, "setting": setting, "truth_ref_chrom": "chr21", "truth_ref_start": str(start), "path_projected": kind, "error_type": error, "truth_reason": reason}


def passing_rows():
    rows = []
    for setting in ("a", "b", "c", "d", "e", "f"):
        for aligner in ("panaligner", "minigraph"):
            rows.append(row(setting=setting, aligner=aligner))
    for index in range(600):
        rows.append(row("incorrect", "panaligner" if index % 2 else "minigraph", start=index * 100_000, error="wrong_chromosome"))
    rows += [row() for _ in range(5000)]
    return rows


GOOD_KMER = Counter(agree=95, disagree=1, not_assessable=4)


class Gates(unittest.TestCase):
    def test_all_pass_is_go(self):
        report = compute(passing_rows(), GOOD_KMER, True, True)
        self.assertEqual(report["verdict"], "GO")
        self.assertEqual({g["result"] for g in report["gates"].values()}, {"PASS"})

    def test_missing_evidence_is_never_a_pass(self):
        report = compute(passing_rows(), None, None, None)
        self.assertEqual(report["verdict"], "NO_GO")
        for name in ("independence", "runs", "reproducibility"):
            self.assertEqual(report["gates"][name]["result"], "NOT_EVALUATED")

    def test_too_few_errors(self):
        rows = [r for r in passing_rows() if r["truth_state"] != "incorrect"] + [row("incorrect", start=i * 100_000) for i in range(499)]
        report = compute(rows, GOOD_KMER, True, True)
        self.assertEqual(report["gates"]["errors"]["result"], "FAIL")
        self.assertEqual(report["verdict"], "NO_GO")

    def test_errors_in_too_few_windows(self):
        rows = [r for r in passing_rows() if r["truth_state"] != "incorrect"] + [row("incorrect", start=i % 29 * 100_000) for i in range(600)]
        self.assertEqual(compute(rows, GOOD_KMER, True, True)["gates"]["spread"]["result"], "FAIL")

    def test_kmer_disagreement(self):
        report = compute(passing_rows(), Counter(agree=95, disagree=5), True, True)
        self.assertEqual(report["gates"]["independence"]["result"], "FAIL")

    def test_no_assessable_reads_is_not_a_pass(self):
        report = compute(passing_rows(), Counter(not_assessable=100), True, True)
        self.assertEqual(report["gates"]["independence"]["result"], "NOT_EVALUATED")

    def test_path_coverage(self):
        rows = passing_rows() + [row("path_unprojectable", kind="") for _ in range(2000)]
        self.assertEqual(compute(rows, GOOD_KMER, True, True)["gates"]["path_coverage"]["result"], "FAIL")

    def test_truth_coverage_counts_unaligned_reads_too(self):
        rows = passing_rows() + [row("unaligned", kind="", reason="no_single_alignment") for _ in range(1000)]
        self.assertEqual(compute(rows, GOOD_KMER, True, True)["gates"]["truth_coverage"]["result"], "FAIL")

    def test_missing_aligner_setting_fails_runs(self):
        rows = [r for r in passing_rows() if r["setting"] != "f"]
        self.assertEqual(compute(rows, GOOD_KMER, True, True)["gates"]["runs"]["result"], "FAIL")

    def test_per_aligner_failure_alone_does_not_stop(self):
        rows = passing_rows() + [row(aligner="panaligner", setting=s) for s in "abcdef"]
        rows = [r for r in rows if not (r["truth_state"] == "incorrect" and r["aligner"] == "panaligner")]
        rows += [row("incorrect", "minigraph", start=i * 100_000) for i in range(300)]
        report = compute(rows, GOOD_KMER, True, True)
        self.assertEqual(report["gates"]["per_aligner"]["result"], "FAIL")
        self.assertEqual(report["verdict"], "GO")


@unittest.skipIf(np is None, "numpy is not installed")
class KmerCheck(unittest.TestCase):
    def setUp(self):
        from ml.wg_truth import kmer_check
        self.module = kmer_check

    def test_reverse_complement_gives_same_canonical_kmer(self):
        forward = b"ACGTTGCAAGGCTTAACCGGATATCGCGTAG"
        complement = bytes.maketrans(b"ACGT", b"TGCA")
        reverse = forward.translate(complement)[::-1]
        a, valid_a = self.module.canonical_kmers(forward)
        b, valid_b = self.module.canonical_kmers(reverse)
        self.assertEqual(a.tolist(), b.tolist())
        self.assertTrue(valid_a.all() and valid_b.all())

    def test_ambiguous_base_invalidates_windows(self):
        _, valid = self.module.canonical_kmers(b"A" * 31 + b"N" + b"C" * 31)
        self.assertEqual(valid.tolist()[:2], [True, False])
        self.assertTrue(valid[-1])
        self.assertEqual(int(valid.sum()), 2)

    def test_short_sequence(self):
        self.assertEqual(len(self.module.canonical_kmers(b"ACGT")[0]), 0)

    def test_verdict(self):
        self.assertEqual(self.module.verdict([("c", 10)], "c", 0, 100), "not_assessable")
        self.assertEqual(self.module.verdict([("c", 10), ("c", 20), ("d", 1)], "c", 0, 100), "agree")
        self.assertEqual(self.module.verdict([("d", 10), ("d", 20), ("c", 1)], "c", 0, 100), "disagree")
        self.assertEqual(self.module.verdict([("c", 5000)] * 3, "c", 0, 100), "disagree")

    def test_scan_counts_occurrences_on_both_strands(self):
        rng = np.random.default_rng(1)
        sequence = bytes(rng.choice(list(b"ACGT"), 400).tolist())
        unique = sequence[10:41]
        repeated = sequence[100:131]
        complement = bytes.maketrans(b"ACGT", b"TGCA")
        reference = sequence + repeated.translate(complement)[::-1]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "r.fa"
            path.write_bytes(b">chrA\n" + reference + b"\n")
            kmer_a = self.module.canonical_kmers(unique)[0]
            kmer_b = self.module.canonical_kmers(repeated)[0]
            targets = np.unique(np.concatenate([kmer_a, kmer_b]))
            counts, places = self.module.scan_reference(path, targets)
        self.assertEqual(int(counts[np.searchsorted(targets, kmer_a[0])]), 1)
        self.assertEqual(int(counts[np.searchsorted(targets, kmer_b[0])]), 2)
        self.assertEqual(places[int(np.searchsorted(targets, kmer_a[0]))], ("chrA", 10))


class Inputs(unittest.TestCase):
    def test_read_records_uppercases_and_requires_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "a.fa"
            path.write_text(">x desc\nacg\nT\n>y\nGG\n", encoding="ascii")
            self.assertEqual(read_records(path, {"x"}), {"x": "ACGT"})
            with self.assertRaises(SystemExit):
                read_records(path, {"z"})

    def test_md5_is_content_only(self):
        self.assertEqual(sequence_md5("ACGT"), sequence_md5("ACGT"))
        self.assertNotEqual(sequence_md5("ACGT"), sequence_md5("ACGA"))

    def test_seed_differs_by_haplotype_and_setting(self):
        seeds = {setting_seed(11, a, l, h) for a in (0.95, 0.85, 0.75) for l in (10000, 2000) for h in (0, 1)}
        self.assertEqual(len(seeds), 12)


if __name__ == "__main__":
    unittest.main()
