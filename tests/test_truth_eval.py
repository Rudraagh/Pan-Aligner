import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from simulate_truth_reads import build_liftover, parse_truth_name, read_fasta, write_fasta_with_fai
from truth_eval import classify_reads, mapq_curve, overlaps, parse_mapeval, project_stable_path, summarize


REF = "ref|NC_000001.11|:226870184-226927726"


def gaf_line(read, strand, path, path_len, start, end, matches, mapq, tags=("tp:A:P", "s1:i:100", "s2:i:0")):
    return "\t".join([read, "10000", "0", "10000", strand, path, str(path_len), str(start), str(end),
                      str(matches), "10000", str(mapq), *tags]) + "\n"


class LiftoverTests(unittest.TestCase):
    def test_equal_length_haplotype_maps_one_to_one(self):
        lift = build_liftover("ACGTACGT", "ACCTACGT")
        self.assertEqual([lift(x) for x in (0, 3, 8)], [0, 3, 8])

    def test_insertion_shifts_coordinates_after_the_window(self):
        reference = "AAAACCCCGGGG"
        haplotype = "AAAACCTTCCGGGG"  # 2 bp insertion inside the C run
        lift = build_liftover(reference, haplotype)
        self.assertEqual(lift(2), 2)
        self.assertEqual(lift(len(haplotype)), len(reference))
        self.assertEqual(lift(12), 10)

    def test_deletion_shifts_the_other_way(self):
        lift = build_liftover("AAAACCCCGGGG", "AAAACCGGGG")
        self.assertEqual(lift(10), 12)
        self.assertEqual(lift(1), 1)

    def test_widely_separated_indels_are_refused(self):
        reference = "A" * 10 + "C" * 200 + "G" * 10
        haplotype = "A" * 10 + "T" + "C" * 200 + "G" * 10 + "T"  # changes at both ends, 2 bp longer
        with self.assertRaises(ValueError):
            build_liftover(reference, haplotype)


class FastaTests(unittest.TestCase):
    def test_fai_offsets_point_at_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.fa"
            write_fasta_with_fai([("a", "ACGT"), ("bb", "GG")], path)
            data = path.read_bytes()
            for line in Path(f"{path}.fai").read_text().splitlines():
                name, length, offset, _, _ = line.split("\t")
                self.assertEqual(data[int(offset):int(offset) + int(length)].decode(), dict(read_fasta(path))[name])

    def test_truth_name_with_pipes_and_colons(self):
        self.assertEqual(parse_truth_name(f"S1_2!{REF}!1877!11834!-"), ("S1_2", REF, 1877, 11834, "-"))


class ScoringTests(unittest.TestCase):
    def test_overlap_rule_matches_mapeval(self):
        self.assertTrue(overlaps((0, 100), 80, 180))      # 20/180 > 0.1
        self.assertFalse(overlaps((0, 100), 95, 195))     # 5/195 < 0.1
        self.assertFalse(overlaps((0, 100), 100, 200))

    def test_projection_of_forward_and_reverse_stable_paths(self):
        path = f">{REF}:0-11877>PSEN2_H_9:11877-11878>{REF}:11878-57543"
        self.assertEqual(project_stable_path(path, 57543, 15, 57541, REF), (15, 57541, False))
        reverse = f"<{REF}:11878-20000<PSEN2_H_9:11877-11878<{REF}:5000-11877"
        length = 8122 + 1 + 6877
        self.assertEqual(project_stable_path(reverse, length, 100, length - 50, REF), (5050, 19900, True))
        self.assertIsNone(project_stable_path(f">PSEN2_H_9:0-10>{REF}:10-500", 500, 0, 400, REF))
        self.assertIsNone(project_stable_path(f">{REF}:0-10<{REF}:10-500", 500, 0, 400, REF))
        # 1 bp inversion bubble inside a reverse path (seen in the PSEN1 graph): ends decide
        inverted = f"<{REF}:61643-69967>{REF}:61642-61643<{REF}:56445-61642"
        self.assertEqual(project_stable_path(inverted, 13522, 87, 9212, REF), (56445 + 4310, 69967 - 87, True))

    def test_classification_counts_and_mapeval_parity_rules(self):
        reads = [f"r1!{REF}!1000!11000!+", f"r2!{REF}!1000!11000!+", f"r3!{REF}!30000!40000!-", f"r4!{REF}!0!10000!+"]
        lines = [
            gaf_line(reads[0], "+", REF, 57543, 1010, 10990, 9000, 60),
            gaf_line(reads[0], "+", REF, 57543, 40000, 45000, 100, 60, tags=("tp:A:S", "s1:i:5")),  # skipped like mapeval
            gaf_line(reads[1], "+", REF, 57543, 40000, 50000, 9000, 3),                              # wrong place
            gaf_line(reads[2], "-", f">{REF}:0-35000>PSEN2_H_9:35000-35001>{REF}:35001-57543", 57544, 30010, 39990, 9000, 60),
        ]                                                                                             # r4 unaligned
        with tempfile.TemporaryDirectory() as tmp:
            gaf = Path(tmp) / "x.gaf"
            gaf.write_text("".join(lines), encoding="utf-8")
            outcomes = classify_reads(reads, gaf, REF)
        summary = summarize(outcomes)
        self.assertEqual((summary["aligned"], summary["unaligned"]), (3, 1))
        self.assertEqual((summary["mapeval_correct"], summary["mapeval_wrong"]), (1, 2))
        self.assertEqual((summary["projected_correct"], summary["projected_wrong"]), (2, 1))
        self.assertEqual(summary["off_reference_path"], 1)
        curve = {row["mapq_threshold"]: row for row in mapq_curve(outcomes)}
        self.assertEqual((curve[0]["aligned"], curve[0]["wrong"]), (3, 2))
        self.assertEqual((curve[5]["aligned"], curve[5]["wrong"]), (2, 1))

    def test_parse_mapeval_rows(self):
        rows = parse_mapeval("Q\t60\t12\t0\t0.000000000\t12\nQ\t3\t2\t1\t0.071428571\t14\n")
        self.assertEqual([(r["mapq"], r["reads"], r["wrong"], r["cumulative_reads"]) for r in rows], [(60, 12, 0, 12), (3, 2, 1, 14)])


if __name__ == "__main__":
    unittest.main()
