"""Read-origin truth for the whole-genome-graph pilot: interval projection, path projection and labels."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from ml.wg_truth.label_alignments import BedCover, judge, label_read, overlap_ratio
from ml.wg_truth.project_paths import PathProjection, project_path
from ml.wg_truth.project_truth import RefInterval, parse_paf_line, parse_sam_line, project_interval

REFS = {"chr1", "chr21"}


def paf(strand="+", qstart=0, qend=1000, qlen=1000, tstart=5000, cigar="1000M", tp="P", query="hap", target="chr21"):
    tend = tstart + 1000
    return f"{query}\t{qlen}\t{qstart}\t{qend}\t{strand}\t{target}\t100000\t{tstart}\t{tend}\t900\t1000\t60\ttp:A:{tp}\tcg:Z:{cigar}"


class IntervalProjection(unittest.TestCase):
    def test_forward_identity(self):
        record = parse_paf_line(paf())
        self.assertEqual(project_interval([record], 100, 300, "+"), (RefInterval("chr21", 5100, 5300, "+"), "ok"))

    def test_reverse_alignment_flips_interval_and_strand(self):
        record = parse_paf_line(paf(strand="-"))
        self.assertEqual(project_interval([record], 100, 300, "+"), (RefInterval("chr21", 5700, 5900, "-"), "ok"))
        self.assertEqual(project_interval([record], 100, 300, "-")[0].strand, "+")

    def test_deletion_in_haplotype_shifts_later_bases(self):
        record = parse_paf_line(paf(cigar="500M100D500M"))
        self.assertEqual(project_interval([record], 400, 600, "+")[0], RefInterval("chr21", 5400, 5700, "+"))

    def test_insertion_in_haplotype_shifts_later_bases(self):
        record = parse_paf_line(paf(cigar="500M100I400M"))
        self.assertEqual(project_interval([record], 400, 700, "+")[0], RefInterval("chr21", 5400, 5600, "+"))

    def test_end_inside_small_insertion_moves_inward(self):
        record = parse_paf_line(paf(cigar="500M100I400M"))
        self.assertEqual(project_interval([record], 200, 550, "+")[0], RefInterval("chr21", 5200, 5500, "+"))

    def test_end_inside_long_insertion_is_refused(self):
        record = parse_paf_line(paf(qend=3000, qlen=3000, cigar="500M2000I500M"))
        self.assertEqual(project_interval([record], 100, 1600, "+"), (None, "end_in_insertion"))

    def test_large_structural_difference_is_refused(self):
        record = parse_paf_line(paf(cigar="500M5000D500M"))
        self.assertEqual(project_interval([record], 400, 600, "+"), (None, "large_structural_difference"))

    def test_outside_alignment_or_across_two_is_refused(self):
        first = parse_paf_line(paf(qend=1000, qlen=3000))
        second = parse_paf_line(paf(qstart=1000, qend=2000, qlen=3000))
        self.assertEqual(project_interval([first, second], 900, 1100, "+"), (None, "no_single_alignment"))
        self.assertEqual(project_interval([first, second], 2500, 2700, "+"), (None, "no_single_alignment"))

    def test_two_containing_alignments_are_refused(self):
        record = parse_paf_line(paf())
        self.assertEqual(project_interval([record, record], 100, 300, "+"), (None, "multiple_alignments"))

    def test_paf_secondary_ignored_and_cg_required(self):
        self.assertIsNone(parse_paf_line(paf(tp="S")))
        with self.assertRaises(ValueError):
            parse_paf_line(paf().replace("cg:Z:1000M", "NM:i:0"))

    def test_reverse_with_query_offset(self):
        record = parse_paf_line(paf(strand="-", qstart=200, qend=1200, qlen=2000))
        self.assertEqual(record.query_start, 200)
        self.assertEqual(record.query_end, 1200)
        self.assertEqual(project_interval([record], 200, 300, "+")[0], RefInterval("chr21", 5900, 6000, "-"))

    def test_sam_matches_paf(self):
        sam_forward = parse_sam_line("hap\t0\tchr21\t5001\t60\t1000M\t*\t0\t0\t*\t*")
        self.assertEqual(project_interval([sam_forward], 100, 300, "+")[0], RefInterval("chr21", 5100, 5300, "+"))
        sam_reverse = parse_sam_line("hap\t16\tchr21\t5001\t60\t200S1000M800H\t*\t0\t0\t*\t*")
        paf_reverse = parse_paf_line(paf(strand="-", qstart=800, qend=1800, qlen=2000))
        self.assertEqual(sam_reverse.query_start, paf_reverse.query_start)
        self.assertEqual(project_interval([sam_reverse], 900, 1000, "+"), project_interval([paf_reverse], 900, 1000, "+"))

    def test_sam_secondary_ignored_supplementary_kept(self):
        self.assertIsNone(parse_sam_line("hap\t256\tchr21\t1\t0\t10M\t*\t0\t0\t*\t*"))
        self.assertIsNone(parse_sam_line("@SQ\tSN:chr21"))
        self.assertIsNotNone(parse_sam_line("hap\t2048\tchr21\t1\t0\t10M\t*\t0\t0\t*\t*"))


class PathProjectionTests(unittest.TestCase):
    def test_single_contig(self):
        self.assertEqual(project_path("-", "chr21", 10, 90, REFS), PathProjection("chr21", 10, 90, "-", "clean"))
        self.assertIsNone(project_path("+", "HG1#1#c1", 10, 90, REFS))

    def test_forward_reference_ends(self):
        path = ">chr21:100-200>chr21:200-300"
        self.assertEqual(project_path("+", path, 10, 190, REFS), PathProjection("chr21", 110, 290, "+", "clean"))

    def test_reverse_flips_strand(self):
        path = "<chr21:200-300<chr21:100-200"
        self.assertEqual(project_path("+", path, 10, 190, REFS), PathProjection("chr21", 110, 290, "-", "clean"))

    def test_inner_non_reference_is_clean(self):
        path = ">chr21:100-200>HG1#1#c1:0-50>chr21:300-400"
        self.assertEqual(project_path("+", path, 0, 250, REFS).kind, "clean")

    def test_path_ending_off_reference_is_flank_and_widened(self):
        path = ">chr21:100-200>HG1#1#c1:0-50"
        self.assertEqual(project_path("+", path, 10, 140, REFS), PathProjection("chr21", 110, 240, "+", "flank"))

    def test_path_starting_off_reference_is_flank_and_widened(self):
        path = ">HG1#1#c1:0-50>chr21:100-200"
        self.assertEqual(project_path("+", path, 20, 150, REFS), PathProjection("chr21", 70, 200, "+", "flank"))

    def test_reverse_flank(self):
        path = "<chr21:100-200<HG1#1#c1:0-50"
        self.assertEqual(project_path("+", path, 10, 140, REFS), PathProjection("chr21", 60, 190, "-", "flank"))

    def test_no_reference_segment(self):
        self.assertIsNone(project_path("+", ">HG1#1#c1:0-50>HG2#1#c2:0-50", 0, 100, REFS))

    def test_ends_on_two_chromosomes_or_orientations(self):
        self.assertIsNone(project_path("+", ">chr1:0-50>chr21:0-50", 0, 100, REFS))
        self.assertIsNone(project_path("+", ">chr21:0-50<chr21:100-150", 0, 100, REFS))

    def test_out_of_order_is_refused(self):
        self.assertIsNone(project_path("+", ">chr21:500-600>chr21:100-200", 0, 200, REFS))

    def test_clip_at_chromosome_start(self):
        path = ">HG1#1#c1:0-500>chr21:100-200"
        self.assertEqual(project_path("+", path, 0, 600, REFS).start, 0)


TRUTH = RefInterval("chr21", 1000, 2000, "+")


class Judge(unittest.TestCase):
    def test_overlap_ratio(self):
        self.assertEqual(overlap_ratio(0, 100, 200, 300), 0.0)
        self.assertEqual(overlap_ratio(0, 100, 50, 150), 50 / 150)

    def test_states(self):
        self.assertEqual(judge(TRUTH, PathProjection("chr21", 1000, 2000, "+", "clean"), 0.1), ("correct", ""))
        self.assertEqual(judge(TRUTH, PathProjection("chr1", 1000, 2000, "+", "clean"), 0.1), ("incorrect", "wrong_chromosome"))
        self.assertEqual(judge(TRUTH, PathProjection("chr21", 1000, 2000, "-", "clean"), 0.1), ("incorrect", "wrong_strand"))
        self.assertEqual(judge(TRUTH, PathProjection("chr21", 9000, 9900, "+", "clean"), 0.1), ("incorrect", "same_chromosome_wrong_place"))

    def test_threshold_boundary_is_exclusive(self):
        path = PathProjection("chr21", 1000, 1500, "+", "clean")  # ratio exactly 0.5
        self.assertEqual(judge(TRUTH, path, 0.5)[0], "incorrect")
        self.assertEqual(judge(TRUTH, path, 0.1)[0], "correct")


def gaf_fields(name, strand="+", path="chr21", path_start=5100, path_end=5300, matches=190):
    return [name, "200", "0", "200", strand, path, "100000", str(path_start), str(path_end), str(matches), "200", "60", "NM:i:10", "s1:i:100", "s2:i:10", "cm:i:20", "dv:f:0.05"]


class LabelRead(unittest.TestCase):
    NAME = "r1!hap!100!300!+"

    def setUp(self):
        self.alignments = {"hap": [parse_paf_line(paf())]}
        self.empty = BedCover(None)

    def label(self, candidates, name=None, ambiguous=None):
        return label_read(name or self.NAME, candidates, self.alignments, ambiguous or self.empty, self.empty, REFS)

    def test_correct_with_features(self):
        row = self.label([gaf_fields(self.NAME)])
        self.assertEqual(row["truth_state"], "correct")
        self.assertEqual(row["mapq"], 60)
        self.assertEqual(row["path_projected"], "clean")

    def test_unaligned(self):
        self.assertEqual(self.label([])["truth_state"], "unaligned")

    def test_incorrect_elsewhere(self):
        row = self.label([gaf_fields(self.NAME, path_start=50000, path_end=50200)])
        self.assertEqual((row["truth_state"], row["error_type"]), ("incorrect", "same_chromosome_wrong_place"))

    def test_longest_alignment_decides(self):
        short = gaf_fields(self.NAME, path_start=50000, path_end=50200, matches=10)
        self.assertEqual(self.label([short, gaf_fields(self.NAME)])["truth_state"], "correct")

    def test_truth_unprojectable_wins_over_alignment(self):
        name = "r1!other!100!300!+"
        self.assertEqual(self.label([gaf_fields(name)], name=name)["truth_state"], "truth_unprojectable")

    def test_path_unprojectable(self):
        row = self.label([gaf_fields(self.NAME, path="HG1#1#c1", path_start=0, path_end=200)])
        self.assertEqual(row["truth_state"], "path_unprojectable")

    def test_strict_state_differs(self):
        row = self.label([gaf_fields(self.NAME, path_start=5240, path_end=5440)])  # overlap 60 of union 340
        self.assertEqual((row["truth_state"], row["truth_state_strict"]), ("correct", "incorrect"))

    def test_ambiguity_uses_source_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "a.bed"
            bed.write_text("hap\t100\t280\nhap\t500\t600\n", encoding="utf-8")
            row = self.label([gaf_fields(self.NAME)], ambiguous=BedCover(bed))
        self.assertEqual(row["truth_state"], "truth_ambiguous")
        self.assertEqual(row["ambiguous_fraction"], 0.9)

    def test_partly_difficult_read_is_still_judged(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "a.bed"
            bed.write_text("hap\t100\t200\n", encoding="utf-8")
            self.assertEqual(self.label([gaf_fields(self.NAME)], ambiguous=BedCover(bed))["truth_state"], "correct")


class BedCoverTests(unittest.TestCase):
    def test_merges_overlaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "a.bed"
            bed.write_text("c\t0\t10\nc\t5\t20\nc\t30\t40\n", encoding="utf-8")
            cover = BedCover(bed)
        self.assertEqual(cover.fraction("c", 0, 40), 30 / 40)
        self.assertEqual(cover.fraction("c", 25, 35), 0.5)
        self.assertEqual(cover.fraction("absent", 0, 10), 0.0)


if __name__ == "__main__":
    unittest.main()
