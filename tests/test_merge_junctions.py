import os
import unittest
import gzip
from collections import OrderedDict

from tests import FakeArgs
from gdc_rnaseq_tools.merge_junctions import (
    StarJunctionRecord,
    load_junction_file,
    main,
)


class TestMergeStarJunctions(unittest.TestCase):
    star_junctions_1 = os.path.join(
        os.path.dirname(__file__), "etc/test_star_junctions_input_1.tsv.gz"
    )
    star_junctions_2 = os.path.join(
        os.path.dirname(__file__), "etc/test_star_junctions_input_2.tsv.gz"
    )
    exp_star_1 = os.path.join(
        os.path.dirname(__file__), "etc/exp_star_junctions_output_1.tsv.gz"
    )
    exp_star_1_2 = os.path.join(
        os.path.dirname(__file__), "etc/exp_star_junctions_output_1_2.tsv.gz"
    )
    out_test_pfx = os.path.join(
        os.path.dirname(__file__), "etc/test_star_merge_junctions_out"
    )
    to_remove = []

    def test_star_junction_record(self):
        """
        Tests the `StarJunctionRecord` class. 
        """
        line = "\t".join(["chr1", "100", "200", "1", "1", "1", "1", "0", "23"])
        rec = StarJunctionRecord.from_line(line)
        self.assertEqual("chr1", rec.chromosome)
        self.assertEqual(100, rec.intron_first)
        self.assertEqual(200, rec.intron_last)
        self.assertEqual(1, rec.strand)
        self.assertEqual(1, rec.motif)
        self.assertEqual(1, rec.annotation)
        self.assertEqual(1, rec.n_unique_mapped)
        self.assertEqual(0, rec.n_multi_mapped)
        self.assertEqual(23, rec.max_splice_overhang)
        self.assertEqual(("chr1", 100, 200, 1, 1, 1), rec.key)

    def test_star_junction_record_funcs(self):
        """
        Tests the `StarJunctionRecord` class functions. 
        """
        line1 = "\t".join(["chr1", "100", "200", "1", "1", "1", "1", "0", "23"])
        line2 = "\t".join(["chr1", "100", "200", "1", "1", "1", "3", "2", "20"])
        rec1 = StarJunctionRecord.from_line(line1)
        rec2 = StarJunctionRecord.from_line(line2)
        rec1 += rec2
        self.assertEqual("chr1", rec1.chromosome)
        self.assertEqual(100, rec1.intron_first)
        self.assertEqual(200, rec1.intron_last)
        self.assertEqual(1, rec1.strand)
        self.assertEqual(1, rec1.motif)
        self.assertEqual(1, rec1.annotation)
        self.assertEqual(4, rec1.n_unique_mapped)
        self.assertEqual(2, rec1.n_multi_mapped)
        self.assertEqual(23, rec1.max_splice_overhang)
        self.assertEqual(("chr1", 100, 200, 1, 1, 1), rec1.key)

        line3 = "\t".join(["chr1", "200", "300", "1", "1", "1", "1", "0", "23"])
        rec3 = StarJunctionRecord.from_line(line3)
        with self.assertRaises(AssertionError):
            rec1 += rec3

    def test_load_junction_file(self):
        """
        Tests the load_junction_file function.
        """
        dat1 = ["chr1", 100, 200, 1, 1, 1, 1, 0, 23]
        dat2 = ["chr10", 400, 700, 1, 1, 1, 1, 3, 13]
        dic = dict()
        dic = load_junction_file(self.star_junctions_1, dic)
        exp_keys = sorted([tuple(dat1[:6]), tuple(dat2[:6])])
        self.assertEqual(exp_keys, sorted(list(dic.keys())))

        found = dic[tuple(dat1[:6])]
        self.assertEqual(tuple(dat1[:6]), found.key)

        found = dic[tuple(dat2[:6])]
        self.assertEqual(tuple(dat2[:6]), found.key)

    def test_load_multi_junction_file(self):
        """
        Tests the load_junction_file function for multiple inputs.
        """
        dat1_1 = ("chr1", 100, 200, 1, 1, 1)
        dat1_2 = ("chr10", 400, 700, 1, 1, 1)
        dat2_1 = ("chr1", 100, 200, 1, 1, 1)
        dat2_2 = ("chr10", 400, 700, 1, 1, 1)
        dat2_3 = ("chr12", 100, 700, 1, 1, 1)

        dic = dict()
        dic = load_junction_file(self.star_junctions_1, dic)
        dic = load_junction_file(self.star_junctions_2, dic)

        exp_keys = sorted(set([dat1_1, dat1_2, dat2_1, dat2_2, dat2_3]))
        self.assertEqual(exp_keys, sorted(list(dic.keys())))

        self.assertEqual(2, dic[dat1_1].n_unique_mapped)
        self.assertEqual(0, dic[dat1_1].n_multi_mapped)
        self.assertEqual(23, dic[dat1_1].max_splice_overhang)

        self.assertEqual(3, dic[dat1_2].n_unique_mapped)
        self.assertEqual(6, dic[dat1_2].n_multi_mapped)
        self.assertEqual(23, dic[dat1_2].max_splice_overhang)

        self.assertEqual(1, dic[dat2_3].n_unique_mapped)
        self.assertEqual(3, dic[dat2_3].n_multi_mapped)
        self.assertEqual(23, dic[dat2_3].max_splice_overhang)

    def test_full_junction_single(self):
        """
        Tests the the whole main() function of junction merge for
        a single file.
        """
        args = FakeArgs()
        args.input = [self.star_junctions_1]
        args.output = self.out_test_pfx + ".1.tsv.gz"
        self.to_remove.append(args.output)
        main(args)
        with gzip.open(self.exp_star_1, "rt") as fh, gzip.open(
            args.output, "rt"
        ) as ofh:
            exp = fh.read()
            found = ofh.read()
            self.assertEqual(exp, found)
        os.remove(args.output)

    def test_full_junction_merge(self):
        """
        Tests from main() entry for single star file.
        """
        args = FakeArgs()
        args.input = [self.star_junctions_1, self.star_junctions_2]
        args.output = self.out_test_pfx + ".1_2.tsv.gz"
        self.to_remove.append(args.output)
        main(args)
        with gzip.open(self.exp_star_1_2, "rt") as fh, gzip.open(
            args.output, "rt"
        ) as ofh:
            exp = fh.read()
            found = ofh.read()
            self.assertEqual(exp, found)
        os.remove(args.output)

    def setUp(self):
        pass

    def tearDown(self):
        for fil in self.to_remove:
            if os.path.exists(fil):
                os.remove(fil)
