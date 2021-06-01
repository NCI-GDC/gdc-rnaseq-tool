import gzip
import os
import unittest
from collections import OrderedDict

from gdc_rnaseq_tools.merge_counts import load_star_file, main, merge_star_counts
from tests import FakeArgs


class TestMergeStarGeneCounts(unittest.TestCase):
    star_counts_1 = os.path.join(
        os.path.dirname(__file__), "etc/test_star_counts_input_1.tsv.gz"
    )
    star_counts_2 = os.path.join(
        os.path.dirname(__file__), "etc/test_star_counts_input_2.tsv.gz"
    )
    exp_star_1 = os.path.join(
        os.path.dirname(__file__), "etc/exp_star_counts_output_1.tsv.gz"
    )
    exp_star_1_2 = os.path.join(
        os.path.dirname(__file__), "etc/exp_star_counts_output_1_2.tsv.gz"
    )
    out_test_pfx = os.path.join(
        os.path.dirname(__file__), "etc/test_star_merge_gene_counts_out"
    )
    to_remove = []

    def test_merge_star_counts(self):
        """
        Tests the `merge_star_counts` function.
        """
        dic = OrderedDict()
        dic["B"] = [[1, 2, 3], [2, 1, 5]]
        dic["A"] = [[2, 2, 2], [3, 3, 3]]
        giter = merge_star_counts(dic)

        key, counts = next(giter)
        self.assertEqual("B", key)
        self.assertEqual([3, 3, 8], counts)

        key, counts = next(giter)
        self.assertEqual("A", key)
        self.assertEqual([5, 5, 5], counts)

        with self.assertRaises(StopIteration):
            key, counts = next(giter)

    def test_load_star_file(self):
        """
        Tests loading the staged STAR counts file.
        """
        dic = OrderedDict()
        dic = load_star_file(self.star_counts_1, dic)
        expected = OrderedDict()
        expected["ZZZZ"] = [[100, 0, 50]]
        expected["AAAA"] = [[0, 50, 75]]
        expected["CCCC"] = [[10, 10, 10]]
        self.assertEqual(expected, dic)

        dic = load_star_file(self.star_counts_2, dic)
        expected = OrderedDict()
        expected["ZZZZ"] = [[100, 0, 50], [10, 20, 50]]
        expected["AAAA"] = [[0, 50, 75], [0, 5, 705]]
        expected["CCCC"] = [[10, 10, 10], [10, 0, 20]]
        self.assertEqual(expected, dic)

    def test_full_single(self):
        """
        Tests from main() entry for single star file.
        """
        args = FakeArgs()
        args.input = [self.star_counts_1]
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

    def test_full_merge(self):
        """
        Tests from main() entry for single star file.
        """
        args = FakeArgs()
        args.input = [self.star_counts_1, self.star_counts_2]
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
