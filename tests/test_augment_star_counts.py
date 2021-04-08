import gzip
import os
import unittest
from collections import OrderedDict

import pandas as pd

from gdc_rnaseq_tools.augment_star_counts import (
    augment,
    calc_fpkm,
    calc_fpkm_uq,
    calc_tpm,
    get_extras,
    load_table,
    main,
    merge_tables,
    save_result,
    validate_table,
)
from gdc_rnaseq_tools.utils import DataError, DataFormatError, get_logger
from tests import FakeArgs


class TestAugmentStarCounts(unittest.TestCase):

    df1_raw = pd.DataFrame(
        {
            'gene_id': [
                'ENSG00000000001.0',
                'ENSG00000000002.0',
                'ENSG00000000003.0',
                'ENSG00000000004.0',
                'ENSG00000000005.0',
                'ENSG00000000006.0',
                'ENSG00000000007.0',
            ],
            'total_exon_length': [1000, 1000, 1000, 1000, 1000, 1000, 1000],
            'gene_name': ['ONE', 'TWO', 'THREE', 'FOUR', 'FIVE', 'SIX', 'SEVEN'],
            'gene_type': [
                'protein_coding',
                'protein_coding',
                'protein_coding',
                'protein_coding',
                'protein_coding',
                'protein_coding',
                'protein_coding',
            ],
            'chromosome': ['chrX', 'chr1', 'chr22', 'chr3', 'chr5', 'chrM', 'chr10'],
            'unstranded': [1500, 0, 1500, 6000, 7500, 6000, 7500],
            'stranded_first': [10, 10, 10, 10, 10, 10, 10],
            'stranded_second': [10, 10, 10, 10, 10, 10, 10],
        }
    )

    df1_tpm = pd.Series([50000.0, 0.0, 50000.0, 200000.0, 250000.0, 200000.0, 250000.0])

    df1_fpkm = pd.Series(
        [50000.0, 0.0, 50000.0, 200000.0, 250000.0, 200000.0, 250000.0]
    )

    df1_fpkm_uq = pd.Series(
        [40000.0, 0.0, 40000.0, 160000.0, 200000.0, 160000.0, 200000.0]
    )

    ts1_counts_file = os.path.join(
        os.path.dirname(__file__), "etc/test_set_1.counts.tsv.gz"
    )
    ts1_gene_info_file = os.path.join(
        os.path.dirname(__file__), "etc/test_set_1.gene_info.tsv.gz"
    )
    ts1_annotated_file = os.path.join(
        os.path.dirname(__file__), "etc/test_set_1.annotated.tsv.gz"
    )
    ts1_final_file = os.path.join(
        os.path.dirname(__file__), "etc/test_set_1.final.tsv.gz"
    )

    simple_file = os.path.join(os.path.dirname(__file__), "etc/simple_counts.tsv.gz")

    to_remove = []

    logger = get_logger("augment_counts_table.testing")

    def test_calc_fpkm(self):
        """
        Tests the `calc_fpkm` function
        """
        fpkm = calc_fpkm(
            expression=self.df1_raw.unstranded,
            feature_effective_length=self.df1_raw.total_exon_length,
            gene_type=self.df1_raw.gene_type,
        )

        self.assertTrue((fpkm == self.df1_fpkm).all())

    def test_calc_fpkm_uq(self):
        """
        Tests the `calc_fpkm_uq` function
        """
        fpkm_uq = calc_fpkm_uq(
            expression=self.df1_raw.unstranded,
            feature_effective_length=self.df1_raw.total_exon_length,
            gene_type=self.df1_raw.gene_type,
            chromosome=self.df1_raw.chromosome,
        )

        self.assertTrue((fpkm_uq == self.df1_fpkm_uq).all())

    def test_calc_tpm(self):
        """
        Tests the `calc_tpm` function
        """
        tpm = calc_tpm(
            expression=self.df1_raw.unstranded,
            feature_effective_length=self.df1_raw.total_exon_length,
        )

        # print(tpm)
        # print(self.df1_tpm)
        self.assertTrue((tpm == self.df1_tpm).all())

    def test_load_table(self):
        """
        Tests the `load_tables` function
        """
        df1 = load_table(self.simple_file)
        pd.testing.assert_frame_equal(df1, self.df1_raw)

    def test_validate_table_good(self):
        """
        Tests the `validate_tables` function
        """
        validate_table(self.df1_raw, self.df1_raw.columns.tolist())

    def test_validate_table_bad(self):
        """
        Tests the `validate_tables` function with incorrect column names
        """
        self.assertRaises(
            DataFormatError,
            validate_table,
            df=self.df1_raw,
            expected_columns=['not', 'in', 'the', 'table'],
        )

    def test_merge(self):
        """
        Tests the `merge_tables` function
        """
        df1 = pd.DataFrame({'id': [1, 2, 3, 4], 'A': ['one', 'two', 'three', 'four'],})
        df2 = pd.DataFrame({'id': [2, 3, 4, 5], 'B': ['owt', 'eerht', 'ruof', 'evif']})
        res = pd.DataFrame(
            {
                'id': [2, 3, 4],
                'A': ['two', 'three', 'four'],
                'B': ['owt', 'eerht', 'ruof'],
            }
        )
        pd.testing.assert_frame_equal(merge_tables(df1, df2, on='id'), res)

    def test_get_extras(self):
        """
        Tests the `get_extras` function
        """
        test_answer = self.df1_raw.iloc[0:4].copy()
        res = get_extras(self.df1_raw)
        pd.testing.assert_frame_equal(res, test_answer)

    def test_save_result(self):
        """
        Test the `save_results` function
        """
        outfile = 'save_results_output.tsv'
        self.to_remove.append(outfile)

        testdf = pd.DataFrame(
            {'id': [1, 2, 3, 4], 'A': ['one', 'two', 'three', 'four'],}
        )

        save_result(df=testdf, outfile=outfile, gencode_version=36)

        with open(outfile, 'rt') as result:
            res_lines = result.read()

        expected = (
            '# gene-model: GENCODE v36\n'
            'id\tA\n'
            '1\tone\n'
            '2\ttwo\n'
            '3\tthree\n'
            '4\tfour\n'
        )
        self.assertEqual(expected, res_lines)

    def test_augment(self):
        """
        Tests `augment` function
        """
        outfile = 'augment_output.tsv'
        self.to_remove.append(outfile)

        augment(
            counts_file=self.ts1_counts_file,
            gene_info_file=self.ts1_gene_info_file,
            outfile=outfile,
            gencode_version=36,
            logger=self.logger,
        )

    def test_full_run(self):
        """
        Full end-to-end test
        """
        outfile = 'main_output.tsv'
        self.to_remove.append(outfile)

        args = FakeArgs()
        args.input = self.ts1_counts_file
        args.gene_info = self.ts1_gene_info_file
        args.output = outfile
        args.gencode_version = 36
        self.to_remove.append(args.output)

        main(args)

        result = pd.read_table(args.output, comment='#')
        expected = pd.read_table(self.ts1_final_file, comment='#')
        pd.testing.assert_frame_equal(result, expected)

    def setUp(self):
        pass

    def tearDown(self):
        for fil in self.to_remove:
            if os.path.exists(fil):
                os.remove(fil)
