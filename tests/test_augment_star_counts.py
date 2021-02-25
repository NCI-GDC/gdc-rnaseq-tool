import gzip
import os
import unittest
from collections import OrderedDict

import pandas as pd

from gdc_rnaseq_tools.augment_star_counts import (
    calc_fpkm,
    calc_fpkm_uq,
    calc_tpm,
    load_counts_table,
    load_gene_info,
    main,
    prep_data,
)
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
            'Chromosome': ['chrX', 'chr1', 'chr22', 'chr3', 'chr5', 'chrM', 'chr10'],
            'stranded_first': [10, 10, 10, 10, 10, 10, 10],
            'stranded_second': [10, 10, 10, 10, 10, 10, 10],
            'unstranded': [1500, 0, 1500, 6000, 7500, 6000, 7500],
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

    to_remove = []

    def test_calc_fpkm(self):
        """
        Tests the `calc_fpkm` function
        """
        fpkm = calc_fpkm(self.df1_raw)

        # print(fpkm)
        # print(self.df1_fpkm)
        self.assertTrue((fpkm == self.df1_fpkm).all())

    def test_calc_fpkm_uq(self):
        """
        Tests the `calc_fpkm_uq` function
        """
        fpkm_uq = calc_fpkm_uq(self.df1_raw)

        # print(fpkm_uq)
        # print(self.df1_fpkm_uq)
        self.assertTrue((fpkm_uq == self.df1_fpkm_uq).all())

    def test_calc_tpm(self):
        """
        Tests the `calc_tpm` function
        """
        tpm = calc_tpm(self.df1_raw)

        # print(tpm)
        # print(self.df1_tpm)
        self.assertTrue((tpm == self.df1_tpm).all())

    def test_merge(self):
        """
        Tests the prep_data function
        """
        df, extradf = prep_data(self.ts1_counts_file, self.ts1_gene_info_file)
        df = pd.concat([extradf, df]).reset_index(drop=True)

        mdf = pd.read_table(self.ts1_annotated_file)
        pd.testing.assert_frame_equal(df, mdf)

    def test_full_run(self):
        """
        Full end-to-end test
        """
        args = FakeArgs()
        args.input = self.ts1_counts_file
        args.gene_info = self.ts1_gene_info_file
        args.output = "counts_report.tsv"
        args.pragma_line = "# gene-model: GENCODE v36"
        self.to_remove.append(args.output)

        main(args)

        df = pd.read_table(args.output, comment='#')
        tdf = pd.read_table(self.ts1_final_file, comment='#')
        pd.testing.assert_frame_equal(df, tdf)

    def setUp(self):
        pass

    def tearDown(self):
        for fil in self.to_remove:
            if os.path.exists(fil):
                os.remove(fil)
