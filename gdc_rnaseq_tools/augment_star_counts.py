import logging
from argparse import Namespace
from enum import Enum
from typing import List, Optional, Text, Union

import numpy as np
import pandas as pd

from gdc_rnaseq_tools.utils import DataError, DataFormatError, get_logger
from tests import FakeArgs


class ColumnNames(Enum):
    @classmethod
    def cols(cls) -> List[Text]:
        '''
        Returns an ordered list of Enum members. Designed to return the
        columns in a specific order for purposes of reading or printing.

        Args:
            None
        
        Returns:
            list of column names as strings
        '''
        return [member.value for member in cls]


class CountsColumns(ColumnNames):
    GENE_ID = 'gene_id'
    UNSTRANDED = 'unstranded'
    STRANDED_FIRST = 'stranded_first'
    STRANDED_SECOND = 'stranded_second'


class GeneInfoColumns(ColumnNames):
    GENE_ID = 'gene_id'
    TOTAL_EXON_LENGTH = 'total_exon_length'
    GENE_NAME = 'gene_name'
    GENE_TYPE = 'gene_type'
    CHROMOSOME = 'chromosome'


class MergedColumns(ColumnNames):
    GENE_ID = 'gene_id'
    TOTAL_EXON_LENGTH = 'total_exon_length'
    GENE_NAME = 'gene_name'
    GENE_TYPE = 'gene_type'
    CHROMOSOME = 'chromosome'
    UNSTRANDED = 'unstranded'
    STRANDED_FIRST = 'stranded_first'
    STRANDED_SECOND = 'stranded_second'


class FinalColumns(ColumnNames):
    GENE_ID = 'gene_id'
    GENE_NAME = 'gene_name'
    GENE_TYPE = 'gene_type'
    UNSTRANDED = 'unstranded'
    STRANDED_FIRST = 'stranded_first'
    STRANDED_SECOND = 'stranded_second'
    TPM_UNSTRANDED = 'tpm_unstranded'
    FPKM_UNSTRANDED = 'fpkm_unstranded'
    FPKM_UQ_UNSTRANDED = 'fpkm_uq_unstranded'


def load_table(
    table_filename: Text, colnames: Optional[List[Text]] = None
) -> pd.DataFrame:
    '''
    Loads tabular data into a DataFrame. If 

    Args:
        table_filename: file name of the tabular data
        colnames: a list of column names to be used when the table has no 
            column headers.
    
    Returns:
        pandas DataFrame
    '''

    return pd.read_table(table_filename, names=colnames, comment='#')


def validate_table(df: pd.DataFrame, expected_columns: List[Text]) -> None:
    '''
    Verifies that the table has exactly the expected columns

    Args:
        df: the data frame to check
        expected_columns: the list of columns
    
    Returns:
        None
    
    Throws
    '''
    if set(expected_columns) - set(df.columns) != set():
        raise DataFormatError('Expected columns not found')


def merge_tables(df1: pd.DataFrame, df2: pd.DataFrame, on: Text) -> pd.DataFrame:
    '''
    Performs an inner-join on data frames

    Args:
        df1: left data frame
        df2: right data frame
        on: common column name used to join data frames
    
    Returns:
        A pandas.DataFrame created by joining on the key
    '''

    return pd.merge(df1, df2, on=on, how='inner')


def get_extras(df: pd.DataFrame) -> pd.DataFrame:
    '''
    STAR counts have 4 extra lines at the top reporting unmapped,
    multimapping, noFeature, and ambiguous reads. These are excluded
    by the join operation for calculations but must be added back in
    to the final data file. This function extracts those rows.

    Args:
        df: the counts data frame
    
    Returns:
        pandas.DataFrame containing misaligned reads stats
    '''

    return df.iloc[0:4].copy()


def calc_tpm(expression: pd.Series, feature_effective_length: pd.Series) -> pd.Series:
    '''
    Transcripts Per Million

    Calculate the TPM values for all genes

        TPM = RPK x 1e6 / M

            RPK = C X 1e3 / L
            C : count of fragments aligned to this gene
            L : sum of exon lengths in gene where overlapping exons are merged
                otherwise known as union exon length
            M : sum over all genes of RPK values
    
    Args:
        expression: raw counts of aligned reads
        feature_effective_length: lengths of unified exons of each gene

    Returns:
        pandas.Series containing the calculated TPM values
    '''
    # RPK - reads per thousand bp of transcript length
    rpk = expression * 1e3 / feature_effective_length
    # sum of RPK signal
    M = rpk.sum()
    # TPM - transcripts per million
    tpm = rpk * 1e6 / M
    return tpm


def calc_fpkm(
    expression: pd.Series, feature_effective_length: pd.Series, gene_type: pd.Series
) -> pd.Series:
    '''
    Fragments Per Kilobases (of transcript) and Millions (of fragments)

    Calculate the FPKM values for all genes

        FPKM = (C × 1e9)/(NL)

            C : count of fragments aligned to this gene
            N : total fragment count to protein-coding genes
            L : sum of exon lengths in gene where overlapping exons are merged
    
    Args:
        expression: raw counts of aligned reads
        feature_effective_length: lengths of unified exons of each gene 
        gene_type: gene biotypes used for calculating sum of expression of 
            protein coding genes
    '''
    # select protein coding genes
    sel = gene_type == 'protein_coding'
    # get sum of counts in protein coding genes
    N = expression.loc[sel].sum()
    # calculate fpkm
    fpkm = expression * 1e9 / (N * feature_effective_length)

    return fpkm


def calc_fpkm_uq(
    expression: pd.Series,
    feature_effective_length: pd.Series,
    gene_type: pd.Series,
    chromosome: pd.Series,
) -> pd.Series:
    '''
    Upper Quartile normalized FPKM

        FPKM-UQ = (C × 1e9)/(UGL)

            C : count of fragments aligned to this gene
            U : upper quartile of fragment counts to autosomal protein-coding
                genes with count > 0
            G : number of protein-coding genes on autosomes
            L : sum of exon lengths in gene where overlapping exons are merged
    
    Args:
        expression: raw counts of aligned reads
        feature_effective_length: lengths of unified exons of each gene 
        gene_type: gene biotypes used for calculating sum of expression of 
            protein coding genes
        chromosome: chromosome name on which gene is found
        
    '''
    # selections for U and G
    sel_prot = gene_type == 'protein_coding'
    sel_autosomes = ~chromosome.isin(['chrX', 'chrY', 'chrM'])
    sel_nonzero = expression > 0
    # combine selections
    sel_U = sel_prot & sel_autosomes & sel_nonzero
    sel_G = sel_prot & sel_autosomes

    # Calculate U and G
    U = np.quantile(expression.loc[sel_U], 0.75)
    G = len(expression[sel_G])

    # Calculate FPKM-UQ
    fpkm_uq = expression * 1e9 / (U * G * feature_effective_length)
    return fpkm_uq


def save_result(
    df: pd.DataFrame, outfile: Text, gencode_version: Optional[int] = None
) -> None:
    '''
    Write output table as TSV with 4 places of floating point precision

    Args:
        df: final results table
        outfile: output file name
        pragma_line: informational line to be added to top of output file
    '''

    with open(outfile, 'w') as out:
        if gencode_version is not None:
            out.write('# gene-model: GENCODE v{}\n'.format(gencode_version))
        df.to_csv(out, sep='\t', header=True, index=False, float_format='%.4f')


def augment(
    counts_file: Text,
    gene_info_file: Text,
    outfile: Text,
    gencode_version: int,
    logger: logging.Logger,
) -> None:
    '''
    Augment STAR read counts with normalized counts and gene info

    Adds gene names and gene bio-types from data extracted from a GENCODE gtf
    Adds TPM, FPKM, FPKM-UQ normalized counts

    Writes an output file as TSV with a header line and pragma line indicating
        the genome annotation version

    Args:
        counts_file: file name for STAR counts file
        gene_info_file: file name for gene info containing 
            [ gene_id, total_exon_length, gene_name, gene_type, chromosome ]
        outfile: output file name
        pragma_line: free-text string to be added to top of results file
        logger: logging.Logger object used to communicate messages
    '''

    # load data
    logger.info("Reading counts file {}".format(counts_file))
    counts = load_table(counts_file, CountsColumns.cols())
    validate_table(counts, CountsColumns.cols())

    logger.info("Reading gene info file {}".format(gene_info_file))
    gene_info = load_table(gene_info_file)
    validate_table(gene_info, GeneInfoColumns.cols())

    # merge counts with gene info
    logger.info("Merging counts and gene info tables")
    merged = merge_tables(gene_info, counts, on=GeneInfoColumns.GENE_ID.value)
    validate_table(merged, MergedColumns.cols())

    # calculate new normalized counts
    logger.info("Calculating normalized counts")
    # FPKM
    merged[FinalColumns.FPKM_UNSTRANDED.value] = calc_fpkm(
        expression=merged[MergedColumns.UNSTRANDED.value],
        feature_effective_length=merged[MergedColumns.TOTAL_EXON_LENGTH.value],
        gene_type=merged[MergedColumns.GENE_TYPE.value],
    )

    # FPKM-UQ
    merged[FinalColumns.FPKM_UQ_UNSTRANDED.value] = calc_fpkm_uq(
        expression=merged[MergedColumns.UNSTRANDED.value],
        feature_effective_length=merged[MergedColumns.TOTAL_EXON_LENGTH.value],
        gene_type=merged[MergedColumns.GENE_TYPE.value],
        chromosome=merged[MergedColumns.CHROMOSOME.value],
    )

    # TPM
    merged[FinalColumns.TPM_UNSTRANDED.value] = calc_tpm(
        expression=merged[MergedColumns.UNSTRANDED.value],
        feature_effective_length=merged[MergedColumns.TOTAL_EXON_LENGTH.value],
    )

    # add back extra alignment stats
    misalign_stats = get_extras(counts)
    final = pd.concat([misalign_stats, merged], axis=0)
    final = final[FinalColumns.cols()].copy()

    # write output table
    logger.info("Saving results to {}".format(outfile))
    save_result(df=final, outfile=outfile, gencode_version=gencode_version)


def main(args: Union[FakeArgs, Namespace]) -> None:
    """
    Main entrypoint for augment_counts_table. Maps CLI args to function args

    Args:
        args: generally an argparse.Namespace object, but FakeArgs for testing
    """
    logger = get_logger("augment_counts_table")
    logger.info("Augmenting STAR gene counts file {}.".format(args.input))

    augment(
        counts_file=args.input,
        gene_info_file=args.gene_info,
        outfile=args.output,
        gencode_version=args.gencode_version,
        logger=logger,
    )
