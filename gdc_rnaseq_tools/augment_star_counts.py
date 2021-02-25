import numpy as np
import pandas as pd

from gdc_rnaseq_tools.utils import DataError, DataFormatError, get_logger


def load_counts_table(counts):
    '''
    Loads per-gene counts data table as created by STAR


    returns pandas DataFrame with the following columns:
        gene_id
        stranded_first
        stranded_second
        unstranded
    '''
    cdf = pd.read_table(
        counts,
        header=None,
        names=['gene_id', 'stranded_first', 'stranded_second', 'unstranded'],
    )
    return cdf


def load_gene_info(gene_info):
    '''
    Loads table with additional information used to either annotate or 
    calculate values needed to augment counts table.

    file should be TSV and should have the following columns:
        gene_id
        total_exon_length
        gene_name
        gene_type
        Chromosome

    returns pandas DataFrame with all columns
    '''
    # load data
    gdf = pd.read_table(gene_info)

    # check that data has exactly the columns expected
    expected_cols = [
        'gene_id',
        'total_exon_length',
        'gene_name',
        'gene_type',
        'Chromosome',
    ]
    if set(expected_cols) - set(gdf.columns) != set():
        raise DataFormatError('Expected columns not found')

    # would be good to check that file is the correct length, but that's
    # tied to the annotation and I'm not sure how to do that
    return gdf


def prep_data(counts, gene_info):
    '''
    Load data files and merge on gene_id via "inner" join operation.

    STAR counts have 4 extra lines at the top reporting unmapped,
    multimapping, noFeature, and ambiguous reads. These are not included
    in the merged DataFrame, but are returned in the extras DataFrame to
    be merged back into the final output.

    Returns a tuple of:
        merged DataFrame
        extras DataFrame
    '''
    cdf = load_counts_table(counts)
    gdf = load_gene_info(gene_info)
    df = pd.merge(gdf, cdf, on='gene_id', how='inner')

    if len(df) != len(gdf):
        raise DataError(
            'Data length or key mismatch between gene counts and gene '
            + 'annotation tables. Possible file truncation.'
        )

    # save table of mis-aligned read stats also reported by star
    mdf = cdf.iloc[0:4].copy()

    return df, mdf


def calc_tpm(df, excol='unstranded', lencol='total_exon_length'):
    '''
    Transcripts Per Million

    Calculate the TPM values for all genes in DataFrame df

    TPM = RPK x 1e6 / M

        RPK = C X 1e3 / L
        C : count of fragments aligned to this gene
        L : sum of exon lengths in gene where overlapping exons are merged
            otherwise known as union exon length
        M : sum over all genes of RPK values
    '''
    # RPK - reads per thousand bp of transcript length
    rpk = df[excol] * 1e3 / df[lencol]
    # sum of RPK signal
    M = rpk.sum()
    # TPM - transcripts per million
    tpm = rpk * 1e6 / M
    return tpm


def calc_fpkm(df, excol='unstranded', lencol='total_exon_length', typecol='gene_type'):
    '''
    Fragments Per Kilobases (of transcript) and Millions (of fragments)

    Calculate the FPKM values for all genes in DataFrame df

    FPKM = (C × 1e9)/(NL)

        C : count of fragments aligned to this gene
        N : total fragment count to protein-coding genes
        L : sum of exon lengths in gene where overlapping exons are merged
    '''
    # select protein coding genes
    sel = df[typecol] == 'protein_coding'
    # get sum of counts in protein coding genes
    N = df.loc[sel][excol].sum()
    # calculate fpkm
    fpkm = df[excol] * 1e9 / (N * df[lencol])

    return fpkm


def calc_fpkm_uq(
    df,
    excol='unstranded',
    lencol='total_exon_length',
    typecol='gene_type',
    chrcol='Chromosome',
):
    '''
    Upper Quartile normalized FPKM

    FPKM-UQ = (C × 1e9)/(UGL)

        C : count of fragments aligned to this gene
        U : upper quartile of fragment counts to autosomal protein-coding genes with count > 0
        G : number of protein-coding genes on autosomes
        L : sum of exon lengths in gene where overlapping exons are merged
    '''
    # selections for U and G
    sel_prot = df[typecol] == 'protein_coding'
    sel_autosomes = ~df[chrcol].isin(['chrX', 'chrY', 'chrM'])
    sel_nonzero = df[excol] > 0
    sel_U = sel_prot & sel_autosomes & sel_nonzero
    sel_G = sel_prot & sel_autosomes

    # Calculate U and G
    U = np.quantile(df[sel_U][excol], 0.75)
    G = len(df[sel_G])

    # Calculate FPKM-UQ
    fpkm_uq = df[excol] * 1e9 / (U * G * df[lencol])
    return fpkm_uq


def save_result(df, outfile, pragma_line=''):
    '''
    Write output table as TSV
    '''

    out_cols = [
        'gene_id',
        'gene_name',
        'gene_type',
        'stranded_first',
        'stranded_second',
        'unstranded',
        'tpm_unstranded',
        'fpkm_unstranded',
        'fpkm_uq_unstranded',
    ]

    with open(outfile, 'w') as out:
        print(pragma_line, file=out)
        df[out_cols].to_csv(
            out, sep='\t', header=True, index=False, float_format='%.4f'
        )


def augment(args, logger):
    '''
    Augment STAR read counts with normalized counts and gene info

    Adds gene names and gene bio-types from data extracted from a GENCODE gtf
    Adds TPM, FPKM, FPKM-UQ normalized counts
    Adds pragma indicating gene annotation version used

    Writes an output file as TSV with a header line and pragma line
    '''

    # load data
    logger.info("Reading counts file {}".format(args.input))
    logger.info("Reading gene_info file {}".format(args.gene_info))
    # with open(args.input) as counts, open(args.gene_info) as gene_info:
    df, edf = prep_data(args.input, args.gene_info)

    # calculate new normalized counts
    logger.info("Calculating normalized counts")
    df['fpkm_unstranded'] = calc_fpkm(df)
    df['fpkm_uq_unstranded'] = calc_fpkm_uq(df)
    df['tpm_unstranded'] = calc_tpm(df)

    # add back extra alignment stats
    df = pd.concat([edf, df], axis=0)

    # write output table
    logger.info("Saving results to {}".format(args.output))
    save_result(df, outfile=args.output, pragma_line=' '.join(args.pragma_line))


def main(args):
    """
    Main entrypoint for augment_counts_table.
    """
    logger = get_logger("augment_counts_table")
    logger.info("Augmenting STAR gene counts file {}.".format(args.input))

    augment(args, logger)
