"""Main entrypoint for the gdc-rnaseq-tools package.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import argparse

import gdc_rnaseq_tools.augment_star_counts as augment_star_counts
import gdc_rnaseq_tools.merge_counts as merge_star_gene_counts
import gdc_rnaseq_tools.merge_junctions as merge_star_junctions
from gdc_rnaseq_tools.utils import get_logger
from gdc_rnaseq_tools import __version__


def load_args():
    """
    Loads the argument parser object.
    :return: argparse.ArgParser
    """
    parser = argparse.ArgumentParser(
        description="Utility functions for the GDC RNA-Seq workflow"
    )
    parser.add_argument("--version", action="version", version=__version__)
    sp = parser.add_subparsers(description="Select a tool", dest="choice")
    sp.required = True

    # Merge star counts
    gcounts = sp.add_parser(
        "merge_star_gene_counts",
        description="Formats and merges STAR gene " + "counts files.",
    )
    gcounts.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        help="Path to STAR gene counts file. Use one or " + "more times.",
    )
    gcounts.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the merged/formatted output file.",
    )

    # Merge junctions
    jmerge = sp.add_parser(
        "merge_star_junctions",
        description="Formats and merges STAR junction "
        + "count files from the same sample.",
    )
    jmerge.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        help="Path to STAR junction counts file. Use one " + "or more times.",
    )
    jmerge.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the merged/formatted output file.",
    )

    # Augment STAR counts table
    augct = sp.add_parser(
        "augment_star_counts",
        description="Adds FPKM/FPKM-UQ/TPM and gene info columns to STAR"
        + " counts output",
    )
    augct.add_argument(
        "-i",
        "--input",
        required=True,
        default="ReadsPerGene.out.tab",
        help="Path to STAR gene counts file.",
    )
    augct.add_argument(
        "-g",
        "--gene-info",
        required=True,
        default="gene_info.tsv",
        help="Table of gene information with columns: gene_id, "
        + "total_exon_length, gene_name, gene_type, Chromosome",
    )
    augct.add_argument(
        "-o",
        "--output",
        required=False,
        default="counts_report.tsv",
        help="Output file name.",
    )
    augct.add_argument(
        "-v",
        "--gencode-version",
        required=True,
        action='store',
        help="adds a pragma line storing the gencode version to output",
    )

    return parser.parse_args()


def main():
    """Main entry point for CLI"""
    logger = get_logger("gdc-rnaseq-tools")
    args = load_args()

    logger.info("Loading tool {0}".format(args.choice))
    tool = None
    if args.choice == "merge_star_gene_counts":
        tool = merge_star_gene_counts
    elif args.choice == "merge_star_junctions":
        tool = merge_star_junctions
    elif args.choice == "augment_star_counts":
        tool = augment_star_counts

    tool.main(args)
    logger.info("Finished!")


if __name__ == "__main__":
    main()
