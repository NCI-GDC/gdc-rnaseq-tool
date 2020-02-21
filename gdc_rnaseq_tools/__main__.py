"""Main entrypoint for the gdc-rnaseq-tools package.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
import argparse

import gdc_rnaseq_tools.merge_counts as merge_star_gene_counts
import gdc_rnaseq_tools.merge_junctions as merge_star_junctions

from gdc_rnaseq_tools.utils import get_logger


def load_args():
    """
    Loads the argument parser object.
    :return: argparse.ArgParser
    """
    parser = argparse.ArgumentParser(
        description="Utility functions for the GDC RNA-Seq workflow"
    )
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

    tool.main(args)
    logger.info("Finished!")


if __name__ == "__main__":
    main()
