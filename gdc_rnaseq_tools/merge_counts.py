"""A gdc-rnaseq-tools subcommand to format and merge STAR gene counts
files from the same sample.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
from collections import OrderedDict

from gdc_rnaseq_tools.utils import get_logger, get_open_function

COLUMN_NAMES = ["gene", "unstranded", "stranded_first", "stranded_second"]


def process_files(args, logger):
    """
    All the logic for formatting/merging STAR gene counts.
    :param args: argparser
    :param logger: `logging.Logger` instance
    """
    writer = get_open_function(args.output)
    logger.info("Writing outputs to {0}".format(args.output))

    with writer(args.output, 'wt') as o:
        # Write header row as comment
        o.write("#" + "\t".join(COLUMN_NAMES) + "\n")
        if len(args.input) > 1:
            logger.info("Merging {0} STAR gene counts files.".format(
                        len(args.input)))
            # Load
            dic = OrderedDict()
            for fil in args.input:
                dic = load_star_file(fil, dic)

            logger.info("Writing merged STAR gene counts to {0}.".format(
                        args.output))
            # Merge and write 
            for gene, counts in merge_star_counts(dic):
                row = [gene] + [str(i) for i in counts]
                o.write('\t'.join(row) + '\n')

        else:
            logger.info("Only 1 STAR gene counts file provided. " +
                        "A new STAR gene counts file will be produced " +
                        "with a header line.")
            logger.info("Writing formatted STAR gene counts to {0}.".format(
                        args.output))

            fil = args.input[0]
            reader = get_open_function(fil)
            with reader(fil, 'rt') as fh:
                for line in fh:
                    o.write(line)


def load_star_file(fil, dic):
    """
    Load star counts file into a dictionary. 
    :param fil: path to STAR counts file to load
    :param dic: ``OrderedDict`` to load file to
    :returns: updated ``OrderedDict``
    """
    reader = get_open_function(fil)
    with reader(fil, 'rt') as fh:
        for line in fh:
            cols = line.rstrip('\r\n').split('\t')
            key = cols[0]
            if key not in dic:
                dic[key] = []
            counts = list(map(int, cols[1:]))
            dic[key].append(counts)
    return dic


def merge_star_counts(dic):
    """
    Generator of merged star records from the ordered dic.
    :param dic: ``OrderedDict`` of counts
    :return: the gene key and list of merged counts
    """
    for key in dic:
        counts = [0 for i in range(3)]
        for rec in dic[key]:
            for i,v in enumerate(rec):
                counts[i] += v
        yield key, counts


def main(args):
    """
    Main entrypoint for merge_star_gene_counts.
    """
    logger = get_logger("merge_star_gene_counts")
    logger.info("Merging/Formatting {0} STAR gene counts files.".format(
                len(args.input)))

    process_files(args, logger)
