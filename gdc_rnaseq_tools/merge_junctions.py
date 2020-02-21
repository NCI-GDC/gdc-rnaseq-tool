"""A gdc-rnaseq-tools subcommand to format and merge STAR junction counts
files from the same sample.

@author: Kyle Hernandez <kmhernan@uchicago.edu>
"""
from operator import itemgetter

from gdc_rnaseq_tools.utils import get_logger, get_open_function

COLUMN_NAMES = [
    "chromosome",
    "intron_start",
    "intron_end",
    "strand",
    "intron_motif",
    "annotation",
    "n_unique_map",
    "n_multi_map",
    "max_splice_overhang",
]


class StarJunctionRecord:
    """Represents a row in the SJ file"""

    def __init__(
        self,
        chromosome,
        intron_first,
        intron_last,
        strand,
        motif,
        annotation,
        n_unique_mapped,
        n_multi_mapped,
        max_splice_overhang,
    ):
        self.chromosome = chromosome
        self.intron_first = int(intron_first)
        self.intron_last = int(intron_last)
        self.strand = int(strand)
        self.motif = int(motif)
        self.annotation = int(annotation)
        self.n_unique_mapped = int(n_unique_mapped)
        self.n_multi_mapped = int(n_multi_mapped)
        self.max_splice_overhang = int(max_splice_overhang)

    @property
    def key(self):
        """
        The first six columns are the identifiers and are used to match
        between input files.
        """
        return (
            self.chromosome,
            self.intron_first,
            self.intron_last,
            self.strand,
            self.motif,
            self.annotation,
        )

    @classmethod
    def from_line(cls, line):
        """
        Initialize record from a line.
        """
        return cls(*line.rstrip("\r\n").split("\t"))

    def __iadd__(self, other):
        """
        Used to merge overlapping records by adding the appropriate columns and
        getting the max of the splice overhang column.
        """
        assert self.key == other.key

        self.n_unique_mapped += other.n_unique_mapped
        self.n_multi_mapped += other.n_multi_mapped
        self.max_splice_overhang = max(
            self.max_splice_overhang, other.max_splice_overhang
        )
        return self

    def __str__(self):
        return "\t".join(
            map(
                str,
                [
                    self.chromosome,
                    self.intron_first,
                    self.intron_last,
                    self.strand,
                    self.motif,
                    self.annotation,
                    self.n_unique_mapped,
                    self.n_multi_mapped,
                    self.max_splice_overhang,
                ],
            )
        )


def process_files(args, logger):
    """
    All the logic for formatting/merging STAR junction counts.
    :param args: argparser
    :param logger: `logging.Logger` instance
    """
    writer = get_open_function(args.output)
    logger.info("Writing outputs to {0}".format(args.output))

    with writer(args.output, "wt") as o:
        # Write header row as comment
        o.write("#" + "\t".join(COLUMN_NAMES) + "\n")
        if len(args.input) > 1:
            logger.info("Merging {0} STAR gene counts files.".format(len(args.input)))
            # Load
            dic = dict()
            for fil in args.input:
                dic = load_junction_file(fil, dic)

            logger.info(
                "Writing merged STAR junction counts to {0}.".format(args.output)
            )
            # Merge and write
            for key in sorted(dic, key=itemgetter(0, 1, 2)):
                o.write(str(dic[key]) + "\n")

        else:
            logger.info(
                "Only 1 STAR junction counts file provided. "
                + "A new STAR junction counts file will be produced "
                + "with a header line."
            )
            logger.info(
                "Writing formatted STAR junction "
                + "counts to {0}.".format(args.output)
            )

            fil = args.input[0]
            reader = get_open_function(fil)
            with reader(fil, "rt") as fh:
                for line in fh:
                    o.write(line)


def load_junction_file(fil, dic):
    """
    Load star junction file into a dictionary.
    :param fil: path to STAR counts file to load
    :param dic: dict to load file to
    :returns: updated dictionary
    """
    reader = get_open_function(fil)
    with reader(fil, "rt") as fh:
        for line in fh:
            rec = StarJunctionRecord.from_line(line)
            if rec.key not in dic:
                dic[rec.key] = rec
            else:
                dic[rec.key] += rec
    return dic


def main(args):
    """
    Main entrypoint for merge_star_gene_counts.
    """
    logger = get_logger("merge_star_junction_counts")
    logger.info(
        "Merging/Formatting {0} STAR junction counts files.".format(len(args.input))
    )

    process_files(args, logger)
