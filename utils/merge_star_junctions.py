"""Formats and merges STAR junction count files from the same sample"""
import argparse
import gzip
from operator import itemgetter

COLUMN_NAMES = ["chromosome", "intron_start", "intron_end", "strand", 
    "intron_motif", "annotation", "n_unique_map", "n_multi_map",
    "max_splice_overhang"
]

class StarJunctionRecord:
    """Represents a row in the SJ file"""
    def __init__(self, chromosome, intron_first, intron_last,
                 strand, motif, annotation, n_unique_mapped,
                 n_multi_mapped, max_splice_overhang):
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
        return (self.chromosome, self.intron_first, self.intron_last, self.strand, self.motif, self.annotation)

    @classmethod
    def from_line(cls, line):
        """
        Initialize record from a line.
        """
        return cls(*line.rstrip('\r\n').split('\t')) 

    def __iadd__(self, other):
        """
        Used to merge overlapping records by adding the appropriate columns and
        getting the max of the splice overhang column.
        """
        assert self.key == other.key

        self.n_unique_mapped += other.n_unique_mapped 
        self.n_multi_mapped += other.n_multi_mapped 
        self.max_splice_overhang = max(self.max_splice_overhang, other.max_splice_overhang)
        return self

    def __str__(self):
        return "\t".join(map(str, [
            self.chromosome, self.intron_first, self.intron_last, self.strand,
            self.motif, self.annotation, self.n_unique_mapped, self.n_multi_mapped,
            self.max_splice_overhang]))

def merge_files(args):
    """Format and merge the files."""
    writer = get_handler(args.output)
    with writer(args.output, 'wt') as o:
        # Write out header as a comment
        o.write('#' + '\t'.join(COLUMN_NAMES) + '\n')
        if len(args.input) > 1:
            dic = dict() 
            for fil in args.input:
                reader = get_handler(fil) 
                for line in reader(fil, 'rt'):
                   rec = StarJunctionRecord.from_line(line) 
                   if rec.key not in dic:
                       dic[rec.key] = rec
                   else:
                       dic[rec.key] += rec

            for key in sorted(dic, key=itemgetter(0, 1, 2)):
                o.write(str(dic[key]) + '\n') 
        else:
            fil = args.input[0]
            reader = get_handler(fil) 
            for line in reader(fil, 'rt'):
                o.write(line) 

def get_handler(fil):
    """Gets open function based on file name"""
    if fil.endswith('.gz'):
        return gzip.open
    else:
        return open

def load_args():
    """Load arg parser"""
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', action='append', required=True,
        help='Path to junction tab file. Use one or more times')
    p.add_argument('-o', '--output', required=True,
        help='Path to output file') 
    return p.parse_args()

if __name__ == '__main__':
    args = load_args()

    merge_files(args)
