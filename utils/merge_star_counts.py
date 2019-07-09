"""Formats and merges STAR gene counts files."""
import argparse
import gzip

from collections import OrderedDict

COLUMN_NAMES = ["gene", "unstranded", "stranded_first", "stranded_second"]

def merge_files(args):
    """
    Make the merge/formatted outputs.
    """
    writer = get_handler(args.output)
    with writer(args.output, 'wt') as o:
        # Write header row as comment
        o.write('#' + '\t'.join(COLUMN_NAMES) + '\n')
        if len(args.input) > 1:
            dic = OrderedDict()
            for fil in args.input:
                reader = get_handler(fil) 
                for line in reader(fil, 'rt'):
                    cols = line.rstrip('\r\n').split('\t') 
                    key = cols[0]
                    if key not in dic:
                        dic[key] = []
                    counts = list(map(int, cols[1:]))
                    dic[key].append(counts) 

            for key in dic:
                counts = [0 for i in range(3)]
                for rec in dic[key]:
                    for i,v in enumerate(rec):
                        counts[i] += v
                row = [key] + [str(i) for i in counts] 
                o.write('\t'.join(row) + '\n') 
        else:
            fil = args.input[0]
            reader = get_handler(fil) 
            for line in reader(fil, 'rt'):
                o.write(line) 

def get_handler(fil):
    """
    Gets open function based on filename.
    """
    if fil.endswith('.gz'):
        return gzip.open
    else:
        return open

def load_args():
    """
    Load argparser.
    """
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', action='append', required=True,
        help='Path to count file. Use one or more times')
    p.add_argument('-o', '--output', required=True,
        help='Path to output file') 
    return p.parse_args()

if __name__ == '__main__':
    args = load_args()

    merge_files(args)
