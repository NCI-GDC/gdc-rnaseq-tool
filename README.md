# gdc-rnaseq-tool
Utility scripts for GDC RNA-seq workflows. The docker file also installs
[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) and
[fqvendorfail](https://github.com/kmhernan/fqvendorfail.git).

## Requirements

* `python>=3.5`

## Merge/Format STAR gene counts

Takes one or more STAR gene counts files from the same sample, adds a header row,
and merges (only if more than 1 is provided) counts by summing across files.

```
usage: gdc-rnaseq-tools merge_star_gene_counts [-h] -i INPUT -o OUTPUT

Formats and merges STAR gene counts files.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to STAR gene counts file. Use one or more times.
  -o OUTPUT, --output OUTPUT
                        Path to the merged/formatted output file.
```

## Merge/Format STAR junctions

Takes one or more STAR junctions files from the same sample, adds a header row,
and merges (only if more than 1 is provided) counts by summing read counts and
taking the max overhang across files.

```
usage: gdc-rnaseq-tools merge_star_junctions [-h] -i INPUT -o OUTPUT

Formats and merges STAR junction count files from the same sample.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to STAR junction counts file. Use one or more
                        times.
  -o OUTPUT, --output OUTPUT
                        Path to the merged/formatted output file.
```
