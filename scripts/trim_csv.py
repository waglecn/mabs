#!/usr/bin/env python
"""
This script counts the number of surviving reads from Trimmomatic output and
formats this data for QC

Example usage:
    trim_csv.py [trimmomatic.log]

INPUT:
    FILE: trimmomatic.log - log file from a trimmomatic run
OUTPUT:
    STDOUT: tsv of data for QC
"""

import re
import sys

in_pairs = re.compile(r'Input Read Pairs: (\d+) ')
both_surv_perc = re.compile(r' Both Surviving: (\d+) \((\d+.\d+)%\)')
fwd_only = re.compile(r' Forward Only Surviving: (\d+) \((\d+.\d+)%\)')
rev_only = re.compile(r' Reverse Only Surviving: (\d+) \((\d+.\d+)%\)')

infile = sys.argv[1]

for l in open(infile, 'r'):
    if l.startswith('Input Read Pairs'):
        pairs = in_pairs.match(l).groups()[0]
        pairs_perc = both_surv_perc.search(l).groups()
        fwd_perc = fwd_only.search(l).groups()
        rev_perc = rev_only.search(l).groups()
        print('{0},{1},{2},{3},{4}'.format(
            pairs, pairs_perc[0], pairs_perc[1],
            int(fwd_perc[0]) + int(rev_perc[0]),
            float(fwd_perc[1]) + float(rev_perc[1])
        ))
        break
