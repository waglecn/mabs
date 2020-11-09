#!/usr/bin/env python
"""
This script summarizes the number of reads surviving from a trimmoatic run.
It is actually a duplicate of another script and should be deleted
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
        print('{0}\t{1}\t{2}\t{3}\t{4}'.format(
            pairs, pairs_perc[0], pairs_perc[1],
            int(fwd_perc[0]) + int(rev_perc[0]),
            float(fwd_perc[1]) + float(rev_perc[1])
        ))
        break
