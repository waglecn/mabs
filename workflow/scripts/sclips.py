#!/usr/bin/env python3
'''
This script calculates softclp bases from an input sam/bam file, removes the
read from the alignment if the fraction of solfclipped bases > 0.20

Example usage: sclips.py STDIN

INPUT:
    STDIN: a bam stream from samtools view
OUTPUT:
    STDOUT: a summary of softclipped bases
'''

import sys
import pysam

infile = pysam.AlignmentFile(sys.argv[2], 'rb')
outfile = None

FILTER = False
INVFILTER  = False
REPORT = False
PLOT = False

outfile = pysam.AlignmentFile("-", "wb", template=infile)

if 'filter' in sys.argv:
    FILTER = True
elif 'invfilter' in sys.argv:
    INVFILTER = True
if 'report' in sys.argv:
    REPORT = True
if 'plot' in sys.argv:
    PLOT = True

FILTER_FRAC = 0.20

soft_clipped_bases = 0
hard_clipped_bases = 0

read_softclipped = []
read_fracsoft = []

read_hardclipped = []
read_frachard = []

for r in infile:
    if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
        sclips = sum([c[1] for c in r.cigartuples if c[0] == 4])
        hclips = sum([c[1] for c in r.cigartuples if c[0] == 5])
        read_softclipped.append(sclips)
        read_hardclipped.append(hclips)
        read_length = r.infer_read_length()
        sfrac = sclips / read_length
        hfrac = hclips / read_length
        read_fracsoft.append(sfrac)
        read_frachard.append(hfrac)
        if FILTER and (sfrac < FILTER_FRAC) and (hfrac < FILTER_FRAC):
            outfile.write(r)
        if INVFILTER and (sfrac > FILTER_FRAC) or (hfrac > FILTER_FRAC):
            outfile.write(r)

if REPORT:
    print("{} softcliped bases\t{} hardclipped bases".format(
        sum(read_softclipped), sum(read_hardclipped)
    ), file=sys.stderr)

# Do other stuff, plots
