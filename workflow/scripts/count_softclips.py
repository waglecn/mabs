#!/usr/bin/env python3
"""
This script uses pysam to count the number of softclip bases in a sam/bam
for QC

example usage:
	count_softclips.py [input.sam/bam]
"""

import sys
import pysam

infile = pysam.AlignmentFile(sys.argv[1],'rb')

soft_clipped_bases = 0
hard_clipped_bases = 0
for r in infile:
    if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
        sclips = sum([c[1] for c in r.cigartuples if c[0] == 4])
        hclips = sum([c[1] for c in r.cigartuples if c[0] == 5])
        soft_clipped_bases += sclips
        hard_clipped_bases += hclips
print(soft_clipped_bases)
