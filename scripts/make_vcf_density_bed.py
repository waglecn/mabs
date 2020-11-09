#!/usr/bin/env python3
"""
This script filters a sample's variants based on density. MINWINDOW defines
the minimum window around a SNV where other SNVS are filtered out.

Example usage:
    make_vcf_density_bed.py [input.vcf]

INPUT:
    FILE: input.vcf - a vcf file
OUTPUT:
    FILE: output.bed - a BED format file indicating positions which need to be
                        filtered out

note: as of 20-11-09 this is not working and needs more development

"""
import sys
from progress import bar
import numpy as np

MINWINDOW = 12

records = [r.strip().split('\t') for r in open(sys.argv[1], 'r') if not r.startswith('#')]
chrs = {}
for v in records:
    if v[0] not in chrs:
        chrs[v[0]] = []
    chrs[v[0]].append(v)

for key in chrs:
    length = len(chrs[key])
    pairs = []
    pbar = bar.Bar(max=length)
    current = 0
    x = current + 1
    while current < length - 1:
        if (
            x < length and
            int(chrs[key][current][1]) - int(chrs[key][x][1]) < MINWINDOW
        ):
            pairs.append((current, x))
            x += 1
        else:
            current += 1
            x = current + 1
            pbar.next()
    pbar.finish()

print(len(pairs))
