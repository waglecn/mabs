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
"""
import sys
import vcf

MINWINDOW = 12

vcf_reader = vcf.Reader(open(sys.argv[1], 'rb'))
chrom = []

snps = []

current = []
for record in [r for r in vcf_reader if r.is_snp]:
    if record.CHROM not in chrom:
        chrom.append(record.CHROM)
    if current == []:
        current.append(record)
    elif abs(current[-1].POS - record.POS) <= MINWINDOW:
        current.append(record)
        # print(
        #    'added: {} -> {}'.format(current[-2], current[-1]),
        #    file=sys.stderr
        # )
    else:
        snps.append(current)
        # print('new window @ {} bp'.format(record.POS, file=sys.stderr)
        current = []
        current.append(record)
snps.append(current)

if len(chrom) > 1:
    print(
        'WARNING, more than one CHROM present: {}'.format(chrom),
        file=sys.stderr
    )


for window in snps:
    if len(window) > 1:
        positions = [i.POS for i in window]
        start = min(positions)
        end = max(positions)
        print(
            "{}\t{}\t{}\t{} SNPS".format(
                # BED is 0-based, vcf is 1 based, BED is non-inclusive
                chrom[0], start - 1, end, len(window)
            )
        )
