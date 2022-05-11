#!/usr/bin/env python3

import os
import pathlib
import sys
from Bio import SeqIO


def main():
    # input gubbins recombination embl
    embl = sys.argv[1]
    ref_file = sys.argv[2]
    records = [r for r in SeqIO.parse(str(ref_file), 'fasta')]
    id = records[0].id
    print(
        'Found: {} -- genome accession to use for bed: {}'.format(
            ref_file, id
        ),
        file=sys.stderr
    )

    try:
        lines = open(embl, 'r').readlines()
        lines = [l.strip() for l in lines if 'misc_feature' in l]
        lines =[l.split(' ')[-1].split('.') for l in lines]
        print("Found: {}".format(embl), file=sys.stderr)
    except FileNotFoundError:
        print("Not Found: {} -- blank bed".format(infile), file=sys.stderr)
        lines = []

    for l in lines:
        print('{}\t{}\t{}\tgubbins'.format(
            id, l[0], l[-1]
        ))


if __name__ == '__main__':
    main()
