#!/usr/bin/env python3

import os
import pathlib
import sys
from Bio import SeqIO


def main():
    # input CFML importation_status txt
    cfml = sys.argv[1]
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
        lines = [l.strip().split('\t') for l in open(cfml, 'r').readlines()]
        print("Found: {}".format(cfml), file=sys.stderr)
    except FileNotFoundError:
        print("Not Found: {} -- blank bed".format(cfml), file=sys.stderr)
        lines = []

    for l in lines:
        print('{}\t{}\t{}\tgubbins'.format(
            id, l[1], l[2]
        ))


if __name__ == '__main__':
    main()
