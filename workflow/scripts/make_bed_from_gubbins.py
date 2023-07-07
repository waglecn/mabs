#!/usr/bin/env python3
import os
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
    assert os.path.exists(sys.argv[3]), exit(
        f"sentinel {sys.argv[3]} not found"
    )

    try:
        lines = open(embl, 'r').readlines()
        lines = [ln.strip() for ln in lines if 'misc_feature' in ln]
        lines = [ln.split(' ')[-1].split('.') for ln in lines]
        print("Found: {}".format(embl), file=sys.stderr)
    except FileNotFoundError:
        print("Not Found: {} -- blank bed".format(embl), file=sys.stderr)
        lines = []

    for ln in lines:
        print('{}\t{}\t{}\tgubbins'.format(
            id, ln[0], ln[-1]
        ))


if __name__ == '__main__':
    main()
