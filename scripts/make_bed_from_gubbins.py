#!/usr/bin/env python3

import os
import pathlib
import sys
from Bio import SeqIO


def main():
    # input gubbins recombination embl
    sentinel = sys.argv[1]
    ref_file = os.path.split(sentinel)[-1].split('.')[0]
    
    if len(sys.argv) > 2:
        id = sys.argv[2]
    else:
        id = 'NC_010397.1'
        base_file = pathlib.Path(sentinel).resolve().parent.parent.parent / 'resources' / \
            'alignment_references' / f"{ref_file}.fasta"
        # base_file = os.path.join(
        #    '..',  
        # )
        records = [r for r in SeqIO.parse(str(base_file), 'fasta')]
        id = records[0].id
        print(
            'Found: {} -- genome accession to use for bed: {}'.format(
                base_file, id
            ), file=sys.stderr
        )


    infile = os.path.join(
        os.path.split(sentinel)[0],
        ref_file + '.concatenated.recombination_predictions.embl'
    )
    

    try:
        lines = open(infile, 'r').readlines()
        lines = [l.strip() for l in lines if 'misc_feature' in l]
        lines =[l.split(' ')[-1].split('.') for l in lines]
        print("Found: {}".format(infile), file=sys.stderr)
    except FileNotFoundError:
        print("Not Found: {} -- blank bed".format(infile), file=sys.stderr)
        lines = []
    
    for l in lines:
        print('{}\t{}\t{}\tgubbins'.format(
            id, l[0], l[-1]
        ))

if __name__ == '__main__':
    main()
