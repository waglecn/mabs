#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

def main(infile=sys.argv[1], name=sys.argv[2]):
    # head, tail = os.path.split(infile)
    # outname = tail.split('.')[0]
    record = [r for r in SeqIO.parse(infile, 'fasta')]
    # this step explicitly dropps the plasmid record[1] in mabs
    print(">{0}\n{1}\n".format(
        name, record[0].seq
    ))

if __name__ == '__main__':
    main()
