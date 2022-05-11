#!/usr/bin/env python3
"""
This script orders the sample names in the project lexicographically, first
numerically by patient then by sample in the series

Example usage:
    order_samples.py [samples.csv]

INPUT:
    FILE: samples.csv - a csv file where the first column is sample name
OUTPUT:
    STDOUT: the lines in sample.csv ordered on sample name

"""

import sys

def order(infile):
    lines = open(infile, 'r')

    samples = [l.strip().split(',')[0] for l in lines if len(l.strip()) > 0]
    samps = []
    for i, s in enumerate(samples):
        s = s.split('-')
        new_s = (int(s[0]), s[1])
        samps.append((new_s, i))

    return sorted(samps, key=lambda s: s[0][0])


def main(filename):
    samples = order(filename)
    lines = [l.strip() for l in open(filename, 'r')]
    for s in samples:
        print(lines[s[1]])



if __name__ == "__main__":
    main(sys.argv[1])