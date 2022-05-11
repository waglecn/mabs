#!/usr/bin/env python3
"""
This script summarizes the kraken outout for QC

example usage:
    kraken_summary.py [kraken.report]

"""

import sys

infile = sys.argv[1]
rank = 'S'
if sys.argv[2]:
    rank = sys.argv[2]

report = open(infile, 'r')
report = [l.strip().split('\t') for l in report]
report = [
    [
        float(l[0]), int(l[1]), int(l[2]), l[3], int(l[4]), l[5].strip()
    ] for l in report if l[3] == 'S'
]
sorted_report = sorted(report, key=lambda i: i[1])
print('{},{},{},{},{},'.format(
    sorted_report[-1][5],
    sorted_report[-1][0],
    sorted_report[-1][1],
    sorted_report[-1][2],
    sorted_report[-1][3],
), end='')
print('{},{},{},{},{}'.format(
    sorted_report[-2][5],
    sorted_report[-2][0],
    sorted_report[-2][1],
    sorted_report[-2][2],
    sorted_report[-2][3],
))

