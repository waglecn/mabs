#!/usr/bin/env python3
"""
This script summarises the depth for QC

example usage:
	depth_summary.py [genome.depth]
"""
import sys

depths = []
for line in open(sys.argv[1], 'r'):
    line = line.strip().split('\t')
    if line[0] == 'NC_010397.1':
        depths.append(int(line[2]))
    # depths.append(int([2]))
total_sites = len(depths)

total_bases = sum(depths)
mean_depth = total_bases / total_sites

num_zero = len([i for i in depths if i == 0])
num_gt_zero = len([i for i in depths if i > 0])
num_gte_40 = len([i for i in depths if i >= 40])

frac_zero = num_zero / total_sites
frac_gt_zero = num_gt_zero / total_sites
frac_gte_40 = num_gte_40 / total_sites

print('{0:.2f},{1:.2f},{2:.2f},{3:.2f}'.format(
    mean_depth, frac_zero, frac_gt_zero, frac_gte_40
))
