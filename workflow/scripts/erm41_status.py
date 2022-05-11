#!/usr/bin/env python3
"""
This scriopt summarizes the erm41 blast results for QC

example usage:
	erm41_status.py [blast.out]
"""

import sys

fields = open(sys.argv[1], 'r').readlines()[0].strip().split('\t')

if int(fields[2]) < 100:
    fields.append('Y')
else:
    fields.append('N')
print(','.join(fields))
