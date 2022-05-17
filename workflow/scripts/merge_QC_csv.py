#!/usr/bin/env python3
"""
    merge_QC_csv.py

    Usage: merge_QC_csv.py outfile infile1 [infile2...]

    outfile: output file in csv format

    infile: input file in csv file
"""

import sys
import pandas as pd

outfile = sys.argv[1]
infiles = sys.argv[2:]

output = None

for inf in infiles:
    try:
        df = pd.read_csv(inf)
        if output is None:
            output = df
        else:
            output = output.append(df, ignore_index=True)
    except Exception as e:
        print(inf, e)

output.index = output['sample']
output.to_csv(outfile)
