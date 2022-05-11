#!/usr/bin/env python3

import sys
import os
from pathlib import Path
import subprocess

infiles = sys.argv[2:]
outfile = sys.argv[1]
tempout = Path(outfile).stem

if len(infiles) > 1:
    cmd1 = ['zcat'] + infiles
    with open(tempout, 'w') as outf:
        subprocess.run(cmd1, stdout=outf)
    cmd2 = ['gzip', tempout]
    subprocess.run(cmd2)
    os.rename(tempout, outfile)
    print("Merged: {}".format(cmd2, outfile))
elif len(infiles) == 1:
    os.symlink(infiles[0], outfile)
    print("Created link {} -> {}".format(infiles[0], outfile))
else:
    exit('Error merging input')

