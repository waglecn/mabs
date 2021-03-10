#!/usr/bin/env python3

import sys
import os
import subprocess

infiles = sys.argv[2:]
outfile = sys.argv[1]

if len(infiles) > 1:
    cmd = ['zcat'] + infiles
    with open(outfile, 'w') as outf:
        subprocess.run(cmd, stdout=outf)
    print("Merged eith {}".format(cmd, outfile))
elif len(infiles) == 1:
    os.symlink(infiles[0], outfile)
    print("Created link {} -> {}".format(infiles[0], outfile))
else:
    exit('Error merging input')

