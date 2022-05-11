#!/usr/bin/env python3
"""
This script summarizes the assembly statistics for QC

example usage:
	assembly_summary [assembly.contigs.fasta]

"""

import sys
from Bio import SeqIO
import re

cover = re.compile('cov=(\d+.\d+)')

records = [r for r in SeqIO.parse(sys.argv[1], 'fasta')]
num_contigs = len(records)
len_ass = sum([len(r) for r in records])
cov = []
cov_bp = []
for r in records:
    coverage = float(cover.search(r.description).groups()[0])
    cov.append(coverage)
    cov_bp.append(len(r) * coverage)
sum_cov_bp = sum(cov_bp)
length_weighted_cov = sum_cov_bp / len_ass
print('{},{},{},{:.3f}'.format(
        sys.argv[1], num_contigs, len_ass, length_weighted_cov
    )
)

