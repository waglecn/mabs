#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

ref_name = sys.argv[1]
ref_path = sys.argv[2]
try:
    assert os.path.exists(ref_path)
except AssertionError:
    print("reference: {} file not found".format(ref_path), file=sys.stderr)

# ref_path = os.path.join(
#     '..', 'resources', 'alignment_references', f"{ref_name}.fasta"
# )

records = [r for r in SeqIO.parse(ref_path, 'fasta')]
print('>{}\n{}\n'.format(
    "{}_{}".format(ref_name, records[0].id),
    records[0].seq
))
