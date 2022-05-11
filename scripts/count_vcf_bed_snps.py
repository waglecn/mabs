#! /usr/bin/env python3

import sys
import vcf
import pandas as pd

in_vcf = sys.argv[1]
bed = sys.argv[2:]

vcf_reader = vcf.Reader(open(sys.argv[1], 'rb'))

snps =[r for r in vcf_reader if r.is_snp]
indels = [r for r in vcf_reader if r.is_indel]
num_snps = len(snps)

num_sites = 0
if len(bed) == 0:
    first_id = next(iter(vcf_reader.contigs))
    num_sites = vcf_reader.contigs[first_id][1]
    # get sites from vcf header for first contig
elif len(bed) == 1:
    data = [
        l.strip().split('\t') for l in open(bed[0], 'r') if not
        l.startswith('#')
    ]
    num_sites = 0
    for l in data:
        num_sites += int(l[2]) - int(l[1])

print(num_snps, len(indels), num_sites)
