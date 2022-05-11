#! /usr/bin/env python3

import sys
import vcf

vcf_reader = vcf.Reader(filename=sys.argv[1])
snps =[r for r in vcf_reader if r.is_snp]
print(len(snps))
