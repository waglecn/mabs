#!/usr/bin/env python3
"""
This script extracts the allele information from a merged VCF file and produces
a SNV alignment, one sequence for each sample in the VCF

INPUT:
    FILE: merged.vcf
OUTPUT:
    FILE: alignment.fasta

-iterate through the file to find out where the header ends, as
the last line in the header tells us what samples are merged in the VCF and
their order
-based on the samples, we know how many sequences we need to generate
-for each variant(position) we have to add either the REF or ALT allele for
that sample

note: as of 20-11-09 investigating here for issues making the SNP alignment
"""
import sys


def get_samples(header):
    return header[-1].split('\t')[9:]


def is_SNP(REF, ALT):
    if len(REF) > 1:
        return False

    for a in ALT.split(','):
        if len(a) > 1:
            return False
    return True


def get_allele(ALT, FORMAT):
    ALT = ALT.split(',')
    return ALT[int(FORMAT[0]) - 1]


def main():
    header = []
    body = []
    for line in sys.stdin:
        if line.startswith('#'):
            header.append(line.strip())
        else:
            body.append(line.strip().split('\t'))

    align = {}
    samples = get_samples(header)
    align['ref'] = []
    for s in samples:
        align[s] = []

    for line in body:
        assert len(line[9:]) == len(samples)
        # if not is_SNP(line[3], line[4]):
        #     continue
        if line[7].startswith('INDEL'):
            continue
        align['ref'].append(line[3])
        for i, sample in enumerate(line[9:]):
            if sample.startswith('.'):
                align[samples[i]].append(line[3])
            else:
                align[samples[i]].append(get_allele(line[4], sample))

    for s in align:
        print('>{}'.format(s))
        print(''.join(align[s]))



if __name__ == "__main__":
    main()

