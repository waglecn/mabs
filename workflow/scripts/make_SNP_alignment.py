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
import vcf

try:
    from progress.bar import Bar
    from progress.spinner import Spinner
    from progress.counter import Counter
except ModuleNotFoundError:
    pass


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
    # vcf_reader = [r for r in vcf.Reader(open(sys.argv[1], 'rb'))]
    vcf_reader = vcf.Reader(open(sys.argv[1], 'rb'))
    for c in vcf_reader.contigs:
        default_chrom = c
        break

    align = {s: [] for s in (vcf_reader.samples + ['ref'])}
    try:
        pbar = Counter(suffix='%(index)d%%')
    except Exception:
        pass

    for r in vcf_reader:
        try:
            pbar.next()
        except Exception:
            pass

        
        if r.CHROM == default_chrom and r.is_snp:
            ref_allele = r.REF
            # if 'ref' not in align:
            #     align['ref'] = []
            align['ref'].append(ref_allele)

            for s in r.samples:
                # if s.sample not in align:
                #     align[s.sample] = []
                if s.gt_bases is not None:
                    align[s.sample].append(s.gt_bases[0])
                else:
                    print('{}\t{}\t{}\t{}\t{}'.format(
                        s.sample, r.POS, ref_allele, s.gt_bases, s
                    ), file=sys.stderr)
                    align[s.sample].append('N')
    try:
        pbar.finish()
    except Exception:
        pass
    for s in align:
        print('>{}'.format(s))
        print(''.join(align[s]))


if __name__ == "__main__":
    main()

