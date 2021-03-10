#!/usr/bin/env python
"""
This script counts the number of surviving reads from Trimmomatic output and
formats this data for QC

Example usage:
    trim_csv.py [trimmomatic.log]

INPUT:
    FILE: trimmomatic.log - log file from a trimmomatic run
OUTPUT:
    STDOUT: tsv of data for QC
"""

import re
import sys
import pandas as pd
from Bio import SeqIO
import subprocess
import gzip


def main():

    sample = sys.argv[1]
    # trim QC
    trim_log = sys.argv[2]
    kraken_report = sys.argv[3]
    assembly_fasta = sys.argv[4]
    erm41_file = sys.argv[5]
    depth_bam = sys.argv[6]
    depth_txt = sys.argv[7]
    mlst_txt = sys.argv[8]
    mrca_txt = sys.argv[9]

    trim_qc = trim_QC(sample, trim_log)
    kraken = kraken_QC(sample, kraken_report)
    assembly = assembly_QC(sample, assembly_fasta)
    erm41_QC = erm41_status(sample, erm41_file)
    map_QC = depth_QC(sample, depth_bam, depth_txt)
    mlst_QC = MLST_txt(sample, mlst_txt)
    mrca = MRCA(sample, mrca_txt)

    # print(trim_qc, kraken, assembly, erm41_QC, map_QC, mlst_QC)
    merged = pd.concat(
        [trim_qc, kraken, assembly, erm41_QC, map_QC, mlst_QC, mrca], axis=1
    )
    # print(merged)
    merged.to_csv(sys.stdout)



def trim_QC(sample, trim_log_file):

    in_pairs = re.compile(r'Input Read Pairs: (\d+) ')
    both_surv_perc = re.compile(r' Both Surviving: (\d+) \((\d+.\d+)%\)')
    fwd_only = re.compile(r' Forward Only Surviving: (\d+) \((\d+.\d+)%\)')
    rev_only = re.compile(r' Reverse Only Surviving: (\d+) \((\d+.\d+)%\)')

    for l in open(trim_log_file, 'r'):
        if l.startswith('Input Read Pairs'):
            pairs = in_pairs.match(l).groups()[0]
            pairs_perc = both_surv_perc.search(l).groups()
            fwd_perc = fwd_only.search(l).groups()
            rev_perc = rev_only.search(l).groups()

            df = pd.DataFrame({
                'sample': sample,
                'read_pairs_in': int(pairs),
                'read_pairs_out': int(pairs_perc[0]),
                'read_pairs_frac': float(pairs_perc[1]),
                'single': int(fwd_perc[0]) + int(rev_perc[0]),
                'single_perc': float(fwd_perc[1]) + float(rev_perc[1])
            }, index=[sample])

            break
    return df


def kraken_QC(sample, kraken_report):
    # kraken QC
    rank = 'S'

    report = open(kraken_report, 'r')
    report = [line.strip().split('\t') for line in report]
    report = [
        [
            float(line[0]), int(line[1]), int(line[2]), line[3],
            int(line[4]), line[5].strip()
        ] for line in report if line[3] == rank
    ]
    sorted_report = sorted(report, key=lambda i: i[1])
    df = pd.DataFrame({
        'top_speces': sorted_report[-1][5],
        'top_fraction': sorted_report[-1][0],
        'top_reads': sorted_report[-1][1],
        'top_reads_at_rank': sorted_report[-1][2],
        'top_rank': sorted_report[-1][3],
        'second_species': sorted_report[-2][5],
        'second_fraction': sorted_report[-2][0],
        'second_reads': sorted_report[-2][1],
        'second_reads_at_rank': sorted_report[-2][2],
        'second_rank': sorted_report[-2][3],

    }, index=[sample])
    return df


def assembly_QC(sample, assembly_fasta):
    cover = re.compile(r'cov=(\d+.\d+)')

    records = [r for r in SeqIO.parse(assembly_fasta, 'fasta')]
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

    df = pd.DataFrame({
        'contigs_file': assembly_fasta,
        'num_contigs': num_contigs,
        'assembly_length': len_ass,
        'length_weigthged_cov': length_weighted_cov
    }, index=[sample])
    return df


def erm41_status(sample, erm41_status):
    fields = open(erm41_status, 'r').readlines()[0].strip().split('\t')

    if int(fields[2]) < 100:
        fields.append('Y')
    else:
        fields.append('N')
    df = pd.DataFrame({
        'subject_id': fields[0],
        'id': fields[1],
        'hsp_length': fields[2],
        'mismatches': fields[3],
        'query_start': fields[4],
        'query_end': fields[5],
        'subject_start': fields[6],
        'subject_end': fields[7],
        'erm41_present': fields[8]
    }, index=[sample])
    return df


def depth_QC(sample, input_bam, input_depth):
    cmd = [
        "samtools", "view", "-c", "-F", "2308", input_bam
    ]
    result = int(subprocess.check_output(cmd).decode('utf-8').strip())

    depths = []
    for line in gzip.open(input_depth, 'rb'):
        line = line.decode('utf-8').strip().split('\t')
        if line[0] == 'NC_010397.1':
            depths.append(int(line[2]))
        # depths.append(int([2]))
    total_sites = len(depths)

    total_bases = sum(depths)
    mean_depth = total_bases / total_sites

    num_zero = len([i for i in depths if i == 0])
    num_gt_zero = len([i for i in depths if i > 0])
    num_gte_40 = len([i for i in depths if i >= 40])

    frac_zero = num_zero / total_sites
    frac_gt_zero = num_gt_zero / total_sites
    frac_gte_40 = num_gte_40 / total_sites

    import pysam

    infile = pysam.AlignmentFile(input_bam,'rb')

    soft_clipped_bases = 0
    hard_clipped_bases = 0
    for r in infile:
        if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
            sclips = sum([c[1] for c in r.cigartuples if c[0] == 4])
            hclips = sum([c[1] for c in r.cigartuples if c[0] == 5])
            soft_clipped_bases += sclips
            hard_clipped_bases += hclips

    df = pd.DataFrame({
        'reads_mapped': result,
        'mean_depth': mean_depth,
        'frac_zero': frac_zero,
        'frac_gt_zero': frac_gt_zero,
        'frac_gte_40': frac_gte_40,
        'soft_clipped_bases': soft_clipped_bases
    }, index=[sample])
    return df


def MLST_txt(sample, mlst_file):
    fields = [line.strip().split('\t') for line in open(mlst_file, 'r')][0]
    df = pd.DataFrame({
        'input_mlst_file': fields[0],
        'scheme': fields[1],
        'ST': fields[2],
        'allele1': fields[3],
        'allele2': fields[4],
        'allele3': fields[5],
        'allele4': fields[6],
        'allele5': fields[7],
        'allele6': fields[8],
        'allele7': fields[9],
    }, index=[sample])
    return df


def MRCA(sample, mrca_txt):
    fields = [line.strip().split(',') for line in open(mrca_txt, 'r')][0]
    df = pd.DataFrame({
        'MRCA_ref': fields[0]
    }, index=[sample])
    return df


if __name__ == '__main__':
    main()
