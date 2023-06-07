#!/usr/bin/env python3

import sys
from pathlib import Path
import gzip
from Bio import SeqIO
import numpy as np
import pandas as pd
from collections import OrderedDict

"""
This script is to extract summary statistics from the QC-steps for each sample
The important stats:
raw short read
    num pairs of reads
    length
    num bases
post-trim QC
    num pairs of reads in
    num pairs of reads out
    fract out
    num unpaired out
    frac unpaired out
trim kraken
    # paired
    first species frac
    first species num
    first species taxon
    second species frac
    second species num
    second species taxon
    unclass frac
    unclass num
    unclass taxon
    # singles
    first species frac
    first species num
    first species taxon
    second species frac
    second species num
    second species taxon
    unclass frac
    unclass num
    unclass taxon
short assembly
    num contigs
    length
    N50
long assembly
    num contigs
    length
    N50
long polish assembly
    num contigs
    length
    N50
species-specific mapping
    reads mapped
    bases in ref
    bases not-covered
    bases covered >=1
    bases covered >=50
    bases covered >=100
    bases covered >=1000
    frac not-covered
    frac covered >=1
    frac covered >=50
    frac covered >=100
    frac covered >=1000
    median coverage
MRCA
MRCA MLST
    ST
    allele(s)
ERM41
    perc id
    length
    bits

 "workflow/scripts/make_sample_QC_csv.py {wildcards.sample} {input.tree}"
        "{input.trim} {input.single_kraken2} {input.pair_kraken2} "
        "{input.short_assembly} {input_long_assembly} "
        "{input.long_polish_assembly} {input.erm41} {input.mabs_bam} "
        "{input.mabs_depth} {input.mlst} {input.MRCA} > {output}"
"""


def parse_fastqc_data(infile):
    data = [l.strip() for l in open(infile, 'r')]
    sections = []
    start = 1
    for i, l in enumerate(data[0:]):
        if ">>END_MODULE" in l:
            sections.append(data[start:i + 1])
            start = i + 1
    # Basic Statistics
    # Per base sequence quality
    # Sequence Length Distribution
    BS = []
    PSQS = []
    SLD = []
    for s in sections:
        if '>>Basic Statistics' in s[0]:
            BS = s
        if '>>Per sequence quality scores' in s[0]:
            PSQS = s
        if '>>Sequence Length Distribution' in s[0]:
            SLD = s
    assert len(BS) != 0
    assert len(PSQS) != 0
    assert len(SLD) != 0

    # parse PSQS
    num = 0
    SQ = 0
    PSQS = [a.split('\t') for a in PSQS if a[0] != '>']
    for i in PSQS[1:]:
        SQ += float(i[0]) * float(i[1])
        num += float(i[1])
    MQ = SQ / num

    total_bases = 0
    approx = False
    SLD = [a.split('\t') for a in SLD if a[0] != '>']
    for i in SLD[1:]:
        if '-' in i[0]:
            approx = True
            i[0] = i[0].split('-')
            low = float(i[0][0])
            high = float(i[0][1])
            mean = (high - low) / 2.0
            total_bases += (mean * float(i[1]))
        else:
            total_bases += (float(i[0]) * float(i[1]))
    Mlen = (total_bases / num)

    total_bases = round(total_bases, -6)

    assert num != 0
    assert MQ != 0
    assert total_bases != '0'

    return(num, MQ, total_bases, Mlen)


def fastqc_short_raw(sample, stats):
    base = Path(f"{stats['results_dir']}/{sample}")

    stats['raw_QC_num_pairs'] = 'N/A'
    stats['raw_QC_reads1'] = 'N/A'
    stats['raw_QC_reads2'] = 'N/A'
    stats['raw_QC_length'] = 'N/A'
    stats['raw_QC_length1'] = 'N/A'
    stats['raw_QC_length2'] = 'N/A'
    stats['raw_QC_num_bases'] = 'N/A'
    stats['raw_QC_bases1'] = 'N/A'
    stats['raw_QC_bases2'] = 'N/A'

    for i in ['1', '2']:
        inf = base / 'input' / f'R{i}_fastqc.data.txt'
        try:
            num, MQ, total_bases, Mlen = parse_fastqc_data(inf)
        except Exception:
            print('Problem opening file: {}'.format(inf), file=sys.stderr)
            return stats

        stats[f'raw_QC_reads{i}'] = num
        stats[f'raw_QC_length{i}'] = Mlen
        stats[f'raw_QC_bases{i}'] = total_bases

    stats['raw_QC_num_pairs'] = stats['raw_QC_reads1']
    stats['raw_QC_length'] = (
        stats['raw_QC_length1'] + stats['raw_QC_length2']
    ) / 2.0
    stats['raw_QC_num_bases'] = stats['raw_QC_bases1'] + stats['raw_QC_bases2']

    return stats


def short_raw(sample, stats):
    results_dir = stats['results_dir']
    base = Path(f'{results_dir}/{sample}')
    from Bio import SeqIO
    stats['raw_QC_num_pairs'] = 'N/A'
    stats['raw_QC_reads1'] = 'N/A'
    stats['raw_QC_reads2'] = 'N/A'
    stats['raw_QC_length'] = 'N/A'
    stats['raw_QC_length1'] = 'N/A'
    stats['raw_QC_length2'] = 'N/A'
    stats['raw_QC_num_bases'] = 'N/A'
    stats['raw_QC_bases1'] = 'N/A'
    stats['raw_QC_bases2'] = 'N/A'

    reads1 = 0
    length1 = 0
    bases1 = 0

    for i in ['1', '2']:
        inf = base / 'input' / f'R{i}.fastq.gz'
    try:
        with gzip.open(inf, 'rt') as handle:
            records1 = SeqIO.parse(handle, format='fastq')
            for r in records1:
                reads1 += 1
                length1 = len(r)
                bases1 += length1
    except Exception:
        print('Problem opening file: {}'.format(inf), file=sys.stderr)
        return stats

    reads2 = 0
    length2 = 0
    bases2 = 0

    inf = base / 'input' / 'R2.fastq.gz'
    try:
        with gzip.open(inf, 'rt') as handle:
            records2 = SeqIO.parse(handle, format='fastq')
            for r in records2:
                reads2 += 1
                length2 = len(r)
                bases2 += length2
    except Exception:
        print('Problem opening file: {}'.format(inf), file=sys.stderr)
        return stats

    stats['raw_QC_num_pairs'] = reads1
    stats['raw_QC_reads1'] = reads1
    stats['raw_QC_reads2'] = reads2
    stats['raw_QC_length'] = length1
    stats['raw_QC_length1'] = length1
    stats['raw_QC_length2'] = length2
    stats['raw_QC_num_bases'] = bases1 + bases2
    stats['raw_QC_bases1'] = bases1
    stats['raw_QC_bases2'] = bases2

    # pprint.pprint(stats['raw'])
    return stats


def fastqc_post_trim_QC(sample, stats):
    base = Path(f"{stats['results_dir']}/{sample}")

    stats['post_trim_QC_num_pairs'] = 'N/A'
    stats['post_trim_QC_readsR1'] = 'N/A'
    stats['post_trim_QC_readsR2'] = 'N/A'
    stats['post_trim_QC_readsS1'] = 'N/A'
    stats['post_trim_QC_readsS2'] = 'N/A'
    stats['post_trim_QC_mean_lengthR1'] = 'N/A'
    stats['post_trim_QC_mean_lengthR2'] = 'N/A'
    stats['post_trim_QC_mean_lengthS1'] = 'N/A'
    stats['post_trim_QC_mean_lengthS2'] = 'N/A'
    stats['post_trim_QC_num_paired_bases'] = 'N/A'
    stats['post_trim_QC_num_single_bases'] = 'N/A'
    stats['post_trim_QC_basesR1'] = 'N/A'
    stats['post_trim_QC_basesR2'] = 'N/A'
    stats['post_trim_QC_basesS1'] = 'N/A'
    stats['post_trim_QC_basesS2'] = 'N/A'

    for i in ['R1', 'R2', 'S1', 'S2']:
        inf = base / 'input' / f'{i}.trim_fastqc.data.txt'
        try:
            num, MQ, total_bases, Mlen = parse_fastqc_data(inf)
        except Exception:
            print('Problem opening file: {}'.format(inf), file=sys.stderr)
            return stats

        stats[f'post_trim_QC_reads{i}'] = num
        stats[f'post_trim_QC_mean_length{i}'] = Mlen
        stats[f'post_trim_QC_bases{i}'] = total_bases

    stats['post_trim_QC_num_pairs'] = (
        stats['post_trim_QC_readsR1'] + stats['post_trim_QC_readsR2']
    )
    stats['post_trim_QC_num_paired_bases'] = (
        stats['post_trim_QC_basesR1'] + stats['post_trim_QC_basesR2']
    )
    stats['post_trim_QC_num_single_bases'] = (
        stats['post_trim_QC_basesS1'] + stats['post_trim_QC_basesS2']
    )

    return stats


def post_trim_QC(sample, stats):
    results_dir = stats['results_dir']
    stats['post_trim_QC_num_pairs'] = 'N/A'
    stats['post_trim_QC_R1'] = 'N/A'
    stats['post_trim_QC_R2'] = 'N/A'
    stats['post_trim_QC_S1'] = 'N/A'
    stats['post_trim_QC_S2'] = 'N/A'
    stats['post_trim_QC_mean_lengthR1'] = 'N/A'
    stats['post_trim_QC_mean_lengthR2'] = 'N/A'
    stats['post_trim_QC_mean_lengthS1'] = 'N/A'
    stats['post_trim_QC_mean_lengthS2'] = 'N/A'
    stats['post_trim_QC_num_paired_bases'] = 'N/A'
    stats['post_trim_QC_num_single_bases'] = 'N/A'
    stats['post_trim_QC_basesR1'] = 'N/A'
    stats['post_trim_QC_basesR2'] = 'N/A'
    stats['post_trim_QC_basesS1'] = 'N/A'
    stats['post_trim_QC_basesS2'] = 'N/A'
    base = Path(f'./{results_dir}/{sample}')

    for i in ['R1', 'R2', 'S1', 'S2']:
        inf = base / 'input' / f'{i}.trim.fastq.gz'
        try:
            with gzip.open(inf, 'rt') as handle:
                records = [r for r in SeqIO.parse(handle, format='fastq')]
                reads = len(records)
                length = np.mean([len(r) for r in records])
                bases = sum(len(r) for r in records)

                stats[f'post_trim_QC_{i}'] = reads
                stats[f'post_trim_QC_mean_length{i}'] = length
                stats[f'post_trim_QC_bases{i}'] = bases
        except Exception:
            print('Problem opening file: {}'.format(inf), file=sys.stderr)
            return stats

    stats['post_trim_QC_num_pairs'] = stats['post_trim_QC_R1']
    stats['post_trim_QC_num_paired_bases'] = (
        stats['post_trim_QC_R1'] + stats['post_trim_QC_R2']
    )
    stats['post_trim_QC_num_single_bases'] = (
        stats['post_trim_QC_S1'] + stats['post_trim_QC_S2']
    )

    return stats


def fastqc_long_QC(sample, stats):
    results_dir = stats['results_dir']

    stats['long_QC_num_reads'] = 'N/A'
    stats['long_QC_mean_length'] = 'N/A'
    stats['long_QC_mean_qual'] = 'N/A'
    stats['long_QC_num_bases'] = 'N/A'

    try:
        base = Path(f'./{results_dir}/{sample}')
        inf = base / 'input' / 'long_fastqc.data.txt'
        num, MQ, total_bases, Mlen = parse_fastqc_data(inf)
    except Exception as e:
        print(e, 'Problem opening file a: {}'.format(inf), file=sys.stderr)
        return stats

    stats['long_QC_num_reads'] = num
    stats['long_QC_mean_length'] = Mlen
    stats['long_QC_mean_qual'] = MQ
    stats['long_QC_num_bases'] = total_bases

    return stats


def kraken(sample, stats):
    results_dir = stats['results_dir']

    base = Path(f'./{results_dir}/{sample}')
    for i in ['raw.paired', 'trimmed.paired', 'trimmed.single']:
        stats[f'kraken_{i}_s1_frac'] = 'N/A'
        stats[f'kraken_{i}_s1_num'] = 'N/A'
        stats[f'kraken_{i}_s1_tax'] = 'N/A'
        stats[f'kraken_{i}_s2_frac'] = 'N/A'
        stats[f'kraken_{i}_s2_num'] = 'N/A'
        stats[f'kraken_{i}_s2_tax'] = 'N/A'
        stats[f'kraken_{i}_unc_frac'] = 'N/A'
        stats[f'kraken_{i}_unc_num'] = 'N/A'
        stats[f'kraken_{i}_unc_taxon'] = 'N/A'

        inf = base / 'input' / f'kraken.{i}'

        try:
            data = [line for line in open(inf, 'r')]
        except Exception as e:
            print(e, 'Problem opening file: {}'.format(inf), file=sys.stderr)
            return stats

        data = [line.strip().split('\t') for line in data]

        assert data[0][-1] == 'unclassified'
        stats[f'kraken_{i}_unc_frac'] = float(data[0][0])
        stats[f'kraken_{i}_unc_num'] = int(data[0][1])
        stats[f'kraken_{i}_unc_taxon'] = data[0][-1]

        spec = [d for d in data if d[3] == 'S']
        stats[f'kraken_{i}_s1_frac'] = float(spec[0][0])
        stats[f'kraken_{i}_s1_num'] = int(spec[0][2])
        stats[f'kraken_{i}_s1_tax'] = spec[0][-1].strip()

        stats[f'kraken_{i}_s2_frac'] = float(spec[1][0])
        stats[f'kraken_{i}_s2_num'] = int(spec[1][2])
        stats[f'kraken_{i}_s2_tax'] = spec[1][-1].strip()

    return stats


def short_assembly(sample, stats):
    results_dir = stats['results_dir']
    stats['short_assembly_num_contigs'] = 'N/A'
    stats['short_assembly_length'] = 'N/A'
    stats['short_assembly_N50'] = 'N/A'
    stats['short_assembly_L50'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')
    inf = base / 'shovill_assembly' / 'contigs.fa'
    try:
        records = [r for r in SeqIO.parse(open(inf, 'r'), format='fasta')]
    except Exception as e:
        print(e, "Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    stats['short_assembly_num_contigs'] = len(records)
    stats['short_assembly_length'] = sum([len(r) for r in records])

    # N50 calculation relies on the fact that contigs are length sorted in
    # Shovill output
    cuml = 0
    for i, r in enumerate(records):
        cuml += len(r)
        if (cuml / stats['short_assembly_length']) > 0.5:
            break
    stats['short_assembly_N50'] = len(r)
    stats['short_assembly_L50'] = i
    # pprint.pprint(stats['assembly'])

    return stats


def long_assembly(sample, stats):
    results_dir = stats['results_dir']
    stats['long_assembly_num_contigs'] = 'N/A'
    stats['long_assembly_length'] = 'N/A'
    stats['long_assembly_N50'] = 'N/A'
    stats['long_assembly_L50'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')
    inf = base / 'dflye' / 'contigs.fa'
    try:
        records = [r for r in SeqIO.parse(open(inf, 'r'), format='fasta')]
    except Exception as e:
        print(e, "Problem opening file: {}".format(inf), file=sys.stderr)
        return stats
    stats['long_assembly_num_contigs'] = len(records)
    stats['long_assembly_length'] = sum([len(r) for r in records])

    # N50 calculation relies on the fact that contigs are length sorted in
    # Shovill output
    cuml = 0
    for i, r in enumerate(records):
        cuml += len(r)
        if (cuml / stats['long_assembly_length']) > 0.5:
            break
    stats['long_assembly_N50'] = len(r)
    stats['long_assembly_L50'] = i
    # pprint.pprint(stats['assembly'])

    return stats


def long_polished_assembly(sample, stats):
    results_dir = stats['results_dir']
    stats['long_polish_assembly_num_contigs'] = 'N/A'
    stats['long_polish_assembly_length'] = 'N/A'
    stats['long_polish_assembly_N50'] = 'N/A'
    stats['long_polish_assembly_L50'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')
    inf = base / 'dflye_short_polish' / 'contigs.fa'
    try:
        records = [r for r in SeqIO.parse(open(inf, 'r'), format='fasta')]
    except Exception as e:
        print(e, "Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    stats['long_polish_assembly_num_contigs'] = len(records)
    stats['long_polish_assembly_length'] = sum([len(r) for r in records])

    # N50 calculation relies on the fact that contigs are length sorted in
    # Shovill output
    cuml = 0
    for i, r in enumerate(records):
        cuml += len(r)
        if (cuml / stats['long_polish_assembly_length']) > 0.5:
            break
    stats['long_polish_assembly_N50'] = len(r)
    stats['long_polish_assembly_L50'] = i
    # pprint.pprint(stats['assembly'])

    return stats


def mapping(sample, stats):
    results_dir = stats['results_dir']
    stats['map_mapped_reads'] = 'N/A'
    stats['map_ref_bases'] = 'N/A'
    stats['map_bases_0cov'] = 'N/A'
    stats['map_bases_gte1cov'] = 'N/A'
    stats['map_bases_gte50cov'] = 'N/A'
    stats['map_bases_gte100cov'] = 'N/A'
    stats['map_bases_gte1000cov'] = 'N/A'
    stats['map_frac_0cov'] = 'N/A'
    stats['map_frac_gte1cov'] = 'N/A'
    stats['map_frac_gte50cov'] = 'N/A'
    stats['map_frac_gte100cov'] = 'N/A'
    stats['map_frac_gte1000cov'] = 'N/A'

    import pysam

    base = Path(f'./{results_dir}/{sample}')

    inf = base / 'ref_mapping' / 'mabs' / 'merged.sorted.bam'
    try:
        sam = pysam.AlignmentFile(inf, "rb")
        stats['map_mapped_reads'] = sam.count()
        stats['map_ref_bases'] = sum(sam.lengths)
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    inf = base / 'ref_mapping' / 'mabs' / 'merged.sorted.depth.gz'
    try:
        depths = [line.strip().split('\t') for line in gzip.open(inf, 'rt')]
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    cov = [int(item[-1]) for item in depths]

    # median coverage
    stats['map_bases_0cov'] = len([i for i in cov if i < 1])
    stats['map_bases_gte1cov'] = len([i for i in cov if i > 1])
    stats['map_bases_gte50cov'] = len([i for i in cov if i > 50])
    stats['map_bases_gte100cov'] = len([i for i in cov if i > 100])
    stats['map_bases_gte1000cov'] = len([i for i in cov if i > 1000])
    stats['map_frac_0cov'] = \
        stats['map_bases_0cov'] / stats['map_ref_bases']
    stats['map_frac_gte1cov'] = \
        stats['map_bases_gte1cov'] / stats['map_ref_bases']
    stats['map_frac_gte50cov'] = \
        stats['map_bases_gte50cov'] / stats['map_ref_bases']
    stats['map_frac_gte100cov'] = \
        stats['map_bases_gte100cov'] / stats['map_ref_bases']
    stats['map_frac_gte1000cov'] = \
        stats['map_bases_gte1000cov'] / stats['map_ref_bases']
    stats['map_median_cov'] = np.median(cov)

    # pprint.pprint(stats['map'])
    return stats

def long_mapping(sample, stats):
    results_dir = stats['results_dir']
    stats['long_map_mapped_reads'] = 'N/A'
    stats['long_map_ref_bases'] = 'N/A'
    stats['long_map_bases_0cov'] = 'N/A'
    stats['long_map_bases_gte1cov'] = 'N/A'
    stats['long_map_bases_gte50cov'] = 'N/A'
    stats['long_map_bases_gte100cov'] = 'N/A'
    stats['long_map_bases_gte1000cov'] = 'N/A'
    stats['long_map_frac_0cov'] = 'N/A'
    stats['long_map_frac_gte1cov'] = 'N/A'
    stats['long_map_frac_gte50cov'] = 'N/A'
    stats['long_map_frac_gte100cov'] = 'N/A'
    stats['long_map_frac_gte1000cov'] = 'N/A'

    import pysam

    base = Path(f'./{results_dir}/{sample}')

    inf = base / 'ref_mapping' / 'mabs' / 'longmerged.sorted.bam'
    try:
        sam = pysam.AlignmentFile(inf, "rb")
        stats['long_map_mapped_reads'] = sam.count()
        stats['long_map_ref_bases'] = sum(sam.lengths)
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    inf = base / 'ref_mapping' / 'mabs' / 'longmerged.sorted.depth.gz'
    try:
        depths = [line.strip().split('\t') for line in gzip.open(inf, 'rt')]
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    cov = [int(item[-1]) for item in depths]

    # median coverage
    stats['long_map_bases_0cov'] = len([i for i in cov if i < 1])
    stats['long_map_bases_gte1cov'] = len([i for i in cov if i > 1])
    stats['long_map_bases_gte50cov'] = len([i for i in cov if i > 50])
    stats['long_map_bases_gte100cov'] = len([i for i in cov if i > 100])
    stats['long_map_bases_gte1000cov'] = len([i for i in cov if i > 1000])
    stats['long_map_frac_0cov'] = \
        stats['long_map_bases_0cov'] / stats['long_map_ref_bases']
    stats['long_map_frac_gte1cov'] = \
        stats['long_map_bases_gte1cov'] / stats['long_map_ref_bases']
    stats['long_map_frac_gte50cov'] = \
        stats['long_map_bases_gte50cov'] / stats['long_map_ref_bases']
    stats['long_map_frac_gte100cov'] = \
        stats['long_map_bases_gte100cov'] / stats['long_map_ref_bases']
    stats['long_map_frac_gte1000cov'] = \
        stats['long_map_bases_gte1000cov'] / stats['long_map_ref_bases']
    stats['long_map_median_cov'] = np.median(cov)

    # pprint.pprint(stats['map'])
    return stats

def MRCA(sample, stats):
    results_dir = stats['results_dir']
    stats['MRCA_ref'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.MRCA.csv'
    try:
        data = [l.strip().split(',') for l in open(inf, 'r')]
        stats['MRCA_ref'] = data[0]
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    return stats


def long_MRCA(sample, stats):
    results_dir = stats['results_dir']
    stats['MRCA_long_ref'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.long.MRCA.csv'
    try:
        data = [l.strip().split(',') for l in open(inf, 'r')]
        stats['MRCA_long_ref'] = data[0]
    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    return stats


def MLST(sample, stats):
    results_dir = stats['results_dir']
    stats['pubMLST_scheme'] = 'N/A'
    stats['pubMLST_ST'] = 'N/A'
    stats['pubMLST_allele1'] = 'N/A'
    stats['pubMLST_allele2'] = 'N/A'
    stats['pubMLST_allele3'] = 'N/A'
    stats['pubMLST_allele4'] = 'N/A'
    stats['pubMLST_allele5'] = 'N/A'
    stats['pubMLST_allele6'] = 'N/A'
    stats['pubMLST_allele7'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.MLST.csv'
    try:
        data = [l.strip().split("\t") for l in open(inf, 'r')]
        stats['pubMLST_scheme'] = data[0][1]
        stats['pubMLST_ST'] = data[0][2]
        stats['pubMLST_allele1'] = data[0][3]
        stats['pubMLST_allele2'] = data[0][4]
        stats['pubMLST_allele3'] = data[0][5]
        stats['pubMLST_allele4'] = data[0][6]
        stats['pubMLST_allele5'] = data[0][7]
        stats['pubMLST_allele6'] = data[0][8]
        stats['pubMLST_allele7'] = data[0][9]

    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    return stats


def MLST_long(sample, stats):
    results_dir = stats['results_dir']
    stats['pubMLST_long_scheme'] = 'N/A'
    stats['pubMLST_long_ST'] = 'N/A'
    stats['pubMLST_long_allele1'] = 'N/A'
    stats['pubMLST_long_allele2'] = 'N/A'
    stats['pubMLST_long_allele3'] = 'N/A'
    stats['pubMLST_long_allele4'] = 'N/A'
    stats['pubMLST_long_allele5'] = 'N/A'
    stats['pubMLST_long_allele6'] = 'N/A'
    stats['pubMLST_long_allele7'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.long.MLST.csv'
    try:
        data = [l.strip().split("\t") for l in open(inf, 'r')]
        stats['pubMLST_long_scheme'] = data[0][1]
        stats['pubMLST_long_ST'] = data[0][2]
        stats['pubMLST_long_allele1'] = data[0][3]
        stats['pubMLST_long_allele2'] = data[0][4]
        stats['pubMLST_long_allele3'] = data[0][5]
        stats['pubMLST_long_allele4'] = data[0][6]
        stats['pubMLST_long_allele5'] = data[0][7]
        stats['pubMLST_long_allele6'] = data[0][8]
        stats['pubMLST_long_allele7'] = data[0][9]

    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats

    return stats


def erm41(sample, stats):
    results_dir = stats['results_dir']
    stats['erm(41)_contig'] = 'N/A'
    stats['erm(41)_perc_id'] = 'N/A'
    stats['erm(41)_hit_length'] = 'N/A'
    stats['erm(41)_gaps'] = 'N/A'
    stats['erm(41)_query_start'] = 'N/A'
    stats['erm(41)_query_end'] = 'N/A'
    stats['erm(41)_subject_start'] = 'N/A'
    stats['erm(41)_subject_end'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.erm41.status'
    try:
        data = [l.strip().split("\t") for l in open(inf, 'r')]
        stats['erm(41)_contig'] = data[0][0]
        stats['erm(41)_perc_id'] = data[0][1]
        stats['erm(41)_hit_length'] = data[0][2]
        stats['erm(41)_gaps'] = data[0][3]
        stats['erm(41)_query_start'] = data[0][4]
        stats['erm(41)_query_end'] = data[0][5]
        stats['erm(41)_subject_start'] = data[0][6]
        stats['erm(41)_subject_end'] = data[0][7]

    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats
    return stats


def erm41_long(sample, stats):
    results_dir = stats['results_dir']
    stats['erm(41)_long_contig'] = 'N/A'
    stats['erm(41)_long_perc_id'] = 'N/A'
    stats['erm(41)_long_hit_length'] = 'N/A'
    stats['erm(41)_long_gaps'] = 'N/A'
    stats['erm(41)_long_query_start'] = 'N/A'
    stats['erm(41)_long_query_end'] = 'N/A'
    stats['erm(41)_long_subject_start'] = 'N/A'
    stats['erm(41)_long_subject_end'] = 'N/A'

    base = Path(f'./{results_dir}/{sample}')

    inf = base / f'{sample}.long.erm41.status'
    try:
        data = [l.strip().split("\t") for l in open(inf, 'r')]
        stats['erm(41)_long_contig'] = data[0][0]
        stats['erm(41)_long_perc_id'] = data[0][1]
        stats['erm(41)_long_hit_length'] = data[0][2]
        stats['erm(41)_long_gaps'] = data[0][3]
        stats['erm(41)_long_query_start'] = data[0][4]
        stats['erm(41)_long_query_end'] = data[0][5]
        stats['erm(41)_long_subject_start'] = data[0][6]
        stats['erm(41)_long_subject_end'] = data[0][7]

    except Exception:
        print("Problem opening file: {}".format(inf), file=sys.stderr)
        return stats
    return stats


def main():
    try:
        sample_name = sys.argv[1]
        results_dir = sys.argv[2]
    except IndexError:
        print("Sample name not supplied", file=sys.stderr)
        exit()

    stats = OrderedDict()
    stats['sample'] = sample_name
    stats['results_dir'] = results_dir

    print('Collecting stats...', file=sys.stderr)
    stats = fastqc_short_raw(sample_name, stats)
    print('  finished raw read QC stats...', file=sys.stderr)
    stats = fastqc_post_trim_QC(sample_name, stats)
    print('  finished trimmed read QC stats...', file=sys.stderr)
    stats = fastqc_long_QC(sample_name, stats)
    print('  finished long read QC stats...', file=sys.stderr)
    stats = kraken(sample_name, stats)
    print('  finished kraken QC stats...', file=sys.stderr)
    stats = short_assembly(sample_name, stats)
    stats = long_assembly(sample_name, stats)
    stats = long_polished_assembly(sample_name, stats)
    print('  finished assembly QC stats...', file=sys.stderr)
    stats = mapping(sample_name, stats)
    print('  finished mapping QC stats...', file=sys.stderr)
    stats = MLST(sample_name, stats)
    stats = MLST_long(sample_name, stats)
    stats = MRCA(sample_name, stats)
    stats = long_MRCA(sample_name, stats)
    stats = erm41(sample_name, stats)
    stats = erm41_long(sample_name, stats)

    df = pd.DataFrame(
        stats, columns=stats.keys(), index=[sample_name],
    )
    print(df.to_csv())


if __name__ == '__main__':
    main()
