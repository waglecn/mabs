import os
import pandas as pd
import glob

# configfile: "config.syn.yaml"
if len(config) == 0:
    exit('no configfile specified: use --configfile [file]')


# TODO -- from
"""
from https://snakemake.readthedocs.io/en/v5.22.1/snakefiles/configuration.html#validation
""" # noqa
# validate(config, "config.schema.yaml")
# samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, "samples.schema.yaml")

# TODO -- from
"""
https://charlesreid1.github.io/building-snakemake-command-line-wrappers-for-workflows.html
""" # noqa

exec_dir = os.getcwd()
print(f"executing from: {exec_dir}", file=sys.stderr)

workdir: config['results_dir']
# samples = pandas.read_table(config['sample_sheet'], sep=',', header=None)

# see line above, this is an ugly way to load the sample sheet, pandas
# will definitely provide more validation
samples = [
    line.strip().split(',') for line in open(
        os.path.join(exec_dir, config['sample_sheet']), 'r'
    ) if not
    line.startswith('#')
]
sample_dict = {sample[0]: sample[1:] for sample in samples}
# print(sample_dict)
sample_names = sample_dict.keys()
# sample_names = list(sample_names)[:20]
# print(sample_names, file=sys.stderr)

filt_samples = list(sample_names)
stage2_excluded_samples = config['exclude_samples']
for name in stage2_excluded_samples:
    try:
        filt_samples.remove(name)
    except Exception as e:
        print(
            'to_exclude {} not found in sample_sheet'.format(name),
            file=sys.stderr
        )


def MRCA_mapped_ref_input(sample_name):
    '''
    Use as input function that maps sample_name to reference
        input:
            sample_name: key should be found in sample_dict
        returns:
            ref_name: should match available references
    '''
    from scripts.tree_MRCA import tree_MRCA
    MRCA = tree_MRCA(
        '{}/{}/mashtree/assembly_mashtree.complete.tree'.format(
            exec_dir, config['results_dir']
        ),
        sample_name
    )
    # print(sample_name, MRCA, file=sys.stderr)
    return MRCA


def ref_from_QC(sample_name, QC_summary_file="QC_summary.csv"):
    '''
    Every rule using this function must require QC_summary.csv in input files
    '''
    ref = 'Unknown'
    try:
        QC_data = pd.read_csv(QC_summary_file)
        QC_data.index = QC_data['sample']
        ref = QC_data.loc[sample_name]['MRCA_ref']
    except FileNotFoundError:
        print('QC summary file not found', file=sys.stderr)
    except KeyError:
        print(f'sample {sample_name} not in QC summary', file=sys.stderr)
    return ref


# Rule "all" default catches output of other rules as input in order to
# simplify running the workflow
# - add target rules here as they are implemented

container: "docker://continuumio/miniconda3:4.9.2"
# external data rules + GATK initialization
include: "stage0.smk"
# QC pipeline
include: "stage1.smk"
# variant calling
include: "stage2.smk"

rule all:
    input:


rule stage1:
    input:
        ###############
        # STAGE 1 - QC
        ###############
        "QC_summary.csv"

rule stage1_all_outputs:
    input:
        # pre-trim QC
        expand(
            "pre_trim_QC/{sample}.{R}_fastqc.html",
            sample=sample_names,
            R=["R1", "R2"],
        ),
        # trimmed input
        expand(
            "trimmed_input/{sample}.{read}.fastq.gz",
            sample=sample_names,
            read=["R1", "R2", "S1", "S2"],
        ),
        # post-trim QC
        expand(
            "post_trim_QC/{sample}.{read}_fastqc.html",
            sample=sample_names,
            read=["R1", "R2", "S1", "S2"],
        ),
        # Kraken contamination, paired and single
        expand(
            "trimmed_kraken/{sample}.trimmed.paired",
            sample=sample_names,
        ),
        expand(
            "trimmed_kraken/{sample}.trimmed.single",
            sample=sample_names,
        ),
        # assembly with shovill
        expand(
            "shovill_assembly/{sample}.shovill.contigs.fa",
            sample=sample_names,
        ),
        # kraken2 the assembly
        expand(
            "assembly_kraken/{sample}.assembly",
            sample=sample_names,
        ),
        # MRCA ref per sample
        expand(
            "MRCA/{sample}.csv",
            sample=sample_names,
        ),
        # determine the basic Erm(41) status with TBLASTN
        expand(
            "erm41_status/{sample}.erm41.status",
            sample=sample_names,
        ),
        "mashtree/assembly_mashtree.complete.tree",
        "mashtree/assembly_mashtree.complete.matrix",
        # type assemblies using MRCAs from mashtree
        "QC_summary.csv"

rule stage2:
    input:
        ##################################
        # STAGE 2 - type-specific analyses
        ##################################
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.mpileup".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter.vcf.gz".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter"
            ".failed.vcf.gz".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter"
            ".AD_failed.vcf.gz".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter.hvar"
            ".vcf.gz".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.0cov.bed".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter.hvar_DF.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],

rule stage3:
    input:
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.0cov.bed",
            ref={ref_from_QC(s) for s in filt_samples},
        ),
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.DF.bed",
            ref={ref_from_QC(s) for s in filt_samples},
        ),
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.bed",
            ref={ref_from_QC(s) for s in filt_samples},
        ),
        [
            "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.vcf.gz".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter"
            ".consensus.fasta".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        [
            "gubbins/{ref}/{sample}.RG_SC_RA_bedfilter_gubbins.fasta".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],
        expand(
            "gubbins/{ref}.concatenated.fasta", 
            ref={ref_from_QC(s)  for s in filt_samples}
        ),
        expand(
            "gubbins/{ref}.gubbins.done",
            ref={ref_from_QC(s)  for s in filt_samples}
        ),
        expand(
            "gubbins/{ref}.gubbins.bed",
            ref={ref_from_QC(s) for s in filt_samples}
        ),
        expand(
            "MRCA_ref_mapping/{ref}.merge.vcf.gz",
            ref={ref_from_QC(s) for s in filt_samples}
        ),
        expand(
            "MRCA_ref_mapping/{ref}.merge.vcf.gz.csi",
            ref={ref_from_QC(s) for s in filt_samples}
        ),
        expand(
            "MRCA_ref_mapping/{ref}.merge_gubbins.vcf.gz",
            ref={ref_from_QC(s) for s in filt_samples}
        ),
        expand(
            "SNP_phylo/{ref}.merge.fasta",
            ref={ref_from_QC(s) for s in filt_samples}
        ),
        expand(
            "SNP_phylo/{ref}.merge.fasta.treefile",
            ref={ref_from_QC(s) for s in filt_samples}
        )



######################################################################
#
# Utility RULES
#
######################################################################
rule index_bam:
    threads: 1
    input:
        "{inbam}.bam"
    output:
        "{inbam}.bam.bai"
    shell:
        "samtools index {input} "

rule index_vcf:
    threads: 1
    input:
        "{invcf}.vcf.gz"
    output:
        "{invcf}.vcf.gz.csi"
    shell:
        "htsfile {input} ;"
        "bcftools index {input}"

rule create_fasta_dict_index:
    threads: 1
    params:
        execdir = exec_dir
    conda: "conda_envs/picard.yaml"
    input:
        f"{exec_dir}/resources/alignment_references/{{ref}}.fasta"
    output:
        seqdict = f"{exec_dir}/resources/alignment_references/{{ref}}.dict",
        fai = f"{exec_dir}/resources/alignment_references/{{ref}}.fasta.fai"
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output.seqdict} && "
        "samtools faidx {input} -o {output.fai}"






###############################################################################
# STAGE 3
###############################################################################
# merge all zero-coverage positions for all samples
# def glob_0covbed(ref, samples, step):
rule merge_0cov_bed:
    threads: 1
    input:
        lambda wildcards: [
            "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA.0cov.bed".format(
                ref=ref_from_QC(s), s=s)
            for s in filt_samples if
            ref_from_QC(s) == wildcards.ref
        ]
    output:
        "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.0cov.bed"
    shell:
        "cat {input} | "
        "grep 'NC_010394\.1' -v | sort -k 1,1n -k 2,2n | uniq | "
        "bedtools merge -i - > {output}"

rule merge_DF_bed:
    threads: 1
    input:
        lambda wildcards: [
            "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar_DF.bed".format(
                ref=ref_from_QC(s), s=s)
            for s in filt_samples if
            ref_from_QC(s) == wildcards.ref
        ]
    output:
        "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.DF.bed"
    shell:
        "cat {input} | "
        r"grep 'NC_010394\.1' -v | sort -k 1,1n -k 2,2n | uniq | "
        "bedtools merge -i -> {output}"

rule merge_bed:
    threads: 1
    input:
        zerocov = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.0cov.bed",
        # this is probably covered by gubbins
        # DF = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.DF.bed",
        PEPPE = f"{exec_dir}/resources/alignment_references/{{ref}}.PE_PPE.bed"
    output:
        "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.bed"
    shell:
        "cat {input} | "
        r"grep 'NC_010394\.1' -v | sort -k 1,1n -k 2,2n | cut -f 1-3 | "
        "bedtools merge -i - > {output}"

rule filter_vcf_with_bed:
    threads: 1
    input:
        bed = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.bed",
        vcf = "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter.hvar.vcf.gz",
        csi = "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_filter.hvar.vcf.gz.csi"
    output:
        vcf = "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.vcf.gz"
    shell:
        "bcftools view {input.vcf} -T ^{input.bed} -Oz -o {output.vcf} ;"

rule extract_sample_VCF_SNPs:
    threads: 1
    params:
        execdir = exec_dir
    input:
        vcf = "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.vcf.gz",
        csi = "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.vcf.gz.csi"
    output:
        "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.consensus.fasta",
    shell:
        "bcftools consensus -i 'type=\"SNP\"' "
        "-f {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        " {input.vcf} > {output}; "

rule reheader_gubbins:
    threads: 1
    params:
        execdir = exec_dir
    conda: "conda_envs/phy_plots.yaml"
    input:
        "filtered_vcf/{ref}/{sample}.RG_SC_RA_bedfilter.consensus.fasta"
    output:
        "gubbins/{ref}/{sample}.RG_SC_RA_bedfilter_gubbins.fasta"
    shell:
        "{params.execdir}/scripts/reheader_gubbins.py {input} > {output}"

rule concatenate_gubbins:
    # 10215* 2/1/2021 17:31  cat ./mabscessus/*.fasta > mabs.fasta
    threads: 1
    conda: "conda_envs/phy_plots.yaml"
    params:
        execdir = exec_dir
    input:
        samples = lambda wildcards: [
            "gubbins/{ref}/{s}.RG_SC_RA_bedfilter_gubbins.fasta".format(
                ref=ref_from_QC(s), s=s
            ) for s in filt_samples if ref_from_QC(s) == wildcards.ref
        ],
    output:
        "gubbins/{ref}.concatenated.fasta"
    shell:
        "cat {input} > {output} ; {params.execdir}/scripts/gubbins_ref.py "
        "{wildcards.ref} >> {output}"

rule run_gubbins:
    # 10218* 2/1/2021 17:31  run_gubbins.py --prefix mabs --min_snps 20 --threads 8 mabs.fasta
    conda: "conda_envs/gubbins.yaml"
    threads: 8
    params:
        execdir = exec_dir,
        prefix = "gubbins/{ref}.concatenated"
    input:
        "gubbins/{ref}.concatenated.fasta"
    output:
        # filtered_fasta = "gubbins/{ref}.concatenated."
        # "filtered_polymorphic_sites.fasta",
        # embl = "gubbins/{ref}.concatenated.recombination_predictions.embl",
        sentinel = "gubbins/{ref}.gubbins.done"
    shell:
        "run_gubbins.py --prefix {params.prefix} --min_snps 20 "
        "--threads {threads} {input} || true ; "
        "touch {output.sentinel}"

rule make_gubbins_bed:
    # 10239* 2/1/2021 18:06  ../../analysis/make_bed_from_gubbins.py mabs.recombination_predictions.embl | sort -k 1,1n -k 2,2n > mabs.gubbins.bed
    threads: 1
    conda: "conda_envs/phy_plots.yaml"
    params:
        execdir = exec_dir
    input:
        # sentinel implies this exists
        # embl = "gubbins/{ref}.concatenated.recombination_predictions.embl",
        sentinel = "gubbins/{ref}.gubbins.done"
    output:
        "gubbins/{ref}.gubbins.bed"
    shell:
        "{params.execdir}/scripts/make_bed_from_gubbins.py {input.sentinel} > "
        "{output}"

rule merge_vcf:
    threads: 8
    input:
        vcfs = lambda wildcards: [
            "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.vcf.gz".format(
                ref=ref_from_QC(s), s=s
            ) for s in filt_samples if ref_from_QC(s) == wildcards.ref
        ],
        indices = lambda wildcards: [
            "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.vcf.gz.csi".format(
                ref=ref_from_QC(s), s=s
            ) for s in filt_samples if ref_from_QC(s) == wildcards.ref
        ],
        
    output:
        "MRCA_ref_mapping/{ref}.merge.vcf.gz"
    shell:
        "bcftools merge -0 --threads {threads} {input.vcfs} -Oz -o {output} "

rule filter_merged_vcf:
    threads: 1
    input:
        vcf = "MRCA_ref_mapping/{ref}.merge.vcf.gz",
        merge_bed = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.bed",
        gubbins_bed =  "gubbins/{ref}.gubbins.bed",
    output:
        "MRCA_ref_mapping/{ref}.merge_gubbins.vcf.gz"
    shell:
        "bcftools view {input.vcf} -T ^{input.merge_bed} -Ov | "
        "bcftools view - -T ^{input.gubbins_bed} -Oz -o {output}"

rule make_SNP_alignment:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}.merge_gubbins.vcf.gz"
    output:
        "SNP_phylo/{ref}.merge.fasta",
    shell:
        "{params.execdir}/scripts/make_SNP_alignment.py {input} > {output}"


def snpEff_db(ref):
    if ref == 'mabscessus':
        return "Mycobacterium_abscessus_atcc_19977"
    elif ref == 'mmassiliense':
        return "Mycobacterium_abscessus_subsp_bolletii_ccug_48898_jcm_15300_" \
            "gca_000497265"
    elif ref == 'mbolettii':
        return "Mycobacterium_abscessus_subsp_bolletii_bd"


rule snpEff_annotation:
    conda:
        "conda_envs/snpeff.yaml"
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}.{step}.merge.vcf.gz"
    params:
        db = lambda wildcards: snpEff_db(wildcards.ref)
    output:
        "MRCA_ref_mapping/{ref}.{step}.merge.snpeff.vcf"
    shell:
        "snpEff {params.db} {input} > {output}"


rule MRCA_SNPs_phylogeny:
    conda: "conda_envs/iqtree.yaml"
    threads: 8
    input:
        "SNP_phylo/{ref}.merge.fasta"
    output:
        tree = "SNP_phylo/{ref}.merge.fasta.treefile"
    shell:
        "iqtree -s {input} -nt {threads} -m MFP -B 10000 -redo"
