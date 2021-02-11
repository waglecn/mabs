import os
import pandas
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
mash_ass = config['mash_ref_taxa']
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
    return tree_MRCA(
        '{}/results/mashtree/assembly_mashtree.complete.tree'.format(exec_dir),
        sample_name
    )


# Rule "all" default catches output of other rules as input in order to
# simplify running the workflow
# - add target rules here as they are implemented

# external data rules + GATK initialization
include: "stage0.smk"
# QC pipeline
include: "stage1.smk"

rule all:


rule stage1:
    input:
        ###############
        # STAGE 1 - QC
        ###############
        # "merged reads"
        expand(
            "merged_input/{sample}.{R}.fastq.gz",
            sample=sample_names,
            R=["R1", "R2"],
        ),
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
        expand(
            "trimmed_input/{sample}.csv",
            sample=sample_names,
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
        expand(
            "trimmed_kraken/{sample}.csv",
            sample=sample_names,
        ),
        # assembly with shovill
        expand(
            "shovill_assembly/{sample}.shovill.contigs.fa",
            sample=sample_names,
        ),
        expand(
            "shovill_assembly/{sample}.shovill.csv",
            sample=sample_names,
        ),
        # kraken2 the assembly
        expand(
            "assembly_kraken/{sample}.assembly",
            sample=sample_names,
        ),
        # determine the basic Erm(41) status with TBLASTN
        expand(
            "erm41_status/{sample}.erm41.status",
            sample=sample_names,
        ),
        expand(
            "erm41_status/{sample}.erm41.csv",
            sample=sample_names,
        ),
        # map all samples to mabs to determine basic coverage
        expand(
            "ref_mapping/mabs/{sample}.merged.sorted.bam",
            sample=sample_names,
        ),
        expand(
            "ref_mapping/mabs/{sample}.merged.sorted.csv",
            sample=sample_names,
        ),
        "mashtree/assembly_mashtree.complete.tree",
        "mashtree/assembly_mashtree.complete.matrix",
        # type assemblies using MRCAs from mashtree
        expand(
            "MRCA_MLST/{sample}.mlst.txt",
            sample=sample_names
        ),
        expand(
            "QC_summary/{sample}.QC.csv",
            sample=sample_names
        ),
        "QC_summary.csv"




rule stage2:
    input:
        ##################################
        # STAGE 2 - type-specific analyses
        ##################################
        # mashtree assemblies with external refseq genomes
        # map samples to sub-species type strains
        # filter softclips and add readgroups
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam".format(
                ref=MRCA_mapped_ref_input(s),
                sample=s
            ) for s in filt_samples
        ],
        # gatk realignment
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_MQ30_BAQ.mpileup".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],

rule with_realignment:
    input:
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1_DF.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.failed.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.AD_failed.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.0cov.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        expand(
            "{exec_dir}/resources/alignment_references/{ref}.PE_PPE.bed",
            ref=['mabscessus', 'mbolettii', 'mmassiliense'],
            exec_dir=exec_dir
        ),
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.merge.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.hvar.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter_hvar"
            ".consensus.fasta".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ],
        [
            "gubbins/{ref}/{sample}.RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1_DF_"
            "BedFilter_hvar.consensus.fasta".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in filt_samples
        ]

rule SNP_alignment:
    input:
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.vcf.gz",
            ref=['mabscessus', 'mbolettii', 'mmassiliense']
        ),
        expand(
            "SNP_phylo/{ref}.RG_SC_RA.merge.fasta",
            ref=['mabscessus', 'mbolettii', 'mmassiliense']
        ),
        expand(
            "SNP_phylo/{ref}.RG_SC_RA.merge.fasta.treefile",
            ref=['mabscessus', 'mbolettii', 'mmassiliense']
        )

# need rule to produce SNP alignment with QC fail samples excluded

rule SNP_alignment_without_mb:
    input:
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.vcf.gz",
            ref=['mabscessus', 'mmassiliense']
        ),
        expand(
            "SNP_phylo/{ref}.RG_SC_RA.merge.fasta",
            ref=['mabscessus', 'mmassiliense']
        ),
        expand(
            "SNP_phylo/{ref}.RG_SC_RA.merge.fasta.treefile",
            ref=['mabscessus', 'mmassiliense']
        )

# rule without_realignment:
#     input:
#         [
#             "MRCA_ref_mapping/{ref}/{sample}"
#             ".RG_SC_mq30_baq50_dp20_sp60_ad1.vcf".format(
#                 ref=MRCA_mapped_ref_input(s), sample=s
#             ) for s in sample_names
#         ],
#         [
#             "MRCA_ref_mapping/{ref}/{sample}"
#             ".RG_SC_MQ30_BAQ50_DP20_SP60_AD1.failed.vcf".format(
#                 ref=MRCA_mapped_ref_input(s), sample=s
#             ) for s in sample_names
#         ],
#         [
#             "MRCA_ref_mapping/{ref}/{sample}"
#             ".RG_SC_MQ30_BAQ50_DP20_SP60_AD1.AD_failed.vcf".format(
#                 ref=MRCA_mapped_ref_input(s), sample=s
#             ) for s in sample_names
#         ],
#         [
#             "MRCA_ref_mapping/{ref}/{sample}.RG_SC.0cov.bed".format(
#                 ref=MRCA_mapped_ref_input(s), sample=s
#             ) for s in sample_names
#         ]



######################################################################
#
# STAGE 2 RULES
#
######################################################################

# # This rule puts everything in context of a rough mashtree to determine MRCA
# # for downstream rules - MLST, specific references
# rule mashtree_assemblies:
#     threads: 8
#     conda:
#         "conda_envs/mashtree.yaml"
#     input:
#         # this is ugly - but there are fewer contraints on input
#         # wildcards in input must match wildcards in output, but since
#         # this output rule has no wildcards, needed to something else
#         # internal = glob.glob("shovill_assembly/*.shovill.contigs.fa"),
#         # internal2 = glob.glob("../bakMA/4-shovill_assemble/*.contigs.fa"),
#         # glob.glob("external/refseq/complete/*.gbff.gz")
#         aggregate_refseq,
#         internal = expand(
#             "shovill_assembly/{sample}.shovill.contigs.fa",
#             sample=sample_names
#         )
#     params:
#         # the mashtree tempdir can get large, depending on the number of
#         # species - see config.yaml
#         tempdir = config['mashtree_tempdir']
#     output:
#         tree = "mashtree/assembly_mashtree.{nodeset}.tree",
#         matrix = "mashtree/assembly_mashtree.{nodeset}.matrix"
#     log:
#         "mashtree/assembly_mashtree.{nodeset}.log"
#     shell:
#         "mashtree --mindepth 0 --numcpus {threads} --outtree {output.tree} "
#         "--sort-order random "
#         "--outmatrix {output.matrix} --tempdir {params.tempdir} "
#         "{input} 2>&1 | tee {log}"

# # This rule determines which samples will use which reference or MLST scheme
# rule MRCA_MLST:
#     threads: 1
#     conda:
#         "conda_envs/ngd_phylo.yaml"
#     params:
#         execdir = exec_dir
#     input:
#         tree = "mashtree/assembly_mashtree.complete.tree",
#         assembly = "shovill_assembly/{sample}.shovill.contigs.fa"
#     output:
#         "MRCA_MLST/{sample}.mlst.txt"
#     shell:
#         # "echo $(scripts/tree_MRCA.py {input.tree} {sample}) ;"
#         "mlst --scheme $({params.execdir}/scripts/tree_MRCA.py {input.tree} "
#         "{wildcards.sample}) "
#         "--threads {threads} {input.assembly} > {output}"

# This is modelled after map_to_mabs, except need to sub in the ref
# If the most appropriate ref is to mabs, this should prevent re-calculation
# of the alignment
# TODO - check the names output by the tree_MRCA.py
rule temp_MRCA_ref_alignment:
    threads: 7
    params:
        execdir = exec_dir
    input:
        tree = "mashtree/assembly_mashtree.complete.tree",
        S1 = "trimmed_input/{sample}.S1.fastq.gz",
        R1 = "trimmed_input/{sample}.R1.fastq.gz",
        R2 = "trimmed_input/{sample}.R2.fastq.gz",
    params:
        ref = lambda wildcards: MRCA_mapped_ref_input(wildcards.sample),
    output:
        paired_temp = temp(
            "MRCA_ref_mapping/{ref}/{sample}.paired.sorted.bam"
        ),
        paired_temp_bai = temp(
            "MRCA_ref_mapping/{ref}/{sample}.paired.sorted.bam.bai"
        ),
        single_temp = temp(
            "MRCA_ref_mapping/{ref}/{sample}.single.sorted.bam"
        ),
        single_temp_bai = temp(
            "MRCA_ref_mapping/{ref}/{sample}.single.sorted.bam.bai"
        ),
        merge_temp = temp("MRCA_ref_mapping/{ref}/{sample}.merge.bam"),
        merge_sorted = temp(
            "MRCA_ref_mapping/{ref}/{sample}.merged.sorted.bam"
        ),
        merge_sorted_bai = temp(
            "MRCA_ref_mapping/{ref}/{sample}.merged.sorted.bam.bai"
        )
    shell:
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/{params.ref} "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/{params.ref} "
        "{input.S1} | samtools view -Sbh - | samtools sort > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort > "
        "{output.merge_sorted} ;"
        "samtools index {output.merge_sorted} {output.merge_sorted_bai}"

rule temp_create_sam_read_groups:
    threads: 1
    conda: "conda_envs/picard.yaml"
    input:
        "MRCA_ref_mapping/{ref}/{sample}.merged.sorted.bam"
    output:
        temp("MRCA_ref_mapping/{ref}/{sample}.RG.merged.sorted.bam")
    log:
        "MRCA_ref_mapping/{ref}/{sample}.picard.RG.log"
    shell:
        "picard AddOrReplaceReadGroups -INPUT {input} "
        "-OUTPUT {output} -SORT_ORDER coordinate -RGID {wildcards.sample} "
        "-RGLB unknown -RGPL Illumina -RGSM {wildcards.sample} "
        "-RGPU project -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT "
        " 2>&1 | tee {log}"

rule MRCA_ref_softclip_filter:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}/{sample}.RG.merged.sorted.bam"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam"
    shell:
        "{params.execdir}/scripts/sclips.py filter {input} > {output}"


rule create_fasta_dict_index:
    threads: 1
    params:
        execdir = exec_dir
    conda: "conda_envs/picard.yaml"
    input:
        "{}/resources/alignment_references/{{ref}}.fasta".format(exec_dir)
    output:
        seqdict = "{params.exec}/resources/alignment_references/{ref}.dict",
        fai = "{params.exec}/resources/alignment_references/{ref}.fasta.fai"
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output.seqdict} && "
        "samtools faidx {input} -o {output.fai}"

rule index_bam:
    threads: 1
    input:
        "{inbam}.bam"
    output:
        "{inbam}.bam.bai"
    shell:
        "samtools index {input} "

rule MRCA_ref_gatk_realignment_intervals:
    threads: 1
    conda:
        "conda_envs/gatk3.yaml"
    params:
        execdir = exec_dir
    input:
        bam = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam",
        bai = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam.bai",
        seqdict = "{}/resources/alignment_references/{{ref}}.dict".format(
            exec_dir
        ),
        fai = "{}/resources/alignment_references/{{ref}}.fasta.fai".format(
            exec_dir
        ),
        gatk = "{}/resources/gatk-registered".format(exec_dir)
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals"
    log:
        "MRCA_ref_mapping/{ref}/{sample}.gatk3_intervals.log"
    shell:
        "gatk3 -T RealignerTargetCreator -R "
        "{params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -o {output} 2>&1 | tee {log}"

rule MRCA_ref_gatk_realignment:
    threads: 1
    params:
        execdir = exec_dir
    conda:
        "conda_envs/gatk3.yaml"
    input:
        bam = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam",
        bai = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam.bai",
        intervals = "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals",
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
    log:
        "MRCA_ref_mapping/{ref}/{sample}.gatk_realign.log"
    shell:
        "gatk3 -rf OverclippedRead "
        # "--filter_is_too_short_value 100 -minRead 100 "
        # "--do_not_require_softclips_both_ends -rf ReadLength -maxRead 500 "
        "-T IndelRealigner "
        "-R {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -targetIntervals {input.intervals} -o {output} "
        " 2>&1 | tee {log}"

rule MRCA_make_mpileup:
    threads: 1
    params:
        execdir = exec_dir
    input:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
        "MRCA_ref_mapping/{ref}/{sample}.{step}.bam"
    output:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_mq30_baq.mpileup"
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ.mpileup"
    shell:
        # note that mpileup has moved to bcftools
        # -d max depth
        # -q min mapQ
        # -t in samtools 1.7 -> -a output tags
        # -u in samtools 1.7 -> -O v uncompressed VCF output
        # -g in samtools 1.7 -> bcftools includes genotype likelihoods default
        "bcftools mpileup -d 1000 -q 30 -a DP,AD,ADF,ADR,SP -Ov "
        "-f {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "{input} -o {output}"

rule MRCA_call_and_filter_variants:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1.vcf"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call -Ov -m -v {input} | "
        # filtering
        # MQ mapping quality
        # ID-SP - Phred strand bias p-value
        # ID-ADF[1] - allelic depth forward strand of first alt allele
        # ID-ADR[1] - alleleic depth reverse strand of first alt allele
        # QUAL - Phred scaled ALT quality
        # FORMAT-DP - number of high-quality bases
        # FORMAT-SP - Phred strand bias p-value
        # FORMAT-ADF and ADR as in INFO
        # Note need to specify sample 0:
        # see https://github.com/samtools/bcftools/issues/757
        "bcftools filter -i 'SP<45 & ADF[0:1]>1 & ADR[0:1]>1 & MQ>30 & "
        "QUAL>50 & FORMAT/DP > 10 & SP<45 & ADF[0:1]>1 & ADR[0:1]>1' "
        "-o {output}"

rule MRCA_inverse_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_"
        "DP20_SP60_AD1.failed.vcf"
    shell:
        "bcftools call -O v -m -v {input} | "
        "bcftools filter -i 'SP>=45 || MQ<=30 || FORMAT/DP<=10 || QUAL<=50' "
        "-o {output}"

rule MRCA_inverse_AD_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_"
        "DP20_SP60_AD1.AD_failed.vcf"
    params:
        snp_cutoff = 0.90
    shell:
        "bcftools call -O v -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff}) & (ADF[0:1]<=1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule density_filter_bed:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1.vcf"
    output:
        bed = "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF.bed"
    shell:
        "{params.execdir}/scripts/make_vcf_density_bed.py {input} > "
        "{output.bed}"

rule compress_vcf:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1.vcf"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1.vcf.gz"
    shell:
        "bcftools view -Oz -o {output} {input} ;"
        "htsfile {output} ;"
        "bcftools index {output}"

rule density_filter_vcf_with_bed:
    threads: 1
    input:
        bed = "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1"
        "_DF.bed",
        vcf = "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1"
        ".vcf.gz",
    output:
        vcf = "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf",
    shell:
        "bcftools view {input.vcf} -T ^{input.bed} -Ov -o {output.vcf} "

# sample BED for zero coverage
rule make_bed_0cov:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.bam"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

# make ref BED for PE  PPE genes
rule make_PE_PPE_BED:
    threads: 1
    params:
        execdir = exec_dir
    conda:
        "conda_envs/biopython.yaml"
    input:
        "{}/resources/alignment_references/{{ref}}.gbk".format(exec_dir)
    output:
        "{params.execdir}/resources/alignment_references/{ref}.PE_PPE.bed"
    shell:
        "{params.execdir}/scripts/make_PE_PPE_BED.py {input} > {output}"

# merge all zero-coverage positions for all samples
rule merge_0cov_bed:
    threads: 1
    input:
        lambda wildcards: glob.glob(
            "MRCA_ref_mapping/{ref}/*.{step}.0cov.bed".format(
                ref=wildcards.ref, step=wildcards.step
            )
        )
    output:
        "MRCA_ref_mapping/{ref}.{step}_merge.0cov.bed"
    shell:
        "cat {input} | "
        "grep 'NC_010394\.1' -v | sort -k 1,1n -k 2,2n | uniq | "
        "bedtools merge -i - > {output}"


rule merge_bed:
    threads: 1
    input:
        bed0 = "MRCA_ref_mapping/{ref}.{step}_merge.0cov.bed",
        PEPPEbed = "{}/resources/alignment_references/" \
            "{{ref}}.PE_PPE.bed".format(
                exec_dir
            ),
    output:
        "MRCA_ref_mapping/{ref}.{step}.merge.bed"
    shell:
        "cat {input.bed0} {input.PEPPEbed} | "
        # this is bad, fix this
        "grep 'NC_010394\.1' -v | "
        "sort -k 1,1n -k 2,2n | cut -f 1-3 | bedtools merge -i - > {output}"

rule compress_DF_vcf:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf.gz"
    shell:
        "bcftools view -Oz -o {output} {input} ;"
        "htsfile {output} ;"
        "bcftools index {output}"

rule filter_vcf_with_bed:
    threads: 1
    input:
        bed = "MRCA_ref_mapping/{ref}.{step}.merge.bed",
        vcf = "MRCA_ref_mapping/{ref}/{sample}.{step}_" \
            "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf.gz"
    output:
        vcf = "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.vcf.gz",
        csi = "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.vcf.gz.csi"
    shell:
        "bcftools view {input.vcf} -T ^{input.bed} -Oz -o {output.vcf} ;"
        "bcftools index {output.vcf} -o {output.csi}"

rule filter_hsnps:
    threads: 1
    params:
        snp_cutoff = 0.9
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter.vcf.gz"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input}; "
        "bcftools index {output}"

rule extract_sample_VCF_SNPs:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter.hvar.vcf.gz"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter_hvar.consensus.fasta"
    shell:
        "bcftools consensus -i 'type=\"SNP\"' "
        "-f {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        " {input} > {output}; "

rule reheader_gubbins:
    threads: 1
    params:
        execdir = exec_dir
    conda: "conda_envs/phy_plots.yaml"
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter_hvar.consensus.fasta"
    output:
        "gubbins/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
        "BedFilter_hvar.consensus.fasta"
    shell:
        "{params.execdir}/scripts/reheader_gubbins.py {input} > {output}"

"""
10215* 2/1/2021 17:31  cat ./mabscessus/*.fasta > mabs.fasta
10218* 2/1/2021 17:31  run_gubbins.py --prefix mabs --min_snps 20 --threads 8 mabs.fasta
10226* 2/1/2021 17:53  cat mabs.filtered_polymorphic_sites.fasta
10239* 2/1/2021 18:06  ../../analysis/make_bed_from_gubbins.py mabs.recombination_predictions.embl | sort -k 1,1n -k 2,2n > mabs.gubbins.bed
10243* 2/1/2021 18:10  conda activate mabs
10260* 2/1/2021 18:17  bcftools view ../MRCA_ref_mapping/mabscessus.RG_SC_RA.merge.vcf.gz -T ^mabs.gubbins.bed -Oz -o mabs.merge.gubbins_filtered.vcf.gz ; bcftools index mabs.merge.gubbins_filtered.vcf.gz -o mabs.merge.gubbins_filtered.vcf.gz.csi
10262* 2/1/2021 18:56  bcftools view mabs.merge.gubbins_filtered.vcf.gz | ../../scripts/make_SNP_alignment.py > mabs.test.filtered_gubbins.fasta
10263* 2/1/2021 18:57  conda activate phyplots
10264* 2/1/2021 18:57  snp-dists -b -c mabs.test.filtered_gubbins.fasta > mabs.test.dist.csv
10265* 2/1/2021 18:58  snp-dists -b -c mabs.fasta > /dev/null
10266* 2/1/2021 18:59  ../../scripts/pairSNV_vs_time.py mabs.test.dist.csv /home/nick/Dropbox/MA/sample.csv
""" # noqa
# FIX THJE REFERENCE ACCESSION in BED file

# output of the reheader_gubbins rule to be concatenated into {ref}.fasta
# produced followingoutput, not all is needed
# {ref}.filtered_polymorphic_sites.fasta
# {ref}.filtered_polymorphic_sites.phylip
# {ref}.final_tree.tre
# {ref}.node_labelled.final_tree.tre
# {ref}.per_branch_statistics.csv
# {ref}.recombination_predictions.embl
# {ref}.recombination_predictions.gff
# {ref}.summary_of_snp_distribution.vcf
rule run_gubbins:
    threads: 8
    params:
        execdir = exec_dir,
        prefix = "gubbins/{ref}"
    input:
        "gubbins/{ref}.fasta"
    output:
        filtered_fasta = "gubbins/{ref}.filtered_polymorphic_sites.fasta",
        embl = "gubbins/{ref}.recombination_predictions.embl"
    shell:
        "run_gubbins.py --prefix {params.prefix} --min_snps 20 "
        "--threads {threads} {input}"

rule make_gubbins_bed:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "gubbins/{ref}.recombination_predictions.embl"
    output:
        "gubbins/{ref}.gubbins.bed"
    shell:
        "{params.execdir}/analysis/make_bed_from_gubbins.py {input} > {output}"


rule merge_vcf:
    threads: 8
    input:
        vcfs = lambda wildcards: glob.glob(
            "MRCA_ref_mapping/{ref}/*.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
            "BedFilter.hvar.vcf.gz".format(
                ref=wildcards.ref, step=wildcards.step
            )
        )
    output:
        "MRCA_ref_mapping/{ref}.{step}.merge.vcf.gz"
    shell:
        "bcftools merge --threads {threads} {input.vcfs} | "
        "bcftools view - -Oz -o {output}"

# filter merged vcf with gubbins bed
#"bcftools view {input.vcf} -T ^{input.bed} -Oz -o {output.vcf} ;"
#"bcftools index {output.vcf} -o {output.csi}"

rule make_SNP_alignment:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}.{step}.merge.vcf.gz"
    output:
        "SNP_phylo/{ref}.{step}.merge.fasta",
    shell:
        "bcftools view {input} | "
        "{params.execdir}/scripts/make_SNP_alignment.py > {output}"


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
        "SNP_phylo/{ref}.{step}.merge.fasta"
    output:
        tree = "SNP_phylo/{ref}.{step}.merge.fasta.treefile"
    shell:
        "iqtree -s {input} -nt {threads} -m MFP -B 10000 -redo"
