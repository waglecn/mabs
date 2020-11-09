import pandas
import glob

configfile: "config.yaml"

# samples = pandas.read_table(config['sample_sheet'], sep=',', header=None)

# see line above, this is an ugly way to load the sample sheet, pandas
# will definitely provide more validation
samples = [
    line.strip().split(',') for line in open(config['sample_sheet'], 'r')
]
sample_dict = {sample[0]: sample[1:] for sample in samples}
# print(sample_dict)
sample_names = sample_dict.keys()

filt_samples = list(sample_names)
excluded_samples = [
    '29-3L', '29-3K', '29-2Z', '29-2L', '29-G', '27-D', '25-U', '25-J', '20-E',
    '19-A', '19-B', '17-F_2', '17-E', '17-B_combined', '16-A', '14-U', '14-F',
]
for name in excluded_samples:
    filt_samples.remove(name)


# helper function definitions
def input_files(sample_name, read):
    '''
    This function parses the sample sheet csv and provides the sample_names
    list to the expand functions
        input:
            string: sample_name
            string: 'R1' or 'R2'
        returns:
            list: sorted strings of R1 or R2 filenames for the sample reads

    sys exit if something doesn't match
    '''
    if read == 'R1':
        fastqs = [
            reads for reads in sample_dict[sample_name] if
            reads.endswith('_R1.fastq.gz')
        ]
    elif read == 'R2':
        fastqs = [
            reads for reads in sample_dict[sample_name] if
            reads.endswith('_R2.fastq.gz')
        ]
    else:
        exit('Uh oh.')
    # print(fastqs)
    return sorted([
        os.path.abspath(
            os.path.join(config['data_dir'], fastq)
        ) for fastq in fastqs
    ])


def MRCA_mapped_ref_input(sample_name):
    '''
    Use as input function that maps sample_name to reference
        input:
            sample_name: key should be found in sample_dict
        returns:
            ref_name: should match available references
    '''
    from scripts.tree_MRCA import tree_MRCA
    return tree_MRCA('mashtree/assembly_mashtree.complete.tree', sample_name)


# Rule "all" default catches output of other rules as input in order to
# simplify running the workflow
# - add target rules here as they are implemented
rule all:
    input:
        ###############
        # STAGE 1 - QC
        ###############
        # "merged reads"
        expand(
            "merged_input/{sample}.{R}.fastq.gz", sample=sample_names,
            R=["R1", "R2"]
        ),

rule stage1:
    input:
        # pre-trim QC
        expand(
            "pre_trim_QC/{sample}.{R}_fastqc.html", sample=sample_names,
            R=["R1", "R2"]
        ),
        # trimmed input
        expand(
            "trimmed_input/{sample}.{read}.fastq.gz", sample=sample_names,
            read=["R1", "R2", "S1", "S2"]
        ),
        expand("trimmed_input/{sample}.csv", sample=sample_names),
        # post-trim QC
        expand(
            "post_trim_QC/{sample}.{read}_fastqc.html", sample=sample_names,
            read=["R1", "R2", "S1", "S2"]
        ),
        # Kraken contamination, paired and single
        expand(
            "trimmed_kraken/{sample}.trimmed.paired", sample=sample_names
        ),
        expand(
            "trimmed_kraken/{sample}.trimmed.single", sample=sample_names
        ),
        expand("trimmed_kraken/{sample}.csv", sample=sample_names),
        # assembly with shovill
        expand(
            "shovill_assembly/{sample}.shovill.contigs.fa", sample=sample_names
        ),
        expand("shovill_assembly/{sample}.shovill.csv", sample=sample_names),
        # kraken2 the assembly
        expand("assembly_kraken/{sample}.assembly", sample=sample_names),
        # determine the basic Erm(41) status with TBLASTN
        expand("erm41_status/{sample}.erm41.status", sample=sample_names),
        expand("erm41_status/{sample}.erm41.csv", sample=sample_names),
        # map all samples to mabs to determine basic coverage
        expand(
            "ref_mapping/mabs/{sample}.merged.sorted.bam",
            sample=sample_names
        ),
        expand(
            "ref_mapping/mabs/{sample}.merged.sorted.csv", sample=sample_names
        ),
        "mashtree/assembly_mashtree.complete.tree",
        "mashtree/assembly_mashtree.complete.matrix",
        # type assemblies using MRCAs from mashtree
        expand("MRCA_MLST/{sample}.mlst.txt", sample=sample_names),
        expand("QC_summary/{sample}.QC.csv", sample=sample_names)

rule QC_stats_per_sample:
    input:
        trim = "trimmed_input/{sample}.csv",
        kraken2 = "trimmed_kraken/{sample}.csv",
        assembly = "shovill_assembly/{sample}.shovill.csv",
        mlst = "MRCA_MLST/{sample}.mlst.txt",
        erm41 = "erm41_status/{sample}.erm41.csv",
        mabs_depth = "ref_mapping/mabs/{sample}.merged.sorted.csv",
        contigs = "shovill_assembly/{sample}.shovill.contigs.fa",
        tree = "mashtree/assembly_mashtree.complete.tree",
    output:
        "QC_summary/{sample}.QC.csv"
    params:
        # need lambda function to trigger evaluation of wildcards this node is
        # in the DAG
        MRCA = lambda wildcards: MRCA_mapped_ref_input(wildcards.sample)
    shell:
        "echo -e \"{wildcards.sample},$(cat {input.trim}),$(cat {input.kraken2}),"
        "$(cat {input.assembly}),$(cat {input.mlst} | tr \"\t\" \",\"),"
        "$(cat {input.erm41}),"
        "{params.MRCA},$(cat {input.mabs_depth})\" | tee {output}"

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
            ) for s in sample_names
        ],
        # gatk realignment
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_MQ30_BAQ.mpileup".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],

rule with_realignment:
    input:
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.failed.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}"
            ".RG_SC_RA_MQ30_BAQ50_DP20_SP60_AD1.AD_failed.vcf".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.0cov.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        expand(
            "resources/alignment_references/{ref}.PE_PPE.bed",
            ref=['mabscessus', 'mbolettii', 'mmassiliense']
        ),
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.merge.bed".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],
        [
            "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_"
            "MQ30_BAQ50_DP20_SP60_AD1_DF_BedFilter.hvar.vcf.gz".format(
                ref=MRCA_mapped_ref_input(s), sample=s
            ) for s in sample_names
        ],

rule SNP_alignment:
    input:
        expand(
            "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.vcf.gz",
            ref=['mabscessus', 'mbolettii', 'mmassiliense']
        ),
        expand(
            "SNP_phylo/{ref}.RG_SC_RA.merge.fasta",
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

        # build SNP alignment


        # generate SNP tree for each sub-species type

#########################################################################
# External Data Rules
#########################################################################
checkpoint refseq_outgroup_download:
    conda:
        'conda_envs/ngd_phylo.yaml'
    params:
        outgroup = config['outgroup_assembly']
    output:
        directory("external/refseq/outgroup/")
    log:
        "external/refseq_outgroup_download.log"
    shell:
        "ngd -s refseq -r3 --flat-output -A {params.outgroup} -F genbank "
        "-o {output} bacteria ; "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"

rule download_refseq_taxlist:
    params:
        taxid = config['target_species_taxid']
    output:
        "external/taxid{}.target.taxlist".format(
            config['target_species_taxid']
        )
    shell:
        "python scripts/gimme-taxa.py -j {params.taxid} > {output}"

checkpoint refseq_complete_genome_download:
    conda:
        'conda_envs/ngd_phylo.yaml'
    input:
        "external/taxid{}.target.taxlist".format(
            config['target_species_taxid']
        )
    log:
        "external/refseq_complete_download.log"
    output:
        directory("external/refseq/complete/"),
        # a non-dynamic output may not exist with dynamic output
    shell:
        "ngd -s refseq -r 3 -p 3 -o {output} "
        "-t {input} --flat-output -F genbank "
        "-m external/complete.refseq.meta -l complete bacteria ; "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"

checkpoint refseq_all_genome_download:
    conda:
        "conda_envs/ngd_phylo.yaml"
    input:
        "external/taxid{}.target.taxlist".format(
            config['target_species_taxid']
        )
    log:
        "external/refseq_all_download.log"
    output:
        directory("external/refseq/all/")
    shell:
        "ngd -s refseq -r 3 -p 3 -o {output} "
        "-t {input} --flat-output -F genbank "
        "-m external/all.refseq.meta bacteria "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"


def aggregate_refseq(wildcards):
    if wildcards.nodeset == 'complete':
        checkpoint_set = checkpoints.refseq_complete_genome_download.get(
            **wildcards
        ).output[0]
    elif wildcards.nodeset == 'all':
        checkpoint_set = checkpoints.refseq_all_genome_download.get(
            **wildcards
        ).output[0]
    checkpoint_outgroup = checkpoints.refseq_outgroup_download.get(
        **wildcards
    ).output[0]
    assembly_name = glob.glob(checkpoint_set + "*.gbk.gz")
    outgroup_name = glob.glob(checkpoint_outgroup + "*.gbk.gz")

    return assembly_name + outgroup_name


rule register_gatk:
    conda: "conda_envs/gatk3.yaml"
    output: 'resources/gatk-registered'
    params:
        gatk_jar = "resources/{gatk_jar}".format(gatk_jar=config['gatk3_jar'])
    shell:
        "echo {output} ; gatk3-register {params.gatk_jar} ; touch {output}"

# TODO rule for external short read data from SRA

########################################################################
# Internal Data Rules - PHASE 1
########################################################################

# This rule looks like the merge rule from SIGNAL
# unifies input filenames to make wildcard parsing easier
# TODO - look into symlinks so only the samples with multiple read sets
# need to be merged, symlink others with expected output names in order to save
# disk/time
rule merge_samples:
    conda:
        "conda_envs/fastqc.yaml"
    input:
        lambda wildcards: input_files(wildcards.sample, wildcards.R)
    output:
        "merged_input/{sample}.{R}.fastq.gz",
    shell:
        "zcat -f {input} | gzip > {output}"

# FastQC on raw reads
rule pre_trim_QC:
    threads: 1
    conda:
        "conda_envs/fastqc.yaml"
    input:
        "merged_input/{sample}.{R}.fastq.gz",
    output:
        "pre_trim_QC/{sample}.{R}_fastqc.html"
    shell:
        "fastqc --noextract  -t {threads} -o pre_trim_QC {input}"


# Trimming of raw (merged) reads
rule trim_raw_reads:
    conda:
        "conda_envs/fastqc.yaml"
    threads: 4
    input:
        R1 = "merged_input/{sample}.R1.fastq.gz",
        R2 = "merged_input/{sample}.R2.fastq.gz"
    output:
        tR1 = "trimmed_input/{sample}.R1.fastq.gz",
        tR2 = "trimmed_input/{sample}.R2.fastq.gz",
        tS1 = "trimmed_input/{sample}.S1.fastq.gz",
        tS2 = "trimmed_input/{sample}.S2.fastq.gz",
    params:
        adapters = config['adapters_fasta'],
        min_len = config['QC_min_length']
    log:
        "trimmed_input/{sample}.log"
    shell:
        "trimmomatic PE -threads {threads} {input.R1} {input.R2} "
        "{output.tR1} {output.tS1} {output.tR2} {output.tS2} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:"
        "{params.min_len} 2>&1 | tee {log}"

# Produces QC csv
rule post_trim_csv:
    threads: 1
    input:
        "trimmed_input/{sample}.log"
    output:
        "trimmed_input/{sample}.csv"
    shell:
        "scripts/trim_csv.py {input} > {output}"


# FastQC on trimmed_reads
rule post_trim_QC:
    threads: 4
    conda:
        "conda_envs/fastqc.yaml"
    input:
        "trimmed_input/{sample}.{read}.fastq.gz"
    output:
        "post_trim_QC/{sample}.{read}_fastqc.html"
    shell:
        "fastqc --noextract -t {threads} -o post_trim_QC {input}"

# Kraken2 contamination on the merged raw paired input
rule kraken2_contamination_paired:
    threads: 2
    conda:
        "conda_envs/kraken2.yaml"
    input:
        R1 = "trimmed_input/{sample}.R1.fastq.gz",
        R2 = "trimmed_input/{sample}.R2.fastq.gz"
    output:
        report = "trimmed_kraken/{sample}.trimmed.paired",
    params:
        kdb = config['kraken_db']
    log:
        "trimmed_kraken/{sample}.trimmed.paired.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed --paired {input.R1} {input.R2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# QC summary csv
rule kraken2_summary:
    threads: 1
    input:
        "trimmed_kraken/{sample}.trimmed.paired"
    output:
        "trimmed_kraken/{sample}.csv"
    shell:
        "scripts/kraken_summary.py {input} S > {output}"


# Detection of post-trim contamination in single reads that will be used for
# mapping
rule kraken2_contamination_single:
    threads: 2
    conda:
        "conda_envs/kraken2.yaml"
    input:
        S1 = "trimmed_input/{sample}.S1.fastq.gz",
        S2 = "trimmed_input/{sample}.S2.fastq.gz"
    output:
        report = "trimmed_kraken/{sample}.trimmed.single",
    params:
        kdb = config['kraken_db']
    log:
        "trimmed_kraken/{sample}.trimmed.single.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed {input.S1} {input.S2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# Produces a basic assembly using raw reads for MRCA and species typing
rule shovill_assembly:
    threads: 8
    conda:
        "conda_envs/shovill.yaml"
    input:
        R1 = "merged_input/{sample}.R1.fastq.gz",
        R2 = "merged_input/{sample}.R2.fastq.gz"
    output:
        "shovill_assembly/{sample}.shovill.contigs.fa",
    params:
        shovill_tmp = config['shovill_tempdir'],
        outdir = "shovill_assembly/{sample}"
    log:
        "shovill_assembly/{sample}.shovill.log"
    shell:
        # Note: letting Shovill set memory caps at 8, have gotten strange out
        # of memory errors in that pipeline see:
        # <https://github.com/tseemann/shovill/issues/59>
        "shovill --force --cpus {threads} --R1 {input.R1} --R2 {input.R2} "
        "--outdir {params.outdir} --tmpdir {params.shovill_tmp} --ram 16 "
        " 2>&1 | tee {log} ; cp {params.outdir}/contigs.fa {output}"

# QC on assembly output
rule shovill_assembly_summary:
    threads: 1
    conda:
        "conda_envs/biopython.yaml"
    input:
        "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        "shovill_assembly/{sample}.shovill.csv"
    shell:
        "scripts/assembly_summary.py {input} > {output}"

# Kraken2 contamination of assembly
rule kraken2_contamination_assembly:
    threads: 2
    conda:
        "conda_envs/kraken2.yaml"
    input:
        "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        report = "assembly_kraken/{sample}.assembly"
    params:
        kdb = config["kraken_db"]
    log:
        "assembly_kraken/{sample}.kraken.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "{input} --report {output.report} --output - 2>&1 | tee {log}"

# Basic BLAST detection of Erm(41) - looking for intact or truncated CDS
rule assembly_status_erm41:
    threads: 1
    conda:
        "conda_envs/blast.yaml"
    input:
        "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        "erm41_status/{sample}.erm41.status"
    params:
        query_file = config["erm41_query"]
    shell:
        "tblastn -subject {input} -query {params.query_file} -outfmt 6 "
        "| head -n 1 | cut -f 2,3,4,6,7,8,9,10 > {output}"

# QC csv for erm41 status
rule erm41_csv:
    threads: 1
    input:
        "erm41_status/{sample}.erm41.status"
    output:
        "erm41_status/{sample}.erm41.csv"
    shell:
        "scripts/erm41_status.py {input} > {output}"

# Everything gets mapped to ATCC19977 to determine basic coverage for QC
# Consider splitting this up to use threads better
rule map_to_mabs:
    threads: 7
    input:
        R1 = "trimmed_input/{sample}.R1.fastq.gz",
        R2 = "trimmed_input/{sample}.R2.fastq.gz",
        S1 = "trimmed_input/{sample}.S1.fastq.gz",
    output:
        paired_temp = temp("ref_mapping/{sample}.paired.sorted.bam"),
        paired_temp_bai = temp("ref_mapping/{sample}.paired.sorted.bam.bai"),
        single_temp = temp("ref_mapping/{sample}.single.sorted.bam"),
        single_temp_bai = temp("ref_mapping/{sample}.single.sorted.bam.bai"),
        merge_temp = temp("ref_mapping/{sample}.merge.bam"),
        merge_sorted = "ref_mapping/mabs/{sample}.merged.sorted.bam",
        merge_sorted_bai = "ref_mapping/mabs/{sample}.merged.sorted.bam.bai"
    shell:
        "bwa mem -t {threads} -M resources/alignment_references/mabscessus "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa mem -t {threads} -M resources/alignment_references/mabscessus "
        "{input.S1} | samtools view -Sbh - | samtools sort > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort > "
        "{output.merge_sorted} ;"
        "samtools index {output.merge_sorted} {output.merge_sorted_bai}"

rule depth_map_to_mabs:
    threads: 1
    input:
        "ref_mapping/mabs/{sample}.merged.sorted.bam"
    output:
        csv = "ref_mapping/mabs/{sample}.merged.sorted.csv",
        depth = "ref_mapping/mabs/{sample}.merged.sorted.depth"
    shell:
        "bedtools genomecov -d -ibam {input} -g "
        "resources/alignment_references/mabcessus.fasta > {output.depth} ; "
        "echo -e \"$(samtools view -c -F 2308 {input}),"
        "$(scripts/depth_summary.py {output.depth}),"
        "$(scripts/count_softclips.py {input})\" > {output.csv}"


######################################################################
#
# STAGE 2 RULES
#
######################################################################

# This rule puts everything in context of a rough mashtree to determine MRCA
# for downstream rules - MLST, specific references
rule mashtree_assemblies:
    threads: 8
    conda:
        "conda_envs/mashtree.yaml"
    input:
        # this is ugly - but there are fewer contraints on input
        # wildcards in input must match wildcards in output, but since
        # this output rule has no wildcards, needed to something else
        # internal = glob.glob("shovill_assembly/*.shovill.contigs.fa"),
        # internal2 = glob.glob("../bakMA/4-shovill_assemble/*.contigs.fa"),
        # glob.glob("external/refseq/complete/*.gbff.gz")
        aggregate_refseq,
        internal = expand(
            "shovill_assembly/{sample}.shovill.contigs.fa",
            sample=sample_names
        )
    params:
        # the mashtree tempdir can get large, depending on the number of
        # species - see config.yaml
        tempdir = config['mashtree_tempdir']
    output:
        tree = "mashtree/assembly_mashtree.{nodeset}.tree",
        matrix = "mashtree/assembly_mashtree.{nodeset}.matrix"
    log:
        "mashtree/assembly_mashtree.{nodeset}.log"
    shell:
        "mashtree --mindepth 0 --numcpus {threads} --outtree {output.tree} "
        "--sort-order random "
        "--outmatrix {output.matrix} --tempdir {params.tempdir} "
        "{input} 2>&1 | tee {log}"

# This rule determines which samples will use which reference or MLST scheme
rule MRCA_MLST:
    threads: 1
    conda:
        "conda_envs/ngd_phylo.yaml"
    input:
        tree = "mashtree/assembly_mashtree.complete.tree",
        assembly = "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        "MRCA_MLST/{sample}.mlst.txt"
    shell:
        # "echo $(scripts/tree_MRCA.py {input.tree} {sample}) ;"
        "mlst --scheme $(scripts/tree_MRCA.py {input.tree} "
        "{wildcards.sample}) "
        "--threads {threads} {input.assembly} > {output}"

# This is modelled after map_to_mabs, except need to sub in the ref
# If the most appropriate ref is to mabs, this should prevent re-calculation
# of the alignment
# TODO - check the names output by the tree_MRCA.py
rule MRCA_ref_alignment:
    threads: 7
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
        "bwa mem -t {threads} -M resources/alignment_references/{params.ref} "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa mem -t {threads} -M resources/alignment_references/{params.ref} "
        "{input.S1} | samtools view -Sbh - | samtools sort > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort > "
        "{output.merge_sorted} ;"
        "samtools index {output.merge_sorted} {output.merge_sorted_bai}"

rule create_sam_read_groups:
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
    input:
        "MRCA_ref_mapping/{ref}/{sample}.RG.merged.sorted.bam"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam"
    shell:
        "scripts/sclips.py filter {input} > {output}"


rule create_fasta_dict_index:
    threads: 1
    conda: "conda_envs/picard.yaml"
    input:
        "resources/alignment_references/{ref}.fasta"
    output:
        seqdict = "resources/alignment_references/{ref}.dict",
        fai = "resources/alignment_references/{ref}.fasta.fai"
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
    input:
        bam = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam",
        bai = "MRCA_ref_mapping/{ref}/{sample}.RG_SC.bam.bai",
        seqdict = "resources/alignment_references/{ref}.dict",
        fai = "resources/alignment_references/{ref}.fasta.fai",
        gatk = "resources/gatk-registered"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals"
    log:
        "MRCA_ref_mapping/{ref}/{sample}.gatk3_intervals.log"
    shell:
        "gatk3 -T RealignerTargetCreator -R "
        "resources/alignment_references/{wildcards.ref}.fasta "
        "-I {input.bam} -o {output} 2>&1 | tee {log}"

rule MRCA_ref_gatk_realignment:
    threads: 1
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
        "-R resources/alignment_references/{wildcards.ref}.fasta "
        "-I {input.bam} -targetIntervals {input.intervals} -o {output} "
        " 2>&1 | tee {log}"

rule MRCA_make_mpileup:
    threads: 1
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
        "-f resources/alignment_references/{wildcards.ref}.fasta "
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
        "bcftools call --ploidy 1 -Ov -m -v {input} | "
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
        "bcftools filter -i 'SP<60 & ADF[0:1]>1 & ADR[0:1]>1 & MQ>30 & "
        "QUAL>50 & FORMAT/DP > 20 & SP<60 & ADF[0:1]>1 & ADR[0:1]>1' "
        "-o {output}"

rule MRCA_inverse_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_"
        "DP20_SP60_AD1.failed.vcf"
    shell:
        "bcftools call --ploidy 1 -O v -m -v {input} | "
        "bcftools filter -i 'SP>=60 || MQ<=30 || FORMAT/DP<=20 || QUAL<=50' "
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
        "bcftools call --ploidy 1 -O v -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff}) & (ADF[0:1]< =1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule density_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_MQ30_BAQ50_DP20_SP60_AD1.vcf"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_"
        "MQ30_BAQ50_DP20_SP60_AD1_DF.vcf"
    shell:
        "cp {input} {output}"

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
    conda:
        "conda_envs/biopython.yaml"
    input:
        "resources/alignment_references/{ref}.gbk"
    output:
        "resources/alignment_references/{ref}.PE_PPE.bed"
    shell:
        "scripts/make_PE_PPE_BED.py {input} > {output}"

rule merge_bed:
    threads: 1
    input:
        bed0 = "MRCA_ref_mapping/{ref}/{sample}.{step}.0cov.bed",
        PEPPEbed = "resources/alignment_references/{ref}.PE_PPE.bed"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.merge.bed"
    shell:
        "cat {input.bed0} {input.PEPPEbed} | "
        # this is bad, fix this
        "grep 'NC_010394\.1' -v | "
        "sort -k 1,1n -k 2,2n | bedtools merge -i - > {output}"

rule compress_vcf:
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
        bed = "MRCA_ref_mapping/{ref}/{sample}.{step}.merge.bed",
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
        "{params.snp_cutoff})' -Oz -o {output} {input}; bcftools index {output}"

rule merge_vcf:
    threads: 8
    input:
        vcfs = lambda wildcards: glob.glob(
            "MRCA_ref_mapping/{ref}/*.{step}_MQ30_BAQ50_DP20_SP60_AD1_DF_"
            "BedFilter.hvar.vcf.gz".format(ref=wildcards.ref, step=wildcards.step)
        )
    output:
        "MRCA_ref_mapping/{ref}.{step}.merge.vcf.gz"
    shell:
        "bcftools merge --threads {threads} {input} | "
        "bcftools view - -Oz -o {output}"

rule make_SNP_alignment:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}.{step}.merge.vcf.gz"
    output:
        "SNP_phylo/{ref}.{step}.merge.fasta",
    shell:
        "bcftools view {input} | scripts/make_SNP_alignment.py > {output}"


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
