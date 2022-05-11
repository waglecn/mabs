def input_files(sample_name, read):
    '''
    This function parses the sample sheet and provides the sample_names
    list to the expand functions
        input:
            string: sample_name
            string: 'R1' or 'R2'
        returns:
            list: sorted strings of R1 or R2 filenames for the sample reads

    sys exit if something doesn't match
    '''
    try:
        if read == 'MINON':
            fastqs = samples[sample_name][read]
        elif read == 'R1' or read == 'R2':
            fastqs = samples[sample_name]['MISEQ'][read]
    except IndexError:
        exit('Uh oh.')

    reads = sorted(
        [
            config['data_dir'] + r for r in fastqs
        ]
    )
    if len(reads) == 0:
        exit('no reads found')
    return reads

########################################################################
# Internal Data Rules - PHASE 1
########################################################################

# This rule looks like the merge rule from SIGNAL
# unifies input filenames to make wildcard parsing easier
# TODO - look into symlinks so only the samples with multiple read sets
# need to be merged, symlink others with expected output names in order to save
# disk/time
rule merge_ill_samples_only:
    input:
        expand(
            "{res}/{s}/input/{R}.fastq.gz",
            res=res,
            s=sample_names, R=["R1", "R2"]
        )

rule merge_long_samples_only:
    input:
        expand(
            "{res}/{s}/input/long.fastq.gz",
            res=res,
            s=sample_names, R=['MINION']
        )

rule merge_samples:
    conda:
        "envs/fastqc.yaml"
    params:
        execdir = exec_dir
    input:
        lambda wildcards: input_files(wildcards.sample, wildcards.R)
    output:
        f"{res}/{{sample}}/input/{{R}}.fastq.gz",
    shell:
        "zcat -f {input} | gzip > {output}"



# # FastQC on raw reads
# rule pre_trim_QC:
#     threads: 1
#     conda:
#         "conda_envs/fastqc.yaml"
#     input:
#         "merged_input/{sample}.{R}.fastq.gz",
#     output:
#         "pre_trim_QC/{sample}.{R}_fastqc.html"
#     shell:
#         "fastqc --noextract  -t {threads} -o pre_trim_QC {input}"


# # Trimming of raw (merged) reads
# rule trim_raw_reads:
#     conda:
#         "conda_envs/fastqc.yaml"
#     threads: 4
#     input:
#         R1 = "merged_input/{sample}.R1.fastq.gz",
#         R2 = "merged_input/{sample}.R2.fastq.gz"
#     output:
#         tR1 = "trimmed_input/{sample}.R1.fastq.gz",
#         tR2 = "trimmed_input/{sample}.R2.fastq.gz",
#         tS1 = "trimmed_input/{sample}.S1.fastq.gz",
#         tS2 = "trimmed_input/{sample}.S2.fastq.gz",
#     params:
#         adapters = os.path.join(exec_dir, config['adapters_fasta']),
#         min_len = config['QC_min_length']
#     log:
#         "trimmed_input/{sample}.log"
#     shell:
#         "trimmomatic PE -phred33 -threads {threads} {input.R1} {input.R2} "
#         "{output.tR1} {output.tS1} {output.tR2} {output.tS2} "
#         "ILLUMINACLIP:{params.adapters}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:"
#         "{params.min_len} 2>&1 | tee {log}"

# # FastQC on trimmed_reads
# rule post_trim_QC:
#     threads: 1
#     conda:
#         "conda_envs/fastqc.yaml"
#     input:
#         "trimmed_input/{sample}.{read}.fastq.gz"
#     output:
#         "post_trim_QC/{sample}.{read}_fastqc.html"
#     shell:
#         "fastqc --noextract -t {threads} -o post_trim_QC {input}"

# # Kraken2 contamination on the merged raw paired input
# rule kraken2_contamination_paired:
#     threads: 1
#     conda:
#         "conda_envs/kraken2.yaml"
#     input:
#         R1 = "trimmed_input/{sample}.R1.fastq.gz",
#         R2 = "trimmed_input/{sample}.R2.fastq.gz"
#     output:
#         report = "trimmed_kraken/{sample}.trimmed.paired",
#     params:
#         kdb = os.path.join(exec_dir, config['kraken_db']),
#     log:
#         "trimmed_kraken/{sample}.trimmed.paired.log"
#     shell:
#         "kraken2 --db {params.kdb} --threads {threads} --quick "
#         "--gzip-compressed --paired {input.R1} {input.R2} "
#         "--report {output.report} --output - 2>&1 | tee {log}"

# # Detection of post-trim contamination in single reads that will be used for
# # mapping
# rule kraken2_contamination_single:
#     threads: 1
#     conda:
#         "conda_envs/kraken2.yaml"
#     input:
#         S1 = "trimmed_input/{sample}.S1.fastq.gz",
#         S2 = "trimmed_input/{sample}.S2.fastq.gz"
#     output:
#         report = "trimmed_kraken/{sample}.trimmed.single",
#     params:
#         kdb = os.path.join(exec_dir, config['kraken_db'])
#     log:
#         "trimmed_kraken/{sample}.trimmed.single.log"
#     shell:
#         "kraken2 --db {params.kdb} --threads {threads} --quick "
#         "--gzip-compressed {input.S1} {input.S2} "
#         "--report {output.report} --output - 2>&1 | tee {log}"

# # Produces a basic assembly using raw reads for MRCA and species typing
# rule shovill_assembly:
#     threads: 4
#     conda:
#         "conda_envs/shovill.yaml"
#     input:
#         R1 = "merged_input/{sample}.R1.fastq.gz",
#         R2 = "merged_input/{sample}.R2.fastq.gz"
#     output:
#         "shovill_assembly/{sample}.shovill.contigs.fa",
#     params:
#         shovill_tmp = config['shovill_tempdir'],
#         outdir = "shovill_assembly/{sample}"
#     log:
#         "shovill_assembly/{sample}.shovill.log"
#     shell:
#         # Note: letting Shovill set memory caps at 8, have gotten strange out
#         # of memory errors in that pipeline see:
#         # <https://github.com/tseemann/shovill/issues/59>
#         "shovill --force --cpus {threads} --R1 {input.R1} --R2 {input.R2} "
#         "--outdir {params.outdir} --tmpdir {params.shovill_tmp} --ram 16 "
#         " 2>&1 | tee {log} ; cp {params.outdir}/contigs.fa {output}"

# # Kraken2 contamination of assembly
# rule kraken2_contamination_assembly:
#     threads: 1
#     conda:
#         "conda_envs/kraken2.yaml"
#     input:
#         "shovill_assembly/{sample}.shovill.contigs.fa"
#     output:
#         report = "assembly_kraken/{sample}.assembly"
#     params:
#         kdb = os.path.join(exec_dir, config["kraken_db"]),
#     log:
#         "assembly_kraken/{sample}.kraken.log"
#     shell:
#         "kraken2 --db {params.kdb} --threads {threads} --quick "
#         "{input} --report {output.report} --output - 2>&1 | tee {log}"

# # Basic BLAST detection of Erm(41) - looking for intact or truncated CDS
# rule assembly_status_erm41:
#     threads: 1
#     conda:
#         "conda_envs/blast.yaml"
#     input:
#         "shovill_assembly/{sample}.shovill.contigs.fa"
#     output:
#         "erm41_status/{sample}.erm41.status"
#     params:
#         query_file = os.path.join(exec_dir, config["erm41_query"])
#     shell:
#         "tblastn -subject {input} -query {params.query_file} -outfmt 6 "
#         "| head -n 1 | cut -f 2,3,4,6,7,8,9,10 > {output}"

# # Everything gets mapped to ATCC19977 to determine basic coverage for QC
# # Consider splitting this up to use threads better
# rule map_to_mabs:
#     threads: 4
#     input:
#         R1 = "trimmed_input/{sample}.R1.fastq.gz",
#         R2 = "trimmed_input/{sample}.R2.fastq.gz",
#         S1 = "trimmed_input/{sample}.S1.fastq.gz",
#     params:
#         execdir = exec_dir
#     output:
#         paired_temp = temp("ref_mapping/{sample}.paired.sorted.bam"),
#         paired_temp_bai = temp("ref_mapping/{sample}.paired.sorted.bam.bai"),
#         single_temp = temp("ref_mapping/{sample}.single.sorted.bam"),
#         single_temp_bai = temp("ref_mapping/{sample}.single.sorted.bam.bai"),
#         merge_temp = temp("ref_mapping/{sample}.merge.bam"),
#         merge_sorted = "ref_mapping/mabs/{sample}.merged.sorted.bam",
#         merge_sorted_bai = "ref_mapping/mabs/{sample}.merged.sorted.bam.bai"
#     shell:
#         "bwa mem -t {threads} -M "
#         "{params.execdir}/resources/alignment_references/mabscessus "
#         "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
#         "{output.paired_temp} ;"
#         "bwa mem -t {threads} -M "
#         "{params.execdir}/resources/alignment_references/mabscessus "
#         "{input.S1} | samtools view -Sbh - | samtools sort -@ {threads} > "
#         "{output.single_temp} ;"
#         "samtools index {output.paired_temp} ;"
#         "samtools index {output.single_temp} ;"
#         "samtools merge -f {output.merge_temp} {output.paired_temp} "
#         "{output.single_temp} ;"
#         "samtools view -hF 4 {output.merge_temp} | samtools sort -@ {threads} > "
#         "{output.merge_sorted} ;"
#         "samtools index {output.merge_sorted} {output.merge_sorted_bai}"

# rule depth_map_to_mabs:
#     threads: 1
#     input:
#         "ref_mapping/mabs/{sample}.merged.sorted.bam"
#     params:
#         execdir = exec_dir
#     output:
#         depth = "ref_mapping/mabs/{sample}.merged.sorted.depth.gz"
#     shell:
#         "bedtools genomecov -d -ibam {input} -g "
#         "{params.execdir}/resources/alignment_references/mabcessus.fasta | gzip > "
#         "{output.depth} "

# # This rule puts everything in context of a rough mashtree to determine MRCA
# # for downstream rules - MLST, specific references
# rule mashtree_assemblies:
#     threads: 8
#     conda:
#         "conda_envs/mashtree.yaml"
#     input:
#         # this is ugly - but there are fewer constraints on input
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

# rule MRCA_from_tree:
#     threads: 1
#     conda:
#         "conda_envs/ngd_phylo.yaml"
#     input:
#         "mashtree/assembly_mashtree.complete.tree"
#     output:
#         "MRCA/{sample}.csv"
#     params:
#         MRCA = lambda wildcards: MRCA_mapped_ref_input(wildcards.sample)
#     shell:
#         "echo {params.MRCA} > {output}"


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


# rule QC_stats_per_sample:
#     conda: "conda_envs/phy_plots.yaml"
#     input:
#         tree = "mashtree/assembly_mashtree.complete.tree",
#         trim = "trimmed_input/{sample}.log",
#         kraken2 = "trimmed_kraken/{sample}.trimmed.paired",
#         assembly = "shovill_assembly/{sample}.shovill.contigs.fa",
#         mlst = "MRCA_MLST/{sample}.mlst.txt",
#         erm41 = "erm41_status/{sample}.erm41.status",
#         mabs_bam = "ref_mapping/mabs/{sample}.merged.sorted.bam",
#         mabs_depth = "ref_mapping/mabs/{sample}.merged.sorted.depth.gz",
#         MRCA = "MRCA/{sample}.csv"
#     output:
#         "QCsummary/{sample}.QC.csv"
#     params:
#         execdir = exec_dir
#     shell:
#         "{params.execdir}/scripts/make_sample_QC_csv.py {wildcards.sample} "
#         "{input.trim} {input.kraken2} {input.assembly} {input.erm41} "
#         "{input.mabs_bam} {input.mabs_depth} {input.mlst} {input.MRCA} > {output}"


# rule merge_QC_summary:
#     priority: 10
#     threads: 1
#     input:
#         expand("QCsummary/{sample}.QC.csv", sample=sample_names)
#     output:
#         "QC_summary.csv"
#     params:
#         execdir = exec_dir
#     shell:
#         "{params.execdir}/scripts/merge_QC_csv.py {output} {input}"
