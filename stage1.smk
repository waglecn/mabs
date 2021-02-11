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
            reads.split('_')[-1].endswith('1.fastq.gz')
        ]
    elif read == 'R2':
        fastqs = [
            reads for reads in sample_dict[sample_name] if
            reads.split('_')[-1].endswith('2.fastq.gz')
        ]
    else:
        exit('Uh oh.')

    reads = sorted([
        os.path.abspath(
            os.path.join(exec_dir, config['data_dir'], fastq)
        ) for fastq in fastqs
    ])
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
        adapters = os.path.join(exec_dir, config['adapters_fasta']),
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
    params:
        execdir = exec_dir
    input:
        "trimmed_input/{sample}.log"
    output:
        "trimmed_input/{sample}.csv"
    shell:
        "{params.execdir}/scripts/trim_csv.py {input} > {output}"


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
        kdb = os.path.join(exec_dir, config['kraken_db']),
    log:
        "trimmed_kraken/{sample}.trimmed.paired.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed --paired {input.R1} {input.R2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# QC summary csv
rule kraken2_summary:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "trimmed_kraken/{sample}.trimmed.paired"
    output:
        "trimmed_kraken/{sample}.csv"
    shell:
        "{params.execdir}/scripts/kraken_summary.py {input} S > {output}"


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
        kdb = os.path.join(exec_dir, config['kraken_db'])
    log:
        "trimmed_kraken/{sample}.trimmed.single.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed {input.S1} {input.S2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# Produces a basic assembly using raw reads for MRCA and species typing
rule shovill_assembly:
    threads: 7
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
    params:
        execdir = exec_dir
    input:
        "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        "shovill_assembly/{sample}.shovill.csv"
    shell:
        "{params.execdir}/scripts/assembly_summary.py {input} > {output}"

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
        kdb = os.path.join(exec_dir, config["kraken_db"]),
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
        query_file = os.path.join(exec_dir, config["erm41_query"])
    shell:
        "tblastn -subject {input} -query {params.query_file} -outfmt 6 "
        "| head -n 1 | cut -f 2,3,4,6,7,8,9,10 > {output}"

# QC csv for erm41 status
rule erm41_csv:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "erm41_status/{sample}.erm41.status"
    output:
        "erm41_status/{sample}.erm41.csv"
    shell:
        "{params.execdir}/scripts/erm41_status.py {input} > {output}"

# Everything gets mapped to ATCC19977 to determine basic coverage for QC
# Consider splitting this up to use threads better
rule map_to_mabs:
    threads: 7
    input:
        R1 = "trimmed_input/{sample}.R1.fastq.gz",
        R2 = "trimmed_input/{sample}.R2.fastq.gz",
        S1 = "trimmed_input/{sample}.S1.fastq.gz",
    params:
        execdir = exec_dir
    output:
        paired_temp = temp("ref_mapping/{sample}.paired.sorted.bam"),
        paired_temp_bai = temp("ref_mapping/{sample}.paired.sorted.bam.bai"),
        single_temp = temp("ref_mapping/{sample}.single.sorted.bam"),
        single_temp_bai = temp("ref_mapping/{sample}.single.sorted.bam.bai"),
        merge_temp = temp("ref_mapping/{sample}.merge.bam"),
        merge_sorted = "ref_mapping/mabs/{sample}.merged.sorted.bam",
        merge_sorted_bai = "ref_mapping/mabs/{sample}.merged.sorted.bam.bai"
    shell:
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/mabscessus "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/mabscessus "
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
    params:
        execdir = exec_dir
    output:
        csv = "ref_mapping/mabs/{sample}.merged.sorted.csv",
        depth = "ref_mapping/mabs/{sample}.merged.sorted.depth"
    shell:
        "bedtools genomecov -d -ibam {input} -g "
        "{params.execdir}/resources/alignment_references/mabcessus.fasta > "
        "{output.depth} ; "
        "echo -e \"$(samtools view -c -F 2308 {input}),"
        "$({params.execdir}/scripts/depth_summary.py {output.depth}),"
        "$({params.execdir}/scripts/count_softclips.py {input})\" > "
        "{output.csv}"

# This rule puts everything in context of a rough mashtree to determine MRCA
# for downstream rules - MLST, specific references
rule mashtree_assemblies:
    threads: 7
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
    params:
        execdir = exec_dir
    input:
        tree = "mashtree/assembly_mashtree.complete.tree",
        assembly = "shovill_assembly/{sample}.shovill.contigs.fa"
    output:
        "MRCA_MLST/{sample}.mlst.txt"
    shell:
        # "echo $(scripts/tree_MRCA.py {input.tree} {sample}) ;"
        "mlst --scheme $({params.execdir}/scripts/tree_MRCA.py {input.tree} "
        "{wildcards.sample}) "
        "--threads {threads} {input.assembly} > {output}"



rule QC_stats_per_sample:
    input:
        trim = "trimmed_input/{sample}.csv",
        kraken2 = "trimmed_kraken/{sample}.csv",
        assembly = "shovill_assembly/" \
            "{sample}.shovill.csv",
        mlst = "MRCA_MLST/{sample}.mlst.txt",
        erm41 = "erm41_status/{sample}.erm41.csv",
        mabs_depth = "ref_mapping/mabs/" \
            "{sample}.merged.sorted.csv",
        contigs = "shovill_assembly/" \
            "{sample}.shovill.contigs.fa",
        tree = "mashtree/assembly_mashtree.complete.tree",
    output:
        "QC_summary/{sample}.QC.csv"
    params:
        # need lambda function to trigger evaluation of wildcards this node is
        # in the DAG
        MRCA = lambda wildcards: MRCA_mapped_ref_input(wildcards.sample)
    shell:
        "echo -e \"{wildcards.sample},$(cat {input.trim}),"
        "$(cat {input.kraken2}),"
        "$(cat {input.assembly}),$(cat {input.mlst} | tr \"\t\" \",\"),"
        "$(cat {input.erm41}),"
        "{params.MRCA},$(cat {input.mabs_depth})\" | tee {output}"


rule merge_QC_summary:
    threads: 1
    input:
        expand("QC_summary/{sample}.QC.csv", sample=sample_names)
    output:
        "QC_summary.csv"
    params:
        execdir = exec_dir
    shell:
        "{params.execdir}/scripts/merge_QC_csv.py {output} {input}"
