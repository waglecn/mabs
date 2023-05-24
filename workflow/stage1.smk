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
    reads = []
    try:
        if (
            read == 'long' and
            'MINION' in samples[sample_name] and
            len(samples[sample_name]['MINION']) > 0
        ):
            fastqs = samples[sample_name]['MINION']
            reads = fastqs
        elif read == 'R1' or read == 'R2':
            fastqs = samples[sample_name]['MISEQ'][read]
            reads = sorted(
                [
                    config['data_dir'] + r for r in fastqs
                ]
            )
    except IndexError:
        exit('Uh oh.')

    if len(reads) == 0:
        print('{} {} - no reads found'.format(
            sample_name, read
        ), file=sys.stderr)
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
            s=short_sample_names, R=["R1", "R2"]
        )

rule merge_long_samples_only:
    input:
        expand(
            "{res}/{s}/input/long.fastq.gz",
            res=res,
            s=long_sample_names, R=['MINION']
        )

rule merge_samples:
    conda:
        "envs/fastqc.yaml"
    input:
        lambda wildcards: input_files(wildcards.sample, wildcards.R)
    output:
        f"{res}/{{sample}}/input/{{R,\w+}}.fastq.gz",
    shell:
        "zcat -f {input} | gzip > {output}"



# FastQC on raw reads
rule pre_trim_QC:
    threads: 1
    conda:
        "envs/fastqc.yaml"
    input:
        f"{res}/{{sample}}/input/{{R}}.fastq.gz",
    output:
        zipf = f"{res}/{{sample}}/input/{{R}}_fastqc.zip",
        html = f"{res}/{{sample}}/input/{{R}}_fastqc.html",
        data = f"{res}/{{sample}}/input/{{R}}_fastqc.data.txt"
    wildcard_constraints:
        R ="R\d"
    shell:
        "fastqc --noextract  -t {threads} "
        f"-o {res}/{{wildcards.sample}}/input {{input}} ; "
        f"unzip -p {res}/{{wildcards.sample}}/input/{{wildcards.R}}_fastqc.zip "
        f"{{wildcards.R}}_fastqc/fastqc_data.txt > {{output.data}}"


# Trimming of raw (merged) short reads
rule trim_raw_reads:
    conda:
        "envs/fastqc.yaml"
    threads: 4
    input:
        R1 = f"{res}/{{sample}}/input/R1.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.fastq.gz"
    output:
        tR1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        tR2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz",
        tS1 = f"{res}/{{sample}}/input/S1.trim.fastq.gz",
        tS2 = f"{res}/{{sample}}/input/S2.trim.fastq.gz",
    params:
        adapters = config['adapters_fasta'],
        min_len = config['QC_min_length']
    log:
        f"{res}/{{sample}}/input/trimmomatic.log"
    shell:
        "trimmomatic PE -phred33 -threads {threads} {input.R1} {input.R2} "
        "{output.tR1} {output.tS1} {output.tR2} {output.tS2} "
        "ILLUMINACLIP:{params.adapters}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:"
        "{params.min_len} 2>&1 | tee {log}"

# FastQC on trimmed short reads
rule post_trim_QC:
    threads: 1
    conda:
        "envs/fastqc.yaml"
    input:
        f"{res}/{{sample}}/input/{{R}}.trim.fastq.gz"
    output:
        zipf = f"{res}/{{sample}}/input/{{R}}.trim_fastqc.zip",
        html = f"{res}/{{sample}}/input/{{R}}.trim_fastqc.html",
        data = f"{res}/{{sample}}/input/{{R}}.trim_fastqc.data.txt"
    wildcard_constraints:
        R = "[RS]\d"
    shell:
        "fastqc --noextract -t {threads} "
        f"-o {res}/{{wildcards.sample}}/input {{input}} ; "
        f"unzip -p {res}/{{wildcards.sample}}/input/{{wildcards.R}}.trim_fastqc.zip "
        f"{{wildcards.R}}.trim_fastqc/fastqc_data.txt > {{output.data}}"


rule pre_nano_QC:
    threads: 1
    conda:
        "envs/fastqc.yaml"
    input:
        f"{res}/{{sample}}/input/long.fastq.gz"
    output:
        zipf = f"{res}/{{sample}}/input/long_fastqc.html",
        html = f"{res}/{{sample}}/input/long_fastqc.zip",
        data = f"{res}/{{sample}}/input/long_fastqc.data.txt"
    shell:
        "fastqc --noextract -t {threads} "
        f"-o {res}/{{wildcards.sample}}/input {{input}} ; "
        f"unzip -p {res}/{{wildcards.sample}}/input/long_fastqc.zip "
        f"long_fastqc/fastqc_data.txt > {{output.data}}"



# Kraken2 contamination on the merged raw paired input
rule kraken2_contamination_raw_paired:
    threads: 1
    conda:
        "envs/kraken2.yaml"
    input:
        R1 = f"{res}/{{sample}}/input/R1.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.fastq.gz"
    output:
        report = f"{res}/{{sample}}/input/kraken.raw.paired",
    params:
        kdb = config['kraken_db']
    log:
        f"{res}/{{sample}}/input/kraken.raw.paired.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed --paired {input.R1} {input.R2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# Detection of post-trim contamination in single reads that will be used for
# mapping
# Kraken2 contamination on the merged raw paired input
rule kraken2_contamination_trim_paired:
    threads: 1
    conda:
        "envs/kraken2.yaml"
    input:
        R1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz"
    output:
        report = f"{res}/{{sample}}/input/kraken.trimmed.paired",
    params:
        kdb = config['kraken_db']
    log:
        f"{res}/{{sample}}/input/kraken.trimmed.paired.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed --paired {input.R1} {input.R2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

rule kraken2_contamination_single:
    threads: 1
    conda:
        "envs/kraken2.yaml"
    input:
        S1 = f"{res}/{{sample}}/input/S1.trim.fastq.gz",
        S2 = f"{res}/{{sample}}/input/S2.trim.fastq.gz"
    output:
        report = f"{res}/{{sample}}/input/kraken.trimmed.single",
    params:
        kdb = config['kraken_db']
    log:
        f"{res}/{{sample}}/input/kraken.trimmed.single.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "--gzip-compressed {input.S1} {input.S2} "
        "--report {output.report} --output - 2>&1 | tee {log}"

# Produces a basic assembly using raw reads for MRCA and species typing
rule shovill_assembly:
    threads: 4
    conda:
        "envs/shovill.yaml"
    input:
        R1 = f"{res}/{{sample}}/input/R1.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.fastq.gz"
    output:
        f"{res}/{{sample}}/shovill_assembly/{{sample}}.shovill.contigs.fa",
    params:
        shovill_tmp = config['shovill_tempdir'],
        shovill_mem = config['shovill_mem']
    log:
        f"{res}/{{sample}}/shovill_assembly/{{sample}}.shovill.log"
    shell:
        # Note: letting Shovill set memory caps at 8, have gotten strange out
        # of memory errors in that pipeline see:
        # <https://github.com/tseemann/shovill/issues/59>
        "shovill --force --cpus {threads} --R1 {input.R1} --R2 {input.R2} "
        f"--outdir {res}/{{wildcards.sample}}/shovill_assembly "
        "--tmpdir {params.shovill_tmp} "
        "--ram {params.shovill_mem} "
        " 2>&1 | tee {log} ; "
        f"cp {res}/{{wildcards.sample}}/shovill_assembly/contigs.fa {{output}}"

rule dflye_assembly:
    threads: 8
    conda:
        "envs/dflye.yaml"
    params:
        size = config['dflye_genome_size'],
        medaka = config['_DFLYE_MEDAKA'],
        racon = config['dflye_racon']
    input:
        f"{res}/{{sample}}/input/long.fastq.gz"
    output:
        f"{res}/{{sample}}/dflye/contigs.fa"
    log:
        f"{res}/{{sample}}/dflye/dflye_assembly.log"
    shell:
        "dragonflye --reads {input} --trim --gsize {params.size} --outdir "
        f"{res}/{{wildcards.sample}}/dflye --force "
        "--cpus {threads} {params.medaka} --racon {params.racon} "
        "2>&1 | tee {log}"

rule dflye_short_polish_assembly:
    threads: 8
    conda:
        "envs/dflye.yaml"
    params:
        size = config['dflye_genome_size'],
        medaka = config['_DFLYE_MEDAKA'],
        racon = config['dflye_racon']
    input:
        L = f"{res}/{{sample}}/input/long.fastq.gz",
        R1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz"
    output:
        f"{res}/{{sample}}/dflye_short_polish/contigs.fa"
    log:
        f"{res}/{{sample}}/dflye_short_polish/dflye_assembly.log"
    resources:
        mem_mb=62000
    shell:
        "dragonflye --reads {input.L} --R1 {input.R1} --R2 {input.R2} "
        "--trim --gsize {params.size} --outdir "
        f"{res}/{{wildcards.sample}}/dflye_short_polish "
        "--cpus {threads} --racon {params.racon} {params.medaka} "
        # "--pilon 1 " # ignore for now, due to memory
        "--polypolish 3 --force "
        " 2>&1 | tee {log}"


# Kraken2 contamination of assembly
rule kraken2_contamination_assembly:
    threads: 1
    conda:
        "envs/kraken2.yaml"
    input:
        f"{res}/{{sample}}/shovill_assembly/{{sample}}.shovill.contigs.fa"
    output:
        report = f"{res}/{{sample}}/shovill_assembly/kraken.assembly"
    params:
        kdb = config["kraken_db"]
    log:
        f"{res}/{{sample}}/shovill_assembly/kraken.assembly.log"
    shell:
        "kraken2 --db {params.kdb} --threads {threads} --quick "
        "{input} --report {output.report} --output - 2>&1 | tee {log}"

# Basic BLAST detection of Erm(41) - looking for intact or truncated CDS
rule assembly_status_erm41:
    threads: 1
    conda:
        "envs/blast.yaml"
    input:
        f"{res}/{{sample}}/shovill_assembly/{{sample}}.shovill.contigs.fa"
    output:
        f"{res}/{{sample}}/{{sample}}.erm41.status"
    params:
        query_file = config["erm41_query"]
    shell:
        "tblastn -subject {input} -query {params.query_file} -outfmt 6 "
        "| head -n 1 | cut -f 2,3,4,6,7,8,9,10 > {output}"

# Everything gets mapped to ATCC19977 to determine basic coverage for QC
# Consider splitting this up to use threads better
rule map_to_mabs:
    conda:
        "envs/bwa.yaml"
    threads: 4
    input:
        R1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz",
        S1 = f"{res}/{{sample}}/input/S1.trim.fastq.gz",
        idx = f"workflow/resources/alignment_references/mabscessus.fasta.amb"
    params:
        execdir = exec_dir
    output:
        paired_temp = temp(f"{res}/{{sample}}/ref_mapping/paired.sorted.bam"),
        paired_temp_bai = temp(f"{res}/{{sample}}/ref_mapping/paired.sorted.bam.bai"),
        single_temp = temp(f"{res}/{{sample}}/ref_mapping/single.sorted.bam"),
        single_temp_bai = temp(f"{res}/{{sample}}/ref_mapping/single.sorted.bam.bai"),
        merge_temp = temp(f"{res}/{{sample}}/ref_mapping/merge.bam"),
        merge_sorted = f"{res}/{{sample}}/ref_mapping/mabs/merged.sorted.bam",
        merge_sorted_bai = f"{res}/{{sample}}/ref_mapping/mabs/merged.sorted.bam.bai"
    shell:
        "bwa-mem2 mem -t {threads} -M "
        "workflow/resources/alignment_references/mabscessus.fasta "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort - > "
        "{output.paired_temp} ;"
        "bwa-mem2 mem -t {threads} -M "
        "workflow/resources/alignment_references/mabscessus.fasta "
        "{input.S1} | samtools view -Sbh - | samtools sort -@ {threads} - > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort -@ {threads} - > "
        "{output.merge_sorted} ;"
        "samtools index {output.merge_sorted} {output.merge_sorted_bai}"

rule depth_map_to_mabs:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/ref_mapping/mabs/merged.sorted.bam"
    params:
        execdir = exec_dir
    output:
        depth = f"{res}/{{sample}}/ref_mapping/mabs/merged.sorted.depth.gz"
    shell:
        "bedtools genomecov -d -ibam {input} -g "
        "workflow/resources/alignment_references/mabcessus.fasta | gzip > "
        "{output.depth} "

# This rule puts everything in context of a rough mashtree to determine MRCA
# for downstream rules - MLST, specific references
rule mashtree_assemblies:
    threads: 8
    conda:
        "envs/mashtree.yaml"
    input:
        # this is ugly - but there are fewer constraints on input
        # wildcards in input must match wildcards in output, but since
        # this output rule has no wildcards, needed to use something else
        # internal = glob.glob("shovill_assembly/*.shovill.contigs.fa"),
        # internal2 = glob.glob("../bakMA/4-shovill_assemble/*.contigs.fa"),
        # glob.glob("external/refseq/complete/*.gbff.gz")
        # aggregate_refseq is an input function that triggers the download
        # of complete refseq genomes
        aggregate_refseq,
        internal = expand(
            "{res}/{sample}/shovill_assembly/{sample}.shovill.contigs.fa",
            res=res,
            sample=sample_names
        )
    params:
        # the mashtree tempdir can get large, depending on the number of
        # species - see config.yaml
        tempdir = config['mashtree_tempdir']
    output:
        tree = f"{res}/mashtree/assembly_mashtree.{{nodeset}}.tree",
        matrix = f"{res}/mashtree/assembly_mashtree.{{nodeset}}.matrix"
    log:
        f"{res}/mashtree/assembly_mashtree.{{nodeset}}.log"
    shell:
        "mashtree --mindepth 0 --numcpus {threads} --outtree {output.tree} "
        "--sort-order random "
        "--outmatrix {output.matrix} --tempdir {params.tempdir} "
        "{input} 2>&1 | tee {log}"

# Note:
# see https://www.frontiersin.org/articles/10.3389/fcimb.2022.816615/full for
# a multiplex PCR scheme to differentiate mabs subspecies
rule MRCA_from_tree:
    threads: 1
    conda:
        "envs/ngd_phylo.yaml"
    input:
        f"{res}/mashtree/assembly_mashtree.complete.tree"
    output:
        f"{res}/{{sample}}/{{sample}}.MRCA.csv"
    shell:
        "workflow/scripts/tree_MRCA.py ref {input} {wildcards.sample} > {output}"


# This rule determines which samples will use which reference or MLST scheme
# note that in October 2020 the scheme from pubMLST.org was updated to include
# the clonal complex, and mmassiliense was removed. Versions of mlst beyond
# 2.19.x now use the updated mabscesus scheme
rule MRCA_MLST:
    threads: 1
    conda:
        "envs/ngd_phylo.yaml"
    input:
        tree = f"{res}/mashtree/assembly_mashtree.complete.tree",
        assembly = f"{res}/{{sample}}/shovill_assembly/{{sample}}.shovill.contigs.fa"
    output:
        f"{res}/{{sample}}/{{sample}}.MLST.csv"
    shell:
        "mlst --scheme mabscessus --threads {threads} "
        "{input.assembly} > {output}"


def QC_input(wildcards):
    outf = {}    
    if wildcards.sample in sample_names:
        outf['short_preR1_fastqc'] = f"{res}/{wildcards.sample}/input/R1_fastqc.data.txt",
        outf['short_preR2_fastqc'] = f"{res}/{wildcards.sample}/input/R2_fastqc.data.txt",
        outf['short_tR1_fastqc'] = f"{res}/{wildcards.sample}/input/R1.trim_fastqc.data.txt",
        outf['short_tR2_fastqc'] = f"{res}/{wildcards.sample}/input/R2.trim_fastqc.data.txt",
        outf['short_tS1_fastqc'] = f"{res}/{wildcards.sample}/input/S1.trim_fastqc.data.txt",
        outf['short_tS2_fastqc'] = f"{res}/{wildcards.sample}/input/S2.trim_fastqc.data.txt",
        outf['tree'] = f"{res}/mashtree/assembly_mashtree.complete.tree",
        outf['raw_kraken'] = f"{res}/{wildcards.sample}/input/kraken.raw.paired",
        outf['trim'] = f"{res}/{wildcards.sample}/input/trimmomatic.log",
        outf['pair_kraken2'] = f"{res}/{wildcards.sample}/input/kraken.trimmed.paired",
        outf['single_kraken2'] = f"{res}/{wildcards.sample}/input/kraken.trimmed.single",
        outf['short_assembly'] = f"{res}/{wildcards.sample}/shovill_assembly/{wildcards.sample}.shovill.contigs.fa",
        outf['MRCA'] = f"{res}/{wildcards.sample}/{wildcards.sample}.MRCA.csv",
        outf['mlst'] = f"{res}/{wildcards.sample}/{wildcards.sample}.MLST.csv",
        outf['erm41'] = f"{res}/{wildcards.sample}/{wildcards.sample}.erm41.status",
        outf['mabs_bam'] = f"{res}/{wildcards.sample}/ref_mapping/mabs/merged.sorted.bam",
        outf['mabs_depth'] = f"{res}/{wildcards.sample}/ref_mapping/mabs/merged.sorted.depth.gz",
    if wildcards.sample in long_sample_names:
        outf['long_fastqc'] = f"{res}/{wildcards.sample}/input/long_fastqc.data.txt",
        outf['long_assembly'] = f"{res}/{wildcards.sample}/dflye/contigs.fa",
    if wildcards.sample in both_samples:
        outf['long_polish_assembly'] = f"{res}/{wildcards.sample}/dflye_short_polish/contigs.fa",

    return outf


rule QC_stats_per_sample:
    conda: "envs/phy_plots.yaml"
    params:
        resdir = res
    input:
        unpack(lambda wildcards: QC_input(wildcards))
    output:
        f"{res}/{{sample}}/{{sample}}.QC.csv"
    params:
        execdir = exec_dir
    log:
        f"{res}/{{sample}}/sample.QC.log"
    shell:
        "workflow/scripts/sample_stats.py {wildcards.sample} {params.resdir} "
        " > {output} 2> {log}"

rule prokka_annotate_dflye:
    conda:
        "envs/prokka.yaml"
    threads: 8
    params:
        ref = lambda wildcards: ref_from_QC(
            wildcards.sample, f"{res}/QC_summary.csv"
        )
    input:
        f"{res}/{{sample}}/{{dflye}}/contigs.fa",
        f"{res}/QC_summary.csv"
    output:
        f"{res}/{{sample}}/{{sample}}.{{dflye}}.prokka.gff"
    shell:
        "prokka --prefix {wildcards.sample}.{wildcards.dflye}.prokka "
        f"--outdir {res}/{{wildcards.sample}}/ --force --kingdom Bacteria "
        "--compliant --cpus {threads}  "
        "--proteins workflow/resources/alignment_references/{params.ref}.gbk "
        "{input}"


rule merge_QC_summary:
    priority: 10
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        expand(
            "{res}/{sample}/{sample}.QC.csv",
            res=res,
            sample=sample_names
        )
    output:
        f"{res}/QC_summary.csv"
    shell:
        "workflow/scripts/merge_QC_csv.py {output} {input}"

