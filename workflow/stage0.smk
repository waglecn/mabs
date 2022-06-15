#########################################################################
# External Data Rules
#########################################################################
checkpoint refseq_outgroup_download:
    conda:
        'envs/ngd_phylo.yaml'
    params:
        outgroup = config['outgroup_assembly'],
        execdir = exec_dir
    output:
        directory("workflow/resources/external/refseq/outgroup/")
    log:
        "workflow/resources/external/refseq_outgroup_download.log"
    shell:
        "ngd -s refseq -r3 --flat-output -A {params.outgroup} -F genbank "
        "-o {output} bacteria ; "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"

rule download_refseq_taxlist:
    conda:
        "envs/ngd_phylo.yaml"
    params:
        taxid = config['target_species_taxid'],
        execdir = exec_dir
    output:
        "workflow/resources/external/taxid{taxid}.target.taxlist".format(
            taxid=config['target_species_taxid']
        )
    shell:
        "python workflow/scripts/gimme-taxa.py -j "
        "{params.taxid} > {output}"

checkpoint refseq_complete_genome_download:
    conda:
        'envs/ngd_phylo.yaml'
    params:
        taxid = config['target_species_taxid'],
        execdir = exec_dir
    input:
        "workflow/resources/external/taxid{taxid}.target.taxlist".format(
            taxid=config['target_species_taxid']
        )
    log:
        "workflow/resources/external/refseq_complete_download.log"
    output:
        directory("workflow/resources/external/refseq/complete/"),
        # a non-dynamic output may not exist with dynamic output
    shell:
        "ngd -v -s refseq -r 3 -p 3 -o {output} "
        "-t {input} --flat-output -F genbank "
        "-m workflow/resources/external/complete.refseq.meta "
        "-l complete bacteria ; "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"

checkpoint refseq_all_genome_download:
    conda:
        "envs/ngd_phylo.yaml"
    params:
        execdir = exec_dir,
        taxid = config['target_species_taxid']
    input:
        "workflow/resources/external/taxid{taxid}.target.taxlist".format(
            taxid=config['target_species_taxid']
        )
    log:
        "workflow/resources/external/refseq_all_download.log"
    output:
        directory("workflow/resources/external/refseq/all/")
    shell:
        "ngd -s refseq -r 3 -p 3 -o {output} "
        "-t {input} --flat-output -F genbank "
        "-m workflow/resources/external/all.refseq.meta bacteria "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        "done 2>&1 | tee {log}"

checkpoint refseq_download_landmarks:
    conda:
        "envs/ngd_phylo.yaml"
    params:
        execdir = exec_dir,
        accessions = ",".join(config["mash_ref_taxa"].values())
    log:
        "workflow/resources/external/refseq_landmark_download.log"
    output:
        directory("workflow/resources/external/refseq/landmarks/")
    shell:
        "ngd -r 3 -p 3 -o {output} --flat-output -F genbank "
        "-A {params.accessions} bacteria ; "
        "for i in {output}/*.gbff.gz ; do "
        r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
        "echo rename $i to $new ; mv $i $new ; "
        " done 2>&1 | tee {log}"


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
    checkpoint_landmark = checkpoints.refseq_download_landmarks.get(
        **wildcards
    ).output[0]
    assembly_name = glob.glob(checkpoint_set + "/*.gbk.gz")
    outgroup_name = glob.glob(checkpoint_outgroup + "/*.gbk.gz")
    landmark_name = glob.glob(checkpoint_landmark + "/*.gbk.gz")

    assembly_files = assembly_name + outgroup_name + landmark_name
    return assembly_files

rule get_gatk_jar:
    conda:
        "envs/bzip2.yaml"
    params:
        jar = config['gatk3_jar'],
        url = config['gatk3_jar_url']
    output:
        "workflow/resources/{params.jar}"
    shell:
        "wget {params.url} -P workflow/resources/ ; "
        "bunzip2 -c workflow/resources/GenomeAnalysisTK-*.tar.bz2 > {output}"


rule register_gatk:
    conda: "envs/gatk3.yaml"
    output:
        "workflow/resources/gatk-registered"
    params:
        gatk_jar = "workflow/resources/{gatk_jar}".format(
            gatk_jar=config['gatk3_jar'],
        )
    shell:
        "echo {output} ; gatk3-register {params.gatk_jar} ; touch {output}"


rule prep_references:
    input:
        expand("workflow/resources/alignment_references/{ref}.fasta", ref=[
            'mabscessus', 'mmassiliense', 'mbolletii'
        ])


rule get_alignment_references:
    conda:
        "envs/ngd_phylo.yaml"
    threads: 1
    params:
        acc = lambda wildcards: config['mash_ref_taxa'][wildcards.ref]
    output:
        "workflow/resources/alignment_references/{ref,\w+}.fasta"
    shell:
        "ngd -s refseq -r3 --flat-output -A {params.acc} -F fasta "
        "-o workflow/resources/alignment_references bacteria ; "
        "gunzip -c workflow/resources/alignment_references/{params.acc}*.gz > "
        "{output} "

rule make_mapping_bwa_index:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        "workflow/resources/alignment_references/{ref}.fasta"
    output:
        "workflow/resources/alignment_references/{ref}.fasta.amb"
    shell:
        "bwa-mem2 index {input}"

