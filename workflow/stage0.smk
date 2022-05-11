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
        taxid = config['target_species_taxid'],
        execdir = exec_dir
    output:
        "workflow/resources/external/taxid{taxid}.target.taxlist".format(
            taxid=config['target_species_taxid']
        )
    shell:
        "python {params.execdir}/scripts/gimme-taxa.py -j "
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
        "-m external/complete.refseq.meta "
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
        "external/taxid{taxid}.target.taxlist".format(
            taxid=config['target_species_taxid']
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

checkpoint refseq_download_landmarks:
    conda:
        "envs/ngd_phylo.yaml"
    params:
        execdir = exec_dir,
        accessions = ",".join(config["mash_ref_taxa"])
    log:
        "external/refseq_landmark_download.log"
    output:
        directory("external/refseq/landmarks/")
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


rule register_gatk:
    conda: "envs/gatk3.yaml"
    output: os.path.join(exec_dir, 'resources/gatk-registered')
    params:
        gatk_jar = "resources/{gatk_jar}".format(
            gatk_jar=config['gatk3_jar'],
            execdir=exec_dir
        )
    shell:
        "echo {output} ; gatk3-register {params.gatk_jar} ; touch {output}"

