rule prep_references:
    input:
        expand("workflow/resources/alignment_references/{ref}.{ext}",
            ref=['mabscessus', 'mmassiliense', 'mbolletii'],
            ext=['fasta', 'gbk']
        )


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

# rule download_refseq_taxlist:
#     conda:
#         "envs/ngd_phylo.yaml"
#     params:
#         taxid = config['target_species_taxid'],
#         execdir = exec_dir
#     output:
#         "workflow/resources/external/taxid{taxid}.target.taxlist".format(
#             taxid=config['target_species_taxid']
#         )
#     shell:
#         "python workflow/scripts/gimme-taxa.py -j "
#         "{params.taxid} > {output}"
#
# checkpoint refseq_complete_genome_download:
#     conda:
#         'envs/ngd_phylo.yaml'
#     params:
#         taxid = config['target_species_taxid'],
#         execdir = exec_dir
#     input:
#         "workflow/resources/external/taxid{taxid}.target.taxlist".format(
#             taxid=config['target_species_taxid']
#         )
#     log:
#         "workflow/resources/external/refseq_complete_download.log"
#     output:
#         directory("workflow/resources/external/refseq/complete/"),
#         # a non-dynamic output may not exist with dynamic output
#     shell:
#         "ngd -v -s refseq -r 3 -p 3 -o {output} "
#         "-t {input} --flat-output -F genbank "
#         "-m workflow/resources/external/complete.refseq.meta "
#         "-l complete bacteria ; "
#         "for i in {output}/*.gbff.gz ; do "
#         r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
#         "echo rename $i to $new ; mv $i $new ; "
#         "done 2>&1 | tee {log}"
#
# checkpoint refseq_all_genome_download:
#     conda:
#         "envs/ngd_phylo.yaml"
#     params:
#         execdir = exec_dir,
#         taxid = config['target_species_taxid']
#     input:
#         "workflow/resources/external/taxid{taxid}.target.taxlist".format(
#             taxid=config['target_species_taxid']
#         )
#     log:
#         "workflow/resources/external/refseq_all_download.log"
#     output:
#         directory("workflow/resources/external/refseq/all/")
#     shell:
#         "ngd -s refseq -r 3 -p 3 -o {output} "
#         "-t {input} --flat-output -F genbank "
#         "-m workflow/resources/external/all.refseq.meta bacteria "
#         "for i in {output}/*.gbff.gz ; do "
#         r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
#         "echo rename $i to $new ; mv $i $new ; "
#         "done 2>&1 | tee {log}"
#
# checkpoint refseq_download_landmarks:
#     conda:
#         "envs/ngd_phylo.yaml"
#     params:
#         execdir = exec_dir,
#         accessions = ",".join(config["landmark_taxa"].values())
#     log:
#         "workflow/resources/external/refseq_landmark_download.log"
#     output:
#         directory("workflow/resources/external/refseq/landmarks/")
#     shell:
#         "ngd -r 3 -p 3 -o {output} --flat-output -F genbank "
#         "-A {params.accessions} bacteria ; "
#         "for i in {output}/*.gbff.gz ; do "
#         r"new=$(echo $i | sed -e s/\.gbff\.gz/\.gbk\.gz/) ; "
#         "echo rename $i to $new ; mv $i $new ; "
#         " done 2>&1 | tee {log}"


def aggregate_refseq(wildcards):
    # if wildcards.nodeset == 'complete':
    #     checkpoint_set = checkpoints.refseq_complete_genome_download.get(
    #         **wildcards
    #     ).output[0]
    # elif wildcards.nodeset == 'all':
    #     checkpoint_set = checkpoints.refseq_all_genome_download.get(
    #         **wildcards
    #     ).output[0]
    # landmark_name = []
    # if len(config['landmark_taxa']) > 0:
    #     checkpoint_landmark = checkpoints.refseq_download_landmarks.get(
    #         **wildcards
    #     ).output[0]
    #     landmark_name = glob.glob(checkpoint_landmark + "/*.gbk.gz")
    checkpoint_outgroup = checkpoints.refseq_outgroup_download.get(
        **wildcards
    ).output[0]
    # assembly_name = glob.glob(checkpoint_set + "/*.gbk.gz")
    outgroup_name = glob.glob(checkpoint_outgroup + "/*.gbk.gz")
    ref_name = expand(
        "workflow/resources/alignment_references/{ref}.{ext}",
        ref=['mabscessus', 'mmassiliense', 'mbolletii'],
        ext=['fasta', 'gbk']
    )
    return outgroup_name + ref_name


rule get_alignment_references:
    conda:
        "envs/ngd_phylo.yaml"
    threads: 1
    params:
        acc = lambda wildcards: config['ref_taxa'][wildcards.ref]
    output:
        fn = "workflow/resources/alignment_references/{ref,\w+}.fasta",
        gbk = "workflow/resources/alignment_references/{ref,\w+}.gbk"
    shell:
        "ngd -s refseq -r3 --flat-output -A {params.acc} -F genbank,fasta "
        "-o workflow/resources/alignment_references bacteria ; "
        "gunzip -f workflow/resources/alignment_references/{params.acc}*.gz ;"
        "mv workflow/resources/alignment_references/{params.acc}*.fna {output.fn} ;"
        "mv workflow/resources/alignment_references/{params.acc}*.gbff {output.gbk} ;"

rule make_mapping_bwa_index:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        # "workflow/resources/alignment_references/{ref}.fasta"
        "{ref}.fasta"
    output:
        "{ref}.fasta.amb"
        #"workflow/resources/alignment_references/{ref}.fasta.amb"
    shell:
        "bwa-mem2 index {input}"

# rule mash_references:
#     conda:
#         "envs/mash.yaml"
#     threads: 32
#     input:
#         expand("")

###############################
# GATK stuff
###############################
rule get_gatk_jar:
    conda:
        "envs/bzip2.yaml"
    params:
        url = config['gatk3_jar_url']
    output:
        f"workflow/resources/{config['gatk3_jar']}"
    shell:
        "wget {params.url} -P workflow/resources/ ; "
        "bunzip2 -c workflow/resources/GenomeAnalysisTK-*.tar.bz2 | "
        "tar Oxf - > {output}"


rule register_gatk:
    conda: "envs/gatk3.yaml"
    input:
        f"workflow/resources/{config['gatk3_jar']}"
    output:
        "workflow/resources/gatk-registered"
    shell:
        "echo {output} ; gatk3-register {input} ; touch {output}"
