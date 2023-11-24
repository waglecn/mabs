import yaml
try:
    stage2 = yaml.safe_load(open(f"{res}/stage2test.yaml"))
    print(stage2)
except FileNotFoundError as e:
    print("stage2 not configured yet", file=sys.stderr)
    stage2 = {}

fos = [s for s in filt_short_samples if s in short_only_sample_names]
fol = [s for s in filt_long_samples if s in long_only_sample_names]
bo = [
    s for s in both_samples if 
    (s not in config['exclude_short_samples']) and 
    (s not in config['exclude_long_samples'])
]
slb = fos + fol + bo
######################################
# Stage 3
######################################

rule stage3:
    input:
        ##########################
        # ref only
        ##########################
        expand(
            [
                "{res}/gubbins/{ref}.fasta", 
                "{res}/gubbins/{ref}.gubbins.done",
                "{res}/gubbins/{ref}.gubbins.bed",
                "{res}/map/{ref}.consensus.fasta",
                "{res}/map/{ref}.consensus.fasta.treefile",
                "{res}/map/{ref}.consensus.fasta.snpdists.csv"
            ], 
            res=res, ref=list(stage2.keys())
        ),
        #########################
        # sample by ref (short)
        #########################
        [
            expand(
                [
                    "{res}/samples/{s}/map/{ref}/short_final.bed",
                    "{res}/samples/{s}/map/{ref}/short.combined.consensus.fasta",
                ],
                res=res,
                s=[i for i in stage2[sref] if i in filt_short_samples],
                ref=sref
            ) for sref in list(stage2.keys())
        ],
        #########################
        # sample by ref (long)
        #########################
        [
            expand(
                [
                    "{res}/samples/{s}/map/{ref}/long_final.bed",
                    "{res}/samples/{s}/map/{ref}/long.combined.consensus.fasta",
                ],
                res=res,
                s=[i for i in stage2[sref] if i in filt_long_samples],
                ref=sref
            ) for sref in list(stage2.keys())
        ],
        


rule merge_consensus_fasta:
    threads: 1
    input:
        lambda wildcards: expand(
            "{res}/samples/{s}/map/{ref}/short.consensus.fa",
            res=res,
            s=[i for i in stage2[wildcards.ref] if i in filt_short_samples],
            ref=wildcards.ref 
        ) + expand(
                "{res}/samples/{s}/map/{ref}/long.consensus.fa",
                res=res,
                s=[i for i in stage2[wildcards.ref] if i in fol],
                ref=wildcards.ref
        )
    output:
        f"{res}/gubbins/{{ref}}.fasta"
    shell:
        "for i in {input} ; do cat $i >> {output} ; done"

rule run_gubbins:
    conda:
        "envs/gubbins.yaml"
    threads: 32
    params:
        prefix = f"{res}/gubbins/{{ref}}"
    input:
        f"{res}/gubbins/{{ref}}.fasta"
    output:
        f"{res}/gubbins/{{ref}}.gubbins.done"
    log:
        f"{res}/gubbins/{{ref}}.gubbins.log"
    shell:
        "run_gubbins.py --prefix {params.prefix} "
        "--threads {threads} {input} || true  && "
        "touch {output} ; rm -f {wildcards.ref}.* "

rule make_gubbins_bed:
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    params:
        prefix = f"{res}/gubbins/{{ref}}"
    input:
        sentinel = f"{res}/gubbins/{{ref}}.gubbins.done"
    output:
        f"{res}/gubbins/{{ref}}.gubbins.bed"
    shell:
        "workflow/scripts/make_bed_from_gubbins.py "
        "{params.prefix}.recombination_predictions.embl "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.sentinel} | bedtools sort -i stdin | bedtools merge -i stdin "
        "-c 4 -o distinct > {output}"

rule merge_per_sample_short_beds:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        gubbinsbed = f"{res}/gubbins/{{ref}}.gubbins.bed",
        PEbed = f"workflow/resources/alignment_references/{{ref}}.PE_PPE.bed",
        maskbed = f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}.mask.bed",
        DFbed = f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}.filter.hvar_DF.bed",
    output:
        f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}_final.bed"
    shell:
        "cat {input.maskbed} {input.DFbed} {input.gubbinsbed} {input.PEbed} | "
        "bedtools sort -i stdin | bedtools merge -i stdin -c 4 "
        "-o distinct_only > {output}"

rule make_short_masked_output:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        ref = "workflow/resources/alignment_references/{ref}.fasta",
        bed = f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}_final.bed",
        vcf = f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}.filter.hvar.vcf.gz"
    output:
        f"{res}/samples/{{s}}/map/{{ref}}/{{long_short}}.combined.consensus.fasta"
    shell:
        "bcftools consensus -p {wildcards.s}.{wildcards.long_short} "
        "-f {input.ref} --mark-del '-' "
        "-m {input.bed} -i 'strlen(REF)>=strlen(ALT) & INFO/MQ >= 20 & FORMAT/DP >= 10' "
        "{input.vcf} | sed \"/^/s/{wildcards.s}.*/{wildcards.s}.{wildcards.long_short}/\" > {output} "

rule concatenate_consensus:
    threads:1
    conda:
        "envs/bwa.yaml"
    input:
        ref = "workflow/resources/alignment_references/{ref}.fasta",
        isolates = lambda wildcards: [
            f"{res}/samples/{i}/map/{wildcards.ref}/short.combined.consensus.fasta" for i in stage2[wildcards.ref] if i in fos
        ] + [
            f"{res}/samples/{i}/map/{wildcards.ref}/long.combined.consensus.fasta" for i in stage2[wildcards.ref] if i in fol
        ]
    output:
        f"{res}/map/{{ref}}.consensus.fasta"
    shell:
        "for i in {input.ref} {input.isolates} ; do cat $i >> {output} ; done"

rule SNPs_phylogeny:
    conda:
        "envs/iqtree.yaml"
    threads: 32
    input:
        f"{res}/map/{{ref}}.consensus.fasta"
    output:
        f"{res}/map/{{ref}}.consensus.fasta.treefile"
    shell:
        "iqtree -s {input} -nt 32 -m MFP -B 10000 -redo --nmax 2000"

rule SNPs_matrix:
    threads: 1
    conda:
        "envs/phy_plots.yaml"
    input:
        f"{res}/map/{{ref}}.consensus.fasta"
    output:
        f"{res}/map/{{ref}}.consensus.fasta.snpdists.csv"
    shell:
        "snp-dists -bc {input} > {output}"

