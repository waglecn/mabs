# need to handle both long and short reads in stage2
import yaml
try:
    stage2 = yaml.safe_load(open(f"{res}/stage2test.yaml"))
    print(stage2)
except FileNotFoundError as e:
    print("stage2 not configured yet", file=sys.stderr)
    stage2 = {}

include: "stage2_short.smk"
include: "stage2_long.smk"

rule copy_short_internal_reference:
    threads: 1
    conda: "envs/bwa.yaml"
    input:
        f"{res}/samples/{{s}}/shovill_assembly/{{s}}.shovill.contigs.fa"
    output:
        "workflow/resources/alignment_references/{s}.short.fasta",
        "workflow/resources/alignment_references/{s}.short.gbk"
    shell:
        "cp {input} {output[0]} ; touch {output[1]}"


rule copy_long_internal_reference:
    threads: 1
    conda: "envs/bwa.yaml"
    input:
        f"{res}/samples/{{s}}/dflye/{{s}}.dflye.contigs.fa"
    output:
        "workflow/resources/alignment_references/{s}.long.fasta",
        "workflow/resources/alignment_references/{s}.long.gbk"
    shell:
        "cp {input} {output[0]} ; touch {output[1]}"


rule copy_long_polished_internal_reference:
    threads: 1
    conda: "envs/bwa.yaml"
    input:
        f"{res}/samples/{{s}}/dflye_short_polish/{{s}}.polish.contigs.fa"
    output:
        "workflow/resources/alignment_references/{s}.polish.fasta",
        "workflow/resources/alignment_references/{s}.polish.gbk"
    shell:
        "cp {input} {output[0]} ; touch {output[1]}"


rule stage2_short_all_outputs:
    input:
        ##################################
        # STAGE 2 SHORT - type-specific analyses
        ##################################
        [expand(
            [
                "{res}/samples/{s}/map/{ref}/short.mpileup",
                "{res}/samples/{s}/map/{ref}/short.filter.vcf.gz",
                "{res}/samples/{s}/map/{ref}/short.filter.failed.vcf.gz",
                "{res}/samples/{s}/map/{ref}/short.mask.bed",
                "{res}/samples/{s}/map/{ref}/short.variants.bed",
                "{res}/samples/{s}/map/{ref}/short.0cov.bed",
                "{res}/samples/{s}/map/{ref}/short.filter.hvar.vcf.gz",
                "{res}/samples/{s}/map/{ref}/short.filter.AD_failed.vcf.gz",
                "{res}/samples/{s}/map/{ref}/short.lowcov.bed",
                "{res}/samples/{s}/map/{ref}/short.consensus.fa",
                "{res}/samples/{s}/map/{ref}/short.filter.hvar_DF.bed",
            ],
            res=res,
            s=[i for i in stage2[ref] if i in filt_short_samples],
            ref=ref
            #     sample=s
            ) for ref in stage2        
        ]

rule stage2_long_all_outputs:
    input:
        ##################################
        # STAGE 2 LONG - type-specific analyses
        ##################################
        [expand(
            [
                "{res}/samples/{s}/map/{ref}/long.mpileup",
                "{res}/samples/{s}/map/{ref}/long.filter.vcf.gz",
                "{res}/samples/{s}/map/{ref}/long.filter.failed.vcf.gz",
                "{res}/samples/{s}/map/{ref}/long.filter.AD_failed.vcf.gz",
                "{res}/samples/{s}/map/{ref}/long.filter.hvar.vcf.gz",
                "{res}/samples/{s}/map/{ref}/long.0cov.bed",
                "{res}/samples/{s}/map/{ref}/long.lowcov.bed",
                "{res}/samples/{s}/map/{ref}/long.variants.bed",
                "{res}/samples/{s}/map/{ref}/long.mask.bed",
                "{res}/samples/{s}/map/{ref}/long.consensus.fa",
                "{res}/samples/{s}/map/{ref}/long.filter.hvar_DF.bed",
            ],
            res=res,
            s=[i for i in stage2[ref] if i in filt_long_samples],
            ref=ref
            ) for ref in stage2
         ]

rule stage2:
    input:
        rules.stage2_short_all_outputs.input,
        rules.stage2_long_all_outputs.input,

