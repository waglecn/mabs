###############################################################################
# STAGE 3
###############################################################################
rule short_merge_consensus_fasta:
    threads: 1
    input:
        in_fasta = lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{ref}/RG_SC.consensus.fa".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s
            ) for s in filt_short_samples if 
            ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ],
    output:
        f"{res}/gubbins/{{ref}}.RG_SC.short_only.fasta"
    shell:
        "cat {input.in_fasta} > {output}"

rule short_long_merge_consensus_fasta:
    threads: 1
    input:
        in_short_fasta = lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{ref}/RG_SC.consensus.fa".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s
            ) for s in filt_short_samples if
            ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ],
        in_long_fasta = lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{ref}/RG_SC.long.consensus.fa".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv", "long"),
                s=s
            ) for s in filt_long_samples if
            ref_from_QC(s, f"{res}/QC_summary.csv", "long") == wildcards.ref
        ],
    output:
        f"{res}/gubbins/{{ref}}.RG_SC.short_long.fasta"
    shell:
        "cat {input.in_short_fasta} {input.in_long_fasta} > {output}"

rule run_gubbins:
    conda:
        "envs/gubbins.yaml"
    threads: 8
    params:
        prefix = f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}"
    input:
        f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.fasta"
    output:
        # embl = f"{res}/gubbins/{{ref}}.{{step}}.short_only.recombination_predictions.embl",
        sentinel = f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.gubbins.done"
    shell:
        "run_gubbins.py --prefix {params.prefix} --min-snps 20 "
        "--threads {threads} {input} || true ; "
        "touch {output.sentinel} "

rule make_gubbins_bed:
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    params:
        prefix = f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}"
    input:
        # sentinel implies this exists
        # embl = f"{res}/gubbins/{{ref}}.{{kind}}.recombination_predictions.embl",
        sentinel = f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.gubbins.done"
    output:
        f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.gubbins.bed"
    shell:
        "workflow/scripts/make_bed_from_gubbins.py "
        "{params.prefix}.recombination_predictions.embl "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.sentinel} | bedtools sort -i stdin | bedtools merge -i stdin "
        "-c 4 -o distinct > {output}"

rule merge_per_sample_per_kind_beds:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        gubbinsbed = f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.gubbins.bed",
        PEbed = f"workflow/resources/alignment_references/{{ref}}.PE_PPE.bed",
        maskbed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.mask.bed",
        DFbed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.hvar_DF.bed",
    output:
        f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.combined.bed"
    shell:
        "cat {input.maskbed} {input.DFbed} {input.gubbinsbed} {input.PEbed} | "
        "bedtools sort -i stdin | bedtools merge -i stdin -c 4 "
        "-o distinct_only > {output}"

rule merge_per_sample_per_kind_long_beds:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        gubbinsbed =f"{res}/gubbins/{{ref}}.{{step}}.{{kind}}.gubbins.bed",
        PEbed = f"workflow/resources/alignment_references/{{ref}}.PE_PPE.bed",
        maskbed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.long.mask.bed",
        DFbed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.long.hvar_DF.bed"
    output:
        f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.long.combined.bed"
    shell:
        "cat {input.maskbed} {input.DFbed} {input.gubbinsbed} {input.PEbed} | "
        "bedtools sort -i stdin | bedtools merge -i stdin -c 4 "
        "-o distinct_only > {output}"

rule make_short_masked_output:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        bed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.combined.bed",
        vcf = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.hvar.vcf.gz",
        ref = f"workflow/resources/alignment_references/{{ref}}.fasta",
    output:
        f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.combined.consensus.fasta"
    shell:
        f"bcftools consensus -p {{wildcards.s}} "
        "-f {input.ref} --mark-del '-' "
        "-m {input.bed} -i 'strlen(REF)>=strlen(ALT) & INFO/MQ >= 20 & FORMAT/DP >= 10' "
        "{input.vcf} | sed \"/^>/s/{wildcards.s}.*/{wildcards.s}/\" > {output} "

rule make_long_masked_output:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        bed = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.long.combined.bed",
        vcf = f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.long.hvar.vcf.gz",
        ref = f"workflow/resources/alignment_references/{{ref}}.fasta",
    output:
        f"{res}/{{s}}/MRCA_ref_mapping/{{ref}}/{{step}}.{{kind}}.long.combined.consensus.fasta"
    shell:
        f"bcftools consensus -p {{wildcards.s}} "
        "-f {input.ref} --mark-del '-' "
        "-m {input.bed} -i 'strlen(REF)>=strlen(ALT) & INFO/MQ >= 20 & FORMAT/DP >= 10' "
        "{input.vcf} | sed \"/^>/s/{wildcards.s}.*/{wildcards.s}/\" > {output} "

rule concatenate_short_only_consensus:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        short = lambda wildcards: [
            (
                "{res}/{s}/MRCA_ref_mapping/{ref}/"
                "{{step}}.short_only.combined.consensus.fasta"
            ).format(
                res=res, s=s, ref=ref_from_QC(s, f"{res}/QC_summary.csv")
            ) for s in filt_short_samples if 
                ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ]
    output:
        f"{res}/SNP_phylo/{{ref}}.{{step}}.short_only.consensus.fasta"
    shell:
        "cat {input.short} > {output}"

rule concatenate_short_long_consensus:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        lambda wildcards: [
            (
                "{res}/{s}/MRCA_ref_mapping/{ref}/"
                "{{step}}.short_long.combined.consensus.fasta"
            ).format(
                res=res, s=s, ref=ref_from_QC(s, f"{res}/QC_summary.csv")
            ) for s in filt_short_samples if 
                ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ] + [
            (
                "{res}/{s}/MRCA_ref_mapping/{ref}/"
                "{{step}}.short_long.long.combined.consensus.fasta"
            ).format(
                res=res, s=s, ref=ref_from_QC(s, f"{res}/QC_summary.csv", "long")
            ) for s in filt_long_samples if 
                ref_from_QC(
                    s, f"{res}/QC_summary.csv", "long"
                ) == wildcards.ref
        ]
    output:
        f"{res}/SNP_phylo/{{ref}}.{{step}}.short_long.consensus.fasta",
    shell:
        "cat {input} > {output}"

rule MRCA_SNPs_phylogeny:
    conda:
        "envs/iqtree.yaml"
    threads: 8
    input:
        f"{res}/SNP_phylo/{{ref}}.{{step}}.{{kind}}.consensus.fasta"
    output:
        tree = f"{res}/SNP_phylo/{{ref}}.{{step}}.{{kind}}.consensus.fasta.treefile"
    shell:
        "iqtree -s {input} -nt {threads} -m MFP -B 10000 -redo --nmax 2000"

rule MRCA_SNPs_matrix:
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        f"{res}/SNP_phylo/{{ref}}.{{step}}.{{kind}}.consensus.fasta"
    output:
        f"{res}/SNP_phylo/{{ref}}.{{step}}.{{kind}}.consensus.fasta.snpdists.csv"
    shell:
        "snp-dists -bc {input} > {output}"
