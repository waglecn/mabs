###############################################################################
# STAGE 3
###############################################################################

# merge all zero-coverage positions for all samples
# def glob_0covbed(ref, samples, step):

rule merge_0cov_bed:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{ref}/RG_SC_RA.0cov.bed".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s)
            for s in filt_samples if
            ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ]
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.0cov.bed"
    shell:
        # the grep command will fail on the input of an empty bed file
        # something to do with bash strict mode?
        "echo '' | "
        "cat - {input} | "
        "grep 'NC_010394' -v | sort -k 1,1n -k 2,2n | uniq | "
        "bedtools merge -i - > {output}"

rule merge_DF_bed:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{ref}/RG_SC_RA_filter.hvar_DF.bed".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s)
            for s in filt_samples if
            ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ]
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.DF.bed"
    shell:
        "echo '' | cat - {input} | "
        "grep 'NC_010394' -v | sort -k 1,1n -k 2,2n | uniq | "
        "bedtools merge -i - > {output}"

rule merge_bed:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        zerocov = f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.0cov.bed",
        # this is probably covered by gubbins
        # DF = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.DF.bed",
        PEPPE = "workflow/resources/alignment_references/{ref}.PE_PPE.bed"
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.bed"
    shell:
        "echo '' | "
        "cat - {input} | "
        r"grep 'NC_010394\.1' -v | sort -k 1,1n -k 2,2n | cut -f 1-3 | "
        "bedtools merge -i - > {output}"

rule filter_vcf_with_bed:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        bed = f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.bed",
        vcf = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA_filter.hvar.vcf.gz",
        csi = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA_filter.hvar.vcf.gz.csi"
    output:
        vcf = f"{res}/{{sample}}/filtered_vcf/{{ref}}/RG_SC_RA_bedfilter.vcf.gz"
    shell:
        "bcftools view {input.vcf} -T ^{input.bed} -Oz -o {output.vcf} ;"

rule extract_sample_VCF_SNPs:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        vcf = f"{res}/{{sample}}/filtered_vcf/{{ref}}/RG_SC_RA_bedfilter.vcf.gz",
        csi = f"{res}/{{sample}}/filtered_vcf/{{ref}}/RG_SC_RA_bedfilter.vcf.gz.csi"
    output:
        f"{res}/{{sample}}/filtered_vcf/{{ref}}/RG_SC_RA_bedfilter.consensus.fasta",
    shell:
        "bcftools consensus -i 'type=\"SNP\"' "
        "-f workflow/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        " {input.vcf} > {output}; "

rule reheader_gubbins:
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        f"{res}/{{sample}}/filtered_vcf/{{ref}}/RG_SC_RA_bedfilter.consensus.fasta"
    output:
        f"{res}/{{sample}}/gubbins/{{ref}}/RG_SC_RA_bedfilter_gubbins.fasta"
    shell:
        "workflow/scripts/reheader_gubbins.py {input} {wildcards.sample} > {output}"

rule concatenate_gubbins:
    # 10215* 2/1/2021 17:31  cat ./mabscessus/*.fasta > mabs.fasta
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        samples = lambda wildcards: [
            "{res}/{s}/gubbins/{ref}/RG_SC_RA_bedfilter.gubbins.fasta".format(
                res=res,ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s
            ) for s in filt_samples if
            ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ],
    output:
        f"{res}/gubbins/{{ref}}.concatenated.fasta"
    shell:
        "cat {input} > {output} ; "
        "workflow/scripts/gubbins_ref.py {wildcards.ref} "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta >> "
        "{output}"

rule run_gubbins:
    # 10218* 2/1/2021 17:31  run_gubbins.py --prefix mabs --min_snps 20 --threads 8 mabs.fasta
    conda:
        "envs/gubbins.yaml"
    threads: 8
    params:
        prefix = f"{res}/gubbins/{{ref}}.concatenated"
    input:
        f"{res}/gubbins/{{ref}}.concatenated.fasta"
    output:
        # filtered_fasta = "gubbins/{ref}.concatenated."
        # "filtered_polymorphic_sites.fasta",
        embl = f"{res}/gubbins/{{ref}}.concatenated.recombination_predictions.embl",
        sentinel = f"{res}/gubbins/{{ref}}.gubbins.done"
    shell:
        "run_gubbins.py --prefix {params.prefix} --min-snps 20 "
        "--threads {threads} {input} || true ; "
        "touch {output.sentinel}"

rule make_gubbins_bed:
    # 10239* 2/1/2021 18:06  ../../analysis/make_bed_from_gubbins.py mabs.recombination_predictions.embl | sort -k 1,1n -k 2,2n > mabs.gubbins.bed
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        # sentinel implies this exists
        embl = f"{res}/gubbins/{{ref}}.concatenated.recombination_predictions.embl",
        sentinel = f"{res}/gubbins/{{ref}}.gubbins.done"
    output:
        f"{res}/gubbins/{{ref}}.gubbins.bed"
    shell:
        "workflow/scripts/make_bed_from_gubbins.py {input.embl} "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.sentinel} > {output}"

rule merge_vcf:
    conda:
        "envs/bwa.yaml"
    threads: 8
    input:
        vcfs = lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{{ref}}/RG_SC_RA_filter.hvar.vcf.gz".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s
            ) for s in filt_samples if ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ],
        indices = lambda wildcards: [
            "{res}/{s}/MRCA_ref_mapping/{{ref}}/RG_SC_RA_filter.hvar.vcf.gz.csi".format(
                res=res, ref=ref_from_QC(s, f"{res}/QC_summary.csv"), s=s
            ) for s in filt_samples if ref_from_QC(s, f"{res}/QC_summary.csv") == wildcards.ref
        ],
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.vcf.gz"
    shell:
        "bcftools merge -0 --threads {threads} {input.vcfs} -Oz -o {output} "

rule filter_merged_with_gubbins_vcf:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        vcf = f"{res}/MRCA_ref_mapping/{{ref}}/merge.vcf.gz",
        merge_bed = f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.bed",
        gubbins_bed =  f"{res}/gubbins/{{ref}}.gubbins.bed",
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.gubbins.vcf.gz"
    shell:
        "bcftools view {input.vcf} -T ^{input.merge_bed} -Ov | "
        "bcftools view - -T ^{input.gubbins_bed} -Oz -o {output}"

rule filter_merged_without_gubbins_vcf:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        vcf = f"{res}/MRCA_ref_mapping/{{ref}}/merge.vcf.gz",
        merge_bed = f"{res}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.merge.bed",
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.nogubbins.vcf.gz"
    shell:
        "bcftools view {input.vcf} -T ^{input.merge_bed} -Oz -o {output}"


rule make_SNP_alignment:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.{{gubbins}}.vcf.gz"
    output:
        f"{res}/SNP_phylo/{{ref}}.merge.{{gubbins,\w+}}.fasta",
    shell:
        "workflow/scripts/make_SNP_alignment.py {input} > {output}"


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
        "envs/snpeff.yaml"
    threads: 1
    input:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.{{gubbins}}.vcf.gz"
    params:
        db = lambda wildcards: snpEff_db(wildcards.ref)
    output:
        f"{res}/MRCA_ref_mapping/{{ref}}/merge.{{gubbins,\w+}}.snpeff.vcf"
    shell:
        "snpEff {params.db} {input} > {output}"


rule MRCA_SNPs_phylogeny:
    conda:
        "envs/iqtree.yaml"
    threads: 8
    input:
        f"{res}/SNP_phylo/{{ref}}.merge.{{gubbins}}.fasta"
    output:
        tree = f"{res}/SNP_phylo/{{ref}}.merge.{{gubbins,\w+}}.fasta.treefile"
    shell:
        "iqtree -s {input} -nt {threads} -m MFP -B 10000 -redo --nmax 2000"

rule MRCA_SNPs_matrix:
    conda:
        "envs/phy_plots.yaml"
    threads: 1
    input:
        f"{res}/SNP_phylo/{{ref}}.merge.{{gubbins}}.fasta"
    output:
        f"{res}/SNP_phylo/{{ref}}.merge.{{gubbins}}.fasta.snpdists.csv"
    shell:
        "snp-dists -bc {input} > {output}"
