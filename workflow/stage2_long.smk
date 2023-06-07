######################################################################
#
# STAGE 2 RULES - long reads
#
######################################################################
rule long_temp_MRCA_ref_alignment:
    priority: 10
    threads: 8
    conda:
        "envs/bwa.yaml"
    input:
        QC = f"{res}/QC_summary.csv",
        reads = f"{res}/{{sample}}/input/long.fastq.gz",
    output:
        bam = temp(
            f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/longmerged.sorted.bam"
        ),
    shell:
        "minimap2 -x map-ont -t {threads} -a "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.reads} | samtools view -Sbh - | samtools sort > "
        "{output.bam} "

rule long_temp_add_read_groups:
    threads: 1
    input:
        bam = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/longmerged.sorted.bam",
        bai = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/longmerged.sorted.bam.bai"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG.longmerged.sorted.bam")
    log:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/longRG.log"
    shell:
        "samtools addreplacerg -r "
        "'@RG\tID:{wildcards.sample}\tSM:??????\tLB:??????\tPL:MINION' "
        "-o {output} {input.bam} | tee 2> {log}"

rule long_MRCA_ref_softclip_filter:
    conda: "envs/bwa.yaml"
    threads: 1
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG.longmerged.sorted.bam"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.longmerged.sorted.bam"),
    shell:
        "workflow/scripts/sclips.py filter {input} > {output} "

rule long_MRCA_make_mpileup:
    threads: 1
    conda:
        "envs/bwa.yaml"    
    input:
        bam=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.longmerged.sorted.bam",
        bai=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.longmerged.sorted.bam.bai"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.long.mpileup"
    shell:
        # note that mpileup has moved to bcftools
        # -d max depth
        # -q min mapQ
        # -t in samtools 1.7 -> -a output tags
        # -u in samtools 1.7 -> -O v uncompressed VCF output
        # -g in samtools 1.7 -> bcftools includes genotype likelihoods default
        "bcftools mpileup -d 1000 -q 30 -a DP,AD,ADF,ADR,SP -Oz --config ont "
        "-f workflow/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "{input.bam} -o {output}"


rule long_MRCA_call_and_filter_variants:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.vcf.gz"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call -P 0.01 --ploidy 1 -Ov -m -v --ploidy 1 {input} | "
        # filtering
        # MQ mapping quality
        # ID-SP - Phred-scaled strand bias p-value
        # ID-ADF[1] - alleleic depth forward strand of first alt allele
        # ID-ADR[1] - alleleic depth reverse strand of first alt allele
        # QUAL - Phred-scaled ALT quality
        # FORMAT-DP - number of high-quality bases
        # FORMAT-SP - Phred-scaled strand bias p-value
        # FORMAT-ADF and ADR as in INFO
        # Note need to specify sample 0:
        # see https://github.com/samtools/bcftools/issues/757
        "bcftools filter -i 'SP<45 & ADF[0:1]>1 & ADR[0:1]>1 & MQ>30 & "
        "QUAL>50 & FORMAT/DP > 10 & SP<45 & ADF[0:1]>1 & ADR[0:1]>1' "
        "-Oz -o {output}"

rule long_MRCA_inverse_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.failed.vcf.gz"
    shell:
        "bcftools call --ploidy 1 -Oz -m -v {input} | "
        "bcftools filter -i 'SP>=45 || MQ<=30 || FORMAT/DP<=10 || QUAL<=50' "
        "-o {output}"

rule long_MRCA_inverse_AD_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.AD_failed.vcf.gz"
    params:
        snp_cutoff = 0.90
    shell:
        "bcftools call --ploidy 1 -Oz -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) < "
        "{params.snp_cutoff}) & (ADF[0:1]<=1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule long_filter_hsnps:
    threads: 1
    conda:
        "envs/bwa.yaml"
    params:
        snp_cutoff = 0.90
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.vcf.gz"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input} "

# make per-sample beds for merged filtering
rule long_make_bed_0cov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.longmerged.sorted.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

rule long_make_bed_lowcov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.longmerged.sorted.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.lowcov.bed"
    shell:
        "bedtools genomecov -bga -ibam {input} | awk '$4 < 10' > {output}"


rule long_make_variants_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar.vcf.gz"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.variants.bed"
    shell:
        "bcftools query -f'%CHROM\t%POS0\t%END\n' {input} > {output}"

rule long_make_mask_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        varbed=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.variants.bed",
        lowbed=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.lowcov.bed"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.mask.bed"
    shell:
        "bedtools subtract -a {input.lowbed} -b {input.varbed} > {output}"

rule long_make_simple_consensus:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        ref = "workflow/resources/alignment_references/{ref}.fasta",
        vcf = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar.vcf.gz",
        ids = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar.vcf.gz.csi",
        mask = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.mask.bed"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.long.consensus.fa"
    shell:
        f"bcftools consensus -p {{wildcards.sample}} "
        "-f {input.ref} --mark-del '-' "
        "-m {input.mask} -i 'INFO/MQ >= 20 & FORMAT/DP >= 10' {input.vcf} | "
        "sed \"/^>/s/{wildcards.sample}.*/{wildcards.sample}/\" > {output}"

rule long_density_filter_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar.vcf.gz"
    output:
        bed = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.long.hvar_DF.bed"
    shell:
        "workflow/scripts/make_vcf_density_bed.py {input} > {output.bed}"

