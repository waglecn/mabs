######################################################################
#
# STAGE 2 RULES - short reads
#
######################################################################
rule short_temp_MRCA_ref_alignment:
    priority: 10
    threads: 8
    conda:
        "envs/bwa.yaml"
    input:
        QC = f"{res}/QC_summary.csv",
        R1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz",
        S1 = f"{res}/{{sample}}/input/S1.trim.fastq.gz",
        idx = "workflow/resources/alignment_references/{ref}.fasta.amb"
    output:
        paired_temp = temp(
            f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/paired.sorted.bam"
        ),
        single_temp = temp(
            f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/single.sorted.bam"
        ),
        merge_temp = temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/merge.bam"),
        merge_sorted = temp(
            f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/temp.merged.sorted.bam"
        ),
    shell:
        "bwa-mem2 mem -t {threads} -M "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa-mem2 mem -t {threads} -M "
        "workflow/resources/alignment_references/{wildcards.ref}.fasta "
        "{input.S1} | samtools view -Sbh - | samtools sort -@ {threads} > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort -@ {threads} > "
        "{output.merge_sorted} ;"

rule short_temp_add_read_groups:
    threads: 1
    input:
        bam = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/temp.merged.sorted.bam",
        bai = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/temp.merged.sorted.bam.bai"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRG.merged.sorted.bam")
    log:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG.log"
    shell:
        "samtools addreplacerg -r "
        "'@RG\tID:{wildcards.sample}\tSM:??????\tLB:??????\tPL:ILLUMINA' "
        "-o {output} {input.bam} | tee 2> {log}"

rule short_MRCA_ref_softclip_filter:
    conda: "envs/bwa.yaml"
    threads: 1
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRG.merged.sorted.bam",
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRG.merged.sorted.bam.bai"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.merged.sorted.bam"),
    shell:
        "workflow/scripts/sclips.py filter {input} > {output} "

rule short_MRCA_make_mpileup:
    threads: 1
    conda:
        "envs/bwa.yaml"    
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC.mpileup"
    shell:
        # note that mpileup has moved to bcftools
        # -d max depth
        # -q min mapQ
        # -t in samtools 1.7 -> -a output tags
        # -u in samtools 1.7 -> -O v uncompressed VCF output
        # -g in samtools 1.7 -> bcftools includes genotype likelihoods default
        "bcftools mpileup -d 1000 -q 30 -a DP,AD,ADF,ADR,SP -Oz "
        "-f workflow/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "{input} -o {output}"


rule short_MRCA_call_and_filter_variants:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.vcf.gz"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call --ploidy 1 -Ov -m -v --ploidy 1 {input} | "
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

rule short_MRCA_inverse_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.failed.vcf.gz"
    shell:
        "bcftools call --ploidy 1 -Oz -m -v {input} | "
        "bcftools filter -i 'SP>=45 || MQ<=30 || FORMAT/DP<=10 || QUAL<=50' "
        "-o {output}"

rule short_MRCA_inverse_AD_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.AD_failed.vcf.gz"
    params:
        snp_cutoff = 0.90
    shell:
        "bcftools call --ploidy 1 -Oz -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) < "
        "{params.snp_cutoff}) & (ADF[0:1]<=1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule short_filter_hsnps:
    threads: 1
    conda:
        "envs/bwa.yaml"
    params:
        snp_cutoff = 0.90
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.vcf.gz"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input} "

# make per-sample beds for merged filtering
rule short_make_bed_0cov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

rule short_make_bed_lowcov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.lowcov.bed"
    shell:
        "bedtools genomecov -bga -ibam {input} | awk '$4 < 10' > {output}"

rule short_make_variants_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.variants.bed"
    shell:
        "bcftools query -f'%CHROM\t%POS0\t%END\n' {input} > {output}"

rule short_make_mask_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        varbed=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.variants.bed",
        lowbed=f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.lowcov.bed"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.mask.bed"
    shell:
        "bedtools subtract -a {input.lowbed} -b {input.varbed} > {output}"

rule short_make_simple_consensus:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        ref = "workflow/resources/alignment_references/{ref}.fasta",
        vcf = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz",
        ids = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz.csi",
        mask = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.mask.bed"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}.consensus.fa"
    shell:
        f"bcftools consensus -p {{wildcards.sample}} "
        "-f {input.ref} --mark-del '-' "
        "-m {input.mask} -i 'INFO/MQ >= 20 & FORMAT/DP >= 10' {input.vcf} | "
        "sed \"/^>/s/{wildcards.sample}.*/{wildcards.sample}/\" > {output}"

rule short_density_filter_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz",
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar.vcf.gz.csi"
    output:
        bed = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step,[A-Z_]+}}_filter.hvar_DF.bed"
    shell:
        "workflow/scripts/make_vcf_density_bed.py {input} > {output.bed}"


