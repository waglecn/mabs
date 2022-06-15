######################################################################
#
# STAGE 2 RULES
#
######################################################################
rule temp_MRCA_ref_alignment:
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

rule temp_add_read_groups:
    threads: 1
    conda:
        "envs/picard.yaml"
    input:
        bam = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/temp.merged.sorted.bam",
        bai = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/temp.merged.sorted.bam.bai"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRG.merged.sorted.bam")
    log:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/picard.RG.log"
    shell:
        "picard AddOrReplaceReadGroups -INPUT {input.bam} "
        "-OUTPUT {output} -SORT_ORDER coordinate -RGID {wildcards.sample} "
        "-RGLB unknown -RGPL Illumina -RGSM {wildcards.sample} "
        "-RGPU project -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT "
        " 2>&1 | tee {log}"

rule MRCA_ref_softclip_filter:
    conda: "envs/bwa.yaml"
    threads: 1
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRG.merged.sorted.bam"
    output:
        temp(f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRGSC.merged.sorted.bam"),
    shell:
        "workflow/scripts/sclips.py filter {input} > {output} "


rule MRCA_ref_gatk_realignment_intervals:
    threads: 1
    conda:
        "envs/gatk3.yaml"
    input:
        bam = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRGSC.merged.sorted.bam",
        bai = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRGSC.merged.sorted.bam.bai",
        seqdict = "workflow/resources/alignment_references/{ref}.dict",
        fai = "workflow/resources/alignment_references/{ref}.fasta.fai",
        # gatk = "workflow/resources/gatk-registered"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.intervals"
    log:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/gatk3_intervals.log"
    shell:
        "gatk3 -nt 1 -T RealignerTargetCreator -R "
        "workflow/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -o {output} 2>&1 | tee {log}"

rule MRCA_ref_gatk_realignment:
    threads: 1
    conda:
        "envs/gatk3.yaml"
    input:
        bam = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRGSC.merged.sorted.bam",
        bai = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/tempRGSC.merged.sorted.bam.bai",
        intervals = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.intervals",
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.bam"
    log:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/gatk_realign.log"
    shell:
        "gatk3 -nt 1 -rf OverclippedRead "
        # "--filter_is_too_short_value 100 -minRead 100 "
        # "--do_not_require_softclips_both_ends -rf ReadLength -maxRead 500 "
        "-T IndelRealigner "
        "-R workflow/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -targetIntervals {input.intervals} -o {output} "
        " 2>&1 | tee {log}"

rule MRCA_make_mpileup:
    threads: 1
    conda:
        "envs/bwa.yaml"    
    input:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.bam"
    output:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_mq30_baq.mpileup"
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/RG_SC_RA.mpileup"
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


rule MRCA_call_and_filter_variants:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.vcf.gz"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call -Ov -m -v {input} | "
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

rule MRCA_inverse_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.failed.vcf.gz"
    shell:
        "bcftools call -Oz -m -v {input} | "
        "bcftools filter -i 'SP>=45 || MQ<=30 || FORMAT/DP<=10 || QUAL<=50' "
        "-o {output}"

rule MRCA_inverse_AD_filter:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.mpileup"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.AD_failed.vcf.gz"
    params:
        snp_cutoff = 0.90
    shell:
        "bcftools call -Oz -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) < "
        "{params.snp_cutoff}) & (ADF[0:1]<=1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule filter_hsnps:
    threads: 1
    conda:
        "envs/bwa.yaml"
    params:
        snp_cutoff = 0.90
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.vcf.gz"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_filter.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input} "

# make per-sample beds for merged filtering
rule make_bed_0cov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.bam"
    output:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

rule density_filter_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}.vcf.gz"
    output:
        bed = f"{res}/{{sample}}/MRCA_ref_mapping/{{ref}}/{{step}}_DF.bed"
    shell:
        "workflow/scripts/make_vcf_density_bed.py {input} > "
        "{output.bed}"

# make ref BED for PE  PPE genes
rule make_PE_PPE_BED:
    threads: 1
    conda:
        "envs/phy_plots.yaml"
    input:
        "workflow/resources/alignment_references/{{ref}}.gbk"
    output:
        "workflow/resources/alignment_references/{ref}.PE_PPE.bed"
    shell:
        "workflow/scripts/make_PE_PPE_BED.py {input} > {output}"