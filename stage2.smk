######################################################################
#
# STAGE 2 RULES
#
######################################################################
rule temp_MRCA_ref_alignment:
    priority: 10
    threads: 8
    params:
        execdir = exec_dir
    input:
        QC = 'QC_summary.csv',
        S1 = "trimmed_input/{sample}.S1.fastq.gz",
        R1 = "trimmed_input/{sample}.R1.fastq.gz",
        R2 = "trimmed_input/{sample}.R2.fastq.gz",
    output:
        paired_temp = temp(
            "MRCA_ref_mapping/{ref}/{sample}.paired.sorted.bam"
        ),
        single_temp = temp(
            "MRCA_ref_mapping/{ref}/{sample}.single.sorted.bam"
        ),
        merge_temp = temp("MRCA_ref_mapping/{ref}/{sample}.merge.bam"),
        merge_sorted = temp(
            "MRCA_ref_mapping/{ref}/temp.merged.{sample}.sorted.bam"
        ),
    shell:
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/{wildcards.ref} "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa mem -t {threads} -M "
        "{params.execdir}/resources/alignment_references/{wildcards.ref} "
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
    conda: "conda_envs/picard.yaml"
    input:
        bam = "MRCA_ref_mapping/{ref}/temp.merged.{sample}.sorted.bam",
        bai = "MRCA_ref_mapping/{ref}/temp.merged.{sample}.sorted.bam.bai"
    output:
        temp("MRCA_ref_mapping/{ref}/tempRG.merged.{sample}.sorted.bam")
    log:
        "MRCA_ref_mapping/{ref}/{sample}.picard.RG.log"
    shell:
        "picard AddOrReplaceReadGroups -INPUT {input.bam} "
        "-OUTPUT {output} -SORT_ORDER coordinate -RGID {wildcards.sample} "
        "-RGLB unknown -RGPL Illumina -RGSM {wildcards.sample} "
        "-RGPU project -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT "
        " 2>&1 | tee {log}"

rule MRCA_ref_softclip_filter:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}/tempRG.merged.{sample}.sorted.bam"
    output:
        temp("MRCA_ref_mapping/{ref}/tempRGSC.merged.{sample}.sorted.bam")
    shell:
        "{params.execdir}/scripts/sclips.py filter {input} > {output}"

rule MRCA_ref_gatk_realignment_intervals:
    threads: 1
    conda:
        "conda_envs/gatk3.yaml"
    params:
        execdir = exec_dir
    input:
        bam = "MRCA_ref_mapping/{ref}/tempRGSC.merged.{sample}.sorted.bam",
        bai = "MRCA_ref_mapping/{ref}/tempRGSC.merged.{sample}.sorted.bam.bai",
        seqdict = f"{exec_dir}/resources/alignment_references/{{ref}}.dict",
        fai = f"{exec_dir}/resources/alignment_references/{{ref}}.fasta.fai",
        gatk = f"{exec_dir}/resources/gatk-registered"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals"
    log:
        "MRCA_ref_mapping/{ref}/{sample}.gatk3_intervals.log"
    shell:
        "gatk3 -T RealignerTargetCreator -R "
        "{params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -o {output} 2>&1 | tee {log}"

rule MRCA_ref_gatk_realignment:
    threads: 1
    params:
        execdir = exec_dir
    conda:
        "conda_envs/gatk3.yaml"
    input:
        bam = "MRCA_ref_mapping/{ref}/tempRGSC.merged.{sample}.sorted.bam",
        bai = "MRCA_ref_mapping/{ref}/tempRGSC.merged.{sample}.sorted.bam.bai",
        intervals = "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.intervals",
    output:
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
    log:
        "MRCA_ref_mapping/{ref}/{sample}.gatk_realign.log"
    shell:
        "gatk3 -rf OverclippedRead "
        # "--filter_is_too_short_value 100 -minRead 100 "
        # "--do_not_require_softclips_both_ends -rf ReadLength -maxRead 500 "
        "-T IndelRealigner "
        "-R {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "-I {input.bam} -targetIntervals {input.intervals} -o {output} "
        " 2>&1 | tee {log}"

rule MRCA_make_mpileup:
    threads: 1
    params:
        execdir = exec_dir
    input:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.bam"
    output:
        # "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA_mq30_baq.mpileup"
        "MRCA_ref_mapping/{ref}/{sample}.RG_SC_RA.mpileup"
    shell:
        # note that mpileup has moved to bcftools
        # -d max depth
        # -q min mapQ
        # -t in samtools 1.7 -> -a output tags
        # -u in samtools 1.7 -> -O v uncompressed VCF output
        # -g in samtools 1.7 -> bcftools includes genotype likelihoods default
        "bcftools mpileup -d 1000 -q 30 -a DP,AD,ADF,ADR,SP -Oz "
        "-f {params.execdir}/resources/alignment_references/"
        "{wildcards.ref}.fasta "
        "{input} -o {output}"


rule MRCA_call_and_filter_variants:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_filter.vcf.gz"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call -Ov -m -v {input} | "
        # filtering
        # MQ mapping quality
        # ID-SP - Phred strand bias p-value
        # ID-ADF[1] - allelic depth forward strand of first alt allele
        # ID-ADR[1] - alleleic depth reverse strand of first alt allele
        # QUAL - Phred scaled ALT quality
        # FORMAT-DP - number of high-quality bases
        # FORMAT-SP - Phred strand bias p-value
        # FORMAT-ADF and ADR as in INFO
        # Note need to specify sample 0:
        # see https://github.com/samtools/bcftools/issues/757
        "bcftools filter -i 'SP<45 & ADF[0:1]>1 & ADR[0:1]>1 & MQ>30 & "
        "QUAL>50 & FORMAT/DP > 10 & SP<45 & ADF[0:1]>1 & ADR[0:1]>1' "
        "-Oz -o {output}"

rule MRCA_inverse_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_filter.failed.vcf.gz"
    shell:
        "bcftools call -Oz -m -v {input} | "
        "bcftools filter -i 'SP>=45 || MQ<=30 || FORMAT/DP<=10 || QUAL<=50' "
        "-o {output}"

rule MRCA_inverse_AD_filter:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.mpileup"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_filter.AD_failed.vcf.gz"
    params:
        snp_cutoff = 0.90
    shell:
        "bcftools call -Oz -m -v {input} | "
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) < "
        "{params.snp_cutoff}) & (ADF[0:1]<=1 || ADR[0:1]<=1)' -o {output}"

# density rule to filter out snps
rule filter_hsnps:
    threads: 1
    params:
        snp_cutoff = 0.90
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_filter.vcf.gz"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}_filter.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input} "

# make per-sample beds for merged filtering
rule make_bed_0cov:
    threads: 1
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.bam"
    output:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

rule density_filter_bed:
    threads: 1
    params:
        execdir = exec_dir
    input:
        "MRCA_ref_mapping/{ref}/{sample}.{step}.vcf.gz"
    output:
        bed = "MRCA_ref_mapping/{ref}/{sample}.{step}_DF.bed"
    shell:
        "{params.execdir}/scripts/make_vcf_density_bed.py {input} > "
        "{output.bed}"

# make ref BED for PE  PPE genes
rule make_PE_PPE_BED:
    threads: 1
    params:
        execdir = exec_dir
    conda:
        "conda_envs/phy_plots.yaml"
    input:
        f"{exec_dir}/resources/alignment_references/{{ref}}.gbk"
    output:
        "{params.execdir}/resources/alignment_references/{ref}.PE_PPE.bed"
    shell:
        "{params.execdir}/scripts/make_PE_PPE_BED.py {input} > {output}"