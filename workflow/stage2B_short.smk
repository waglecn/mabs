######################################################################
#
# STAGE 2B SHORT RULES - mapping short reads to internal reference
#
######################################################################

# TODO think about moving these to Snakefile to collect all high level rules?
rule stage2B_short_all_outputs:
    input:
        f"{res}/internal_reference_map/iref.fasta",
        expand(
            "{res}/{sample}/internal_ref_mapping/temp.merged.sorted.bam",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.merged.sorted.bam",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.mpileup",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC_filter.vcf.gz",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC_filter.hvar.vcf.gz",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC_DF.bed",
            sample=filt_samples, res=res
        ),
        expand(
            "{res}/internal_reference_map/merged.RG_SC.consensus.fasta.snpdists.csv",
            res=res
        ),
        expand(
            "{res}/internal_reference_map/merged.RG_SC.consensus.fasta",
            res=res
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.consensus.fa",
            res=res, sample=filt_samples
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.mask.bed",
            res=res, sample=filt_samples
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.variants.bed",
            res=res, sample=filt_samples
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.lowcov.bed",
            res=res, sample=filt_samples
        ),
        expand(
            "{res}/{sample}/internal_ref_mapping/RG_SC.0cov.bed",
            res=res, sample=filt_samples
        )

# TODO candidate to move to common.smk
rule prepare_internal_reference:
    threads: 1
    input:
        f"{res}/{config['internal_ref']}/dflye_short_polish/contigs.fa"
    output:
        f"{res}/internal_reference_map/iref.fasta"
    shell:
        f"cat {{input}} | sed -e 's/>/>{config['internal_ref']}/' > "
        "{output} ;"

rule short_internal_map_to_reference:
    threads: 8
    conda:
        "envs/bwa.yaml"
    input:
        R1 = f"{res}/{{sample}}/input/R1.trim.fastq.gz",
        R2 = f"{res}/{{sample}}/input/R2.trim.fastq.gz",
        S1 = f"{res}/{{sample}}/input/S1.trim.fastq.gz",
        idx = f"{res}/internal_reference_map/iref.fasta.amb"
    output:
        paired_temp = temp(
            f"{res}/{{sample}}/internal_ref_mapping/paired.sorted.bam"
        ),
        single_temp = temp(
            f"{res}/{{sample}}/internal_ref_mapping/single.sorted.bam"
        ),
        merge_temp = temp(f"{res}/{{sample}}/internal_ref_mapping/merge.bam"),
        merge_sorted = temp(
            f"{res}/{{sample}}/internal_ref_mapping/temp.merged.sorted.bam"
        ),
    shell:
        "bwa-mem2 mem -t {threads} -M "
        f"{res}/internal_reference_map/iref.fasta "
        "{input.R1} {input.R2} | samtools view -Sbh - | samtools sort > "
        "{output.paired_temp} ;"
        "bwa-mem2 mem -t {threads} -M "
        f"{res}/internal_reference_map/iref.fasta "
        "{input.S1} | samtools view -Sbh - | samtools sort -@ {threads} > "
        "{output.single_temp} ;"
        "samtools index {output.paired_temp} ;"
        "samtools index {output.single_temp} ;"
        "samtools merge -f {output.merge_temp} {output.paired_temp} "
        "{output.single_temp} ;"
        "samtools view -hF 4 {output.merge_temp} | samtools sort -@ {threads} > "
        "{output.merge_sorted} ;"

rule short_internal_add_rg_to_reference_bam:
    threads: 2
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/temp.merged.sorted.bam"
    output:
        temp(f"{res}/{{sample}}/internal_ref_mapping/RG.bam")
    shell:
        "samtools addreplacerg -r \"@RG\tID:{wildcards.sample}\" -o {output} "
        "{input}"

rule short_internal_ref_softclip_filter:
    conda:
        "envs/bwa.yaml"
    threads: 1
    input:
        f"{res}/{{sample}}/internal_ref_mapping/RG.bam"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/RG_SC.merged.sorted.bam",
    shell:
        "workflow/scripts/sclips.py filter {input} > {output} "

rule short_internal_make_mpileup:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/RG_SC.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/RG_SC.mpileup"
    shell:
        # note that mpileup has moved to bcftools
        # -d max depth
        # -q min mapQ
        # -t in samtools 1.7 -> -a output tags
        # -u in samtools 1.7 -> -O v uncompressed VCF output
        # -g in samtools 1.7 -> bcftools includes genotype likelihoods default
        "bcftools mpileup -d 1000 -q 30 -a DP,AD,ADF,ADR,SP -Oz "
        f"-f {res}/internal_reference_map/iref.fasta "
        "{input} -o {output}"

rule short_internal_call_and_filter_variants:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.mpileup"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.vcf.gz"
    shell:
        # -Ov output uncompressed vcf
        # -m multiallelic caller
        # -v variants only
        "bcftools call --ploidy 1 -Ov -m -v {input} | "
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

rule short_internal_filter_hsnps:
    threads: 1
    conda:
        "envs/bwa.yaml"
    params:
        snp_cutoff = 0.90
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.vcf.gz"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.hvar.vcf.gz"
    shell:
        "bcftools filter -i '(AD[0:1]/(AD[0:0]+AD[0:1]) > "
        "{params.snp_cutoff})' -Oz -o {output} {input} "

# make per-sample beds for merged filtering
rule short_internal_make_internal_bed_0cov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.0cov.bed"
    shell:
        "bedtools genomecov -ibam {input} -bga | awk '$4==0' > {output}"

rule short_internal_make_bed_lowcov:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.merged.sorted.bam"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.lowcov.bed"
    shell:
        "bedtools genomecov -bga -ibam {input} | awk '$4 < 10' > {output}"

rule short_internal_make_variants_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.hvar.vcf.gz"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.variants.bed"
    shell:
        "bcftools query -f'%CHROM\t%POS0\t%END\n' {input} > {output}"

rule short_internal_make_mask_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        varbed = f"{res}/{{sample}}/internal_ref_mapping/{{step}}.variants.bed",
        lowbed = f"{res}/{{sample}}/internal_ref_mapping/{{step}}.lowcov.bed"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.mask.bed"
    shell:
        "bedtools subtract -a {input.lowbed} -b {input.varbed} > {output}"


rule short_internal_make_simple_consensus:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        ref = f"{res}/internal_reference_map/iref.fasta",
        vcf = f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.hvar.vcf.gz",
        idx = f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.hvar.vcf.gz.csi",
        mask = f"{res}/{{sample}}/internal_ref_mapping/{{step}}.mask.bed"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}.consensus.fa"
    shell:
        f"bcftools consensus -p {{wildcards.sample}} "
        "-f {res}/internal_reference_map/iref.fasta --mark-del '-' "
        "-m {input.mask} -i 'INFO/MQ >= 20 & FORMAT/DP >= 10' {input.vcf} | "
        "sed \"/^>/s/{wildcards.sample}.*/{wildcards.sample}/\" > {output}"

rule short_internal_merge_simple_consensus:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        expand(
            "{res}/{sample}/internal_ref_mapping/{{step}}.consensus.fa",
            res=res, sample=filt_samples
        ),
        f"{res}/internal_reference_map/iref.fasta"
    output:
        f"{res}/internal_reference_map/merged.{{step}}.consensus.fasta"
    shell:
        "cat {input} > {output}"

rule short_internal_filter_consensus:
    threads: 1
    conda:
        "envs/phy_plots.yaml"
    input:
        f"{res}/internal_reference_map/merged.{{step}}.consensus.fasta"
    output:
        f"{res}/internal_reference_map/merged.{{step}}.consensus.fasta.snpdists.csv"
    shell:
        "snp-dists -bc {input} > {output}"


rule short_internal_density_filter_bed:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_filter.hvar.vcf.gz"
    output:
        f"{res}/{{sample}}/internal_ref_mapping/{{step}}_DF.bed"
    shell:
        "workflow/scripts/make_vcf_density_bed.py {input} > {output}"

# rule to get PE_PPE cds from prokka??

# from https://github.com/samtools/bcftools/issues/1386#issuecomment-870134547
# bcftools query -f'%CHROM\t%POS0\t%END\n' ${NAME}.bcftools.vcf.gz > variants.bed
# # get low coverage sites in bedgraph format
# bedtools genomecov -bga -ibam $BAM | awk '$4 < 10' > low_coverage_sites.bed
# bedtools subtract -a low_coverage_sites.bed -b variants.bed > mask.bed
# #generate consensus and rename header
# bcftools consensus -p $PREFIX -f $REFERENCE --mark-del '-' -m mask.bed -i 'FORMAT/VAF > 0.90 & INFO/MQ >= 20 & FORMAT/DP >= 10' ${NAME}.bcftools.vcf.gz | sed "/^>/s/${PREFIX}.*/${PREFIX}/" > ${NAME}.bcftools.consensus.fa
