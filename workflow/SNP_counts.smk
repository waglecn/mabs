# rules for counting SNPS
#

rule count_snps:
    input:
        [
            "SNP_counts/{ref}/{sample}.vcf_counts.csv".format(
                ref=ref_from_QC(s), sample=s
            ) for s in filt_samples
        ],

        # # 1. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 2. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 3. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar_0cov.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.0cov.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 4. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar_DF.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.DF.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 5. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_fo;ter.hvar.0COVDF.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.0covDF.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 6. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.0COVDFPE.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.0covDFPE.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 7. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.0COVPE.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.0covPE.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 8. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar_gubbins.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.gubbins.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],
        # # 9. "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar_0COVPE_gubbins.vcf.gz"
        # [
        #     "SNP_counts/{ref}/{sample}.RG_SC_RA_filter.hvar.0cov_gubbinsPE.count".format(
        #         ref=ref_from_QC(s), sample=s
        #     ) for s in filt_samples
        # ],

rule count_vcf_snps:
    threads: 1
    params:
        execdir = exec_dir
    input:
        base = "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.vcf.gz",
        hvar = "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.hvar.vcf.gz",
        0cov = "MRCA_ref_mapping/{ref}.RG_SC_RA_merge.0cov.bed",
        DF = "MRCA_ref_mapping/{ref}.RG_SC_RA.merge.DF.bed",
        PE = "{params.execdir}/resources/alignment_references/{ref}.PE_PPE.bed",
        gubbins = "gubbins/{ref}.gubbins.bed"
    output:
        "SNP_counts/{ref}/{s}.vcf_counts.csv"
    shell:
        "{params.execdir}/scripts/vcf_count_summary.py "
        "{input.base} {input.hvar} {input.0cov} {input.DF} {input.DF} "
        "{input.PE} {input.gubbins} > {output} "

rule count_1:
    threads: 1
    params:
        execdir = exec_dir
    input:
        base = "MRCA_ref_mapping/{ref}/{s}.RG_SC_RA_filter.vcf.gz",
    output:
        vcf1 = "SNP_counts/{ref}/{s}.1.csv",
    shell:
        "{params.execdir}/scripts/count_vcf_bed_vcfs.py {input} > {output} "
