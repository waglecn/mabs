def ref_from_QC(sample_name, QC_summary_file="QC_summary.csv", data="short"):
    import sys
    '''
    Every rule using this function must require QC_summary.csv in input files
    '''
    ref = 'Unknown'
    try:
        QC_data = pd.read_csv(QC_summary_file)
        QC_data.index = QC_data['sample']
        if data == 'short':
            ref = QC_data.loc[sample_name]['MRCA_ref']
        elif data == 'long':
            ref = QC_data.loc[sample_name]['MRCA_long_ref']
    except FileNotFoundError:
        print('QC summary file not found', file=sys.stderr)
    except KeyError:
        print(f'sample {sample_name} not in QC summary', file=sys.stderr)
    return ref


def get_samples():
    import yaml
    print(
        "looking for samples defined in: {}".format(config['sample_sheet']),
        file=sys.stderr
    )
    with open(config["sample_sheet"], "r") as stream:
        try:
            samplesheet = yaml.safe_load(stream)
            for k in samplesheet.keys():
                if 'minion_path' in samplesheet[k]:
                    samplesheet[k]['MINION'] = []
                    for p in samplesheet[k]['minion_path']:
                        samplesheet[k]['MINION'] += glob.glob(
                            config['data_dir'] + p
                        )
        except yaml.YAMLError as exc:
            print(exc)
        return samplesheet

######################################################################
#
# Utility RULES
#
######################################################################

rule index_bam:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        "{inbam}.bam"
    output:
        temp("{inbam}.bam.bai")
    shell:
        "samtools index {input} "

rule index_vcf:
    threads: 1
    conda:
        "envs/bwa.yaml"
    input:
        "{invcf}.vcf.gz"
    output:
        "{invcf}.vcf.gz.csi"
    shell:
        "htsfile {input} ;"
        "bcftools index {input}"

rule create_fasta_dict_index:
    threads: 1
    conda:
        "envs/picard.yaml"
    input:
        "workflow/resources/alignment_references/{ref}.fasta"
    output:
        seqdict = "workflow/resources/alignment_references/{ref}.dict",
        fai = "workflow/resources/alignment_references/{ref}.fasta.fai"
    shell:
        "picard CreateSequenceDictionary -R {input} -O {output.seqdict} && "
        "samtools faidx {input} -o {output.fai}"

# make ref BED for PE  PPE genes
rule make_PE_PPE_BED:
    threads: 1
    conda:
        "envs/phy_plots.yaml"
    input:
        "workflow/resources/alignment_references/{ref}.gbk"
    output:
        "workflow/resources/alignment_references/{ref}.PE_PPE.bed"
    shell:
        "workflow/scripts/make_PE_PPE_BED.py {input} > {output}"

###########################################
# Cleanup rules
###########################################

rule clean_gubbins_tmp:
    threads: 1
    log:
        f"{res}/logs/clean_gubbins.log"
    shell:
        "rm -Ifr tmp*"

rule clean_mashtree_tmp:
    threads: 1
    log:
        f"{res}/logs/clean_mashtree_tmp.log"
    shell:
        "rm -Ifr mashtree_tmp"

rule cleanup:
    input:
        rules.clean_gubbins_tmp.log,
        rules.clean_mashtree_tmp.log
        
