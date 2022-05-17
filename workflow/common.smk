def ref_from_QC(sample_name, QC_summary_file="QC_summary.csv"):
    import sys
    '''
    Every rule using this function must require QC_summary.csv in input files
    '''
    ref = 'Unknown'
    try:
        QC_data = pd.read_csv(QC_summary_file)
        QC_data.index = QC_data['sample']
        ref = QC_data.loc[sample_name]['MRCA_ref']
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
