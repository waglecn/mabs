# def MRCA_mapped_ref_input(sample_name):
#     '''
#     Use as input function that maps sample_name to reference
#         input:
#             sample_name: key should be found in sample_dict
#         returns:
#             ref_name: should match available references
#     '''
#     from scripts.tree_MRCA import tree_MRCA
#     MRCA = tree_MRCA(
#         '{}/mashtree/assembly_mashtree.complete.tree'.format(
#             config['results_dir']
#         ),
#         sample_name
#     )
#     # print(sample_name, MRCA, file=sys.stderr)
#     return MRCA


# def ref_from_QC(sample_name, QC_summary_file="QC_summary.csv"):
#     '''
#     Every rule using this function must require QC_summary.csv in input files
#     '''
#     ref = 'Unknown'
#     try:
#         QC_data = pd.read_csv(QC_summary_file)
#         QC_data.index = QC_data['sample']
#         ref = QC_data.loc[sample_name]['MRCA_ref']
#     except FileNotFoundError:
#         print('QC summary file not found', file=sys.stderr)
#     except KeyError:
#         print(f'sample {sample_name} not in QC summary', file=sys.stderr)
#     return ref

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
                            config['data_dir'] + "/" + p
                        )
        except yaml.YAMLError as exc:
            print(exc)
        return samplesheet
