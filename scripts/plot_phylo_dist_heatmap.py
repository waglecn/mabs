#!/usr/bin/env python3
"""
This script plots a heatmap of SNV distances extracted from a tree.

Example usage:
    plot_phylo_dist_heatmap.py [input.tree]

INPUT:
    FILE: input.tree - a newick format tree file
OUTPUT:
    FILE: input.tree.csv - a csv file with distances extracted from the tree
    STDOUT: a matplotlib heatmap

"""

import sys
import ete3
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def make_matrix(tree):
    leaves = tree.get_leaves()
    matrix = []
    for i, l_1 in enumerate(leaves):
        matrix.append([])
        for j, l_2 in enumerate(leaves):
            matrix[i].append(l_1.get_distance(l_2))
    matrix = np.array(matrix)
    leaves = [l.name for l in leaves]
    return leaves, matrix


def read_sample_csv(csv_file):
    labels = [
        'sample', 'labwarenumber', 'sputum_collected', 'PH_AFB_date',
        'labwarereport', 'seq prefix', 'Notes', 'Date_Recd', 'Run_Start_Date',
        'reads_paired_in', 'reads_paired_out', 'perc_pair_out', 'single_out',
        'perc_single_out',
        'top_sci_name', 'top_percent', 'top_frag_clade', 'top_frag_taxa',
        'top_rank',
        '2nd_sci_name', '2nd_percent', '2nd_frag_clade', '2nd_frag_taxa',
        '2nd_rank',
        'assembly_file', 'num_contigs', 'assembly_length',
        'length_weighted_cov', 'mlst_file', 'mlst_scheme', 'MLST_ST',
        'allele1', 'allele2', 'allele3', 'allele4', 'allele5', 'allele6',
        'allele7', 'subjectID', 'hit_length', 'gaps', 'q_start', 'q_end',
        's_start', 's_end', 'erm41_truncated', 'MRCA', 'num_reads_paired',
        'mean_depth', 'frac_0x', 'frac_>0x', 'frac_>=40x', 'softclipped_bases',
        'QC',
    ]
    samples = pd.read_csv(csv_file, names=labels)
    return samples


def main():
    tree_file = sys.argv[1]
    tree = ete3.Tree(tree_file)
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize(direction=1)
    labels, matrix = make_matrix(tree)

    for node in tree:
        if node.name == 'ref':
            node.name = os.path.split(tree_file)[-1].split('.')[0]

    matrix = pd.DataFrame(matrix, columns=labels, index=labels)
    matrix.to_csv('{}.distances.csv'.format(tree_file))

    ax = sns.heatmap(
        matrix, xticklabels=labels, yticklabels=labels,
        linewidths=0.05
    )
    plt.show()


if __name__ == '__main__':
    main()
