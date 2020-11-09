#!/usr/bin/env python3
"""
This script plots boxplots of distances from different categories extracted
from a datafram of distances

Example usage:
    plot_wp_distances.py

INPUT:
    FILE: ../analysis/mabscessus.RG_SC_RA.merge.fasta.csv - snp-dists csv file
    FILE: ../analysis/mmassiliense.RG_SC_RA.merge.fasta.csv - snp-dists csv file
OUTPUT:
    STDOUT: a matplotlib boxplot
"""

import sys
import ete3
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from plot_phylo_dist_heatmap import make_matrix, read_sample_csv


def same_sample(s1, s2):
    pt1 = s1.split('-')[0]
    pt2 = s2.split('-')[0]
    if pt1 == pt2:
        return True
    return False


def get_distances(matrix, labels, species):
    distances = []
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            distances.append(
                (species, 'all', int(matrix[i][j]))
            )
            if same_sample(labels[i], labels[j]):
                distances.append(
                    (species, 'within_pt', int(matrix[i][j]))
                )
            else:
                distances.append(
                    (species, 'between_pt', int(matrix[i][j]))
                )
    return distances


def main():
    matrix = [
        l.strip().split(',') for l in
        open('../analysis/mabscessus.RG_SC_RA.merge.fasta.csv', 'r')
    ]
    labels = matrix[0][1:]
    matrix = [row[1:] for row in matrix[1:]]

    matrix2 = [
        l.strip().split(',') for l in
        open('../analysis/mmassiliense.RG_SC_RA.merge.fasta.csv', 'r')
    ]
    labels2 = matrix2[0][1:]
    matrix2 = [row[1:] for row in matrix2[1:]]
    print(
        len(matrix), len(matrix[0]), len(matrix2), len(matrix2[0]),
        labels, labels2
    )

    distances = []
    distances += get_distances(matrix, labels, 'abscessus')
    distances += get_distances(matrix2, labels2, 'massiliense')

    sns.set_palette("Reds")
    df = pd.DataFrame(distances, columns=['subspecies', 'type', 'distance'])
    sns.boxplot(x='subspecies', y='distance', hue='type', data=df)
    ax = plt.gca()
    ax.set_yscale('log')

    print(
        df.loc[df['subspecies'] == 'massiliense'].loc[df['type'] == 'within_pt']['distance'].median()
        
    )

    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_yticks([10, 100, 1000, 10000])

    plt.ylabel('SNV distance')
    plt.show()

if __name__ == '__main__':
    main()
