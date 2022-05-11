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


def same_pt(s1, s2):
    pt1 = s1.split('-')[0]
    pt2 = s2.split('-')[0]
    if pt1 == pt2:
        return True
    return False


def main():
    mabs = pd.read_csv("results/SNP_phylo/mabs.csv", index_col=0)
    mmas = pd.read_csv("results/SNP_phylo/mmas.csv", index_col=0)

    distances = pd.DataFrame(columns=[
        'subspecies', 'A', 'B', 'type', 'distance']
    )
    mabs_pts = mabs.index.values
    for i, x in enumerate(mabs_pts):
        for y in mabs_pts[i + 1:]:
            if same_pt(x, y):
                t = 'within_pt'
            else:
                t = 'between_pt'
            distances = distances.append({
                'subspecies': 'mabscessus',
                'A': x,
                'B': y,
                'type': t,
                'distance': mabs[x][y]
            }, ignore_index=True)
    mmas_pts = mmas.index.values
    for i, x in enumerate(mmas_pts):
        for y in mmas_pts[i + 1:]:
            if same_pt(x, y):
                t = 'within_pt'
            else:
                t = 'between_pt'
            distances = distances.append({
                'subspecies': 'mmassiliense',
                'A': x,
                'B': y,
                'type': t,
                'distance': mmas[x][y]
            }, ignore_index=True)
    mmas_pts = mmas.index.values

    for i in ['mabscessus', 'mmassiliense']:
        for j in ['between_pt', 'within_pt']:
            print(i, j)    
            print(
                distances[
                    (distances['subspecies'] == i) &
                    (distances['type'] == j)
                ]['distance'].quantile([0.25, 0.5, 0.75])
            )
    
    
    sns.set_palette("Reds")
    sns.boxplot(x='subspecies', y='distance', hue='type', data=distances)
    ax = plt.gca()
    ax.set_yscale('log')

    print(
        distances.loc[distances['subspecies'] == 'mabscessus'].loc[distances['type'] == 'within_pt']['distance'].median(),
        distances.loc[distances['subspecies'] == 'mmassiliense'].loc[distances['type'] == 'within_pt']['distance'].median()
    )

    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.yaxis.set_major_formatter(formatter)
    ax.set_yticks([10, 100, 1000, 10000])
    xticks = ax.get_xticks()
    print(xticks)
    plt.xticks(xticks, ['$\it{M. abscesssus}$', '$\it{M. massiliense}$'])

    plt.ylabel('SNVs')
    plt.title('subspecies SNV distance')
    plt.show()

if __name__ == '__main__':
    main()
