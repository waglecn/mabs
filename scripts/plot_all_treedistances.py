#!/usr/bin/env python3
"""
This script will produce a basic histogram of distances extraced from a tree

Example usage:
    plot_all_treedistances.py [input.tree]

INPUT:
    FILE: input.tree - a newick format tree
OUTPUT:
    STDOUT: matplotlib histogram of pairwise distances extracted from the input
"""

import sys
import ete3
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from plot_phylo_dist_heatmap import make_matrix


def main():
    tree_file = sys.argv[1]
    tree = ete3.Tree(tree_file)
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize(direction=1)
    matrix = make_matrix(tree)

    for node in tree:
        if node.name == 'ref':
            node.name = os.path.split(tree_file)[-1].split('.')[0]

    labels = [l.name for l in tree.get_leaves()]
    matrix = pd.DataFrame(matrix, columns=labels, index=labels)
    matrix.to_csv('{}.distances.csv'.format(tree_file))

    keep = np.triu(np.ones(matrix.shape)).astype('bool').reshape(matrix.size)
    distances = matrix.stack()[keep]
    sns.histplot(distances)
    plt.show()

if __name__ == '__main__':
    main()
