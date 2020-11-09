#!/usr/bin/env python3
"""
This script will return which of several marked sequences shares the MRCA with
a target sample in order to assign a lineage to a M abscessus sample

Example usage:
    tree_MRCA.py [input.tree] [target.sample]

INPUT:
    FILE: input.tree - a tree which includes the hardcoded reference sequences
                        from M. abscessus and the target.sample leaf
OUTPUT:
    STDOUT: MRCA - the name of the reference subspecies that the sample is
                    assigned to
"""

import sys
import os
import ete3


def tree_MRCA(in_tree, sample):
    if not os.path.exists(in_tree):
        print('no tree')
        return None

    mbtree = ete3.Tree(in_tree)
    leaves = mbtree.get_leaves()

    # need to assume something about the name
    sample_node = [l for l in leaves if l.name.startswith(sample)]
    # print(sample, 'leaves', sample_node)
    if len(sample_node) == 0:
        return None

    targets = {
        'mbol_target': None,
        'mmas_target': None,
        'mabs_target': None
    }
    for t in leaves:
        if 'GCF_003609715.1' in t.name:
            targets['mbol_target'] = t
        #if 'GCF_000239055.1' in t.name:
        if 'GCF_000497265.2' in t.name:
            targets['mmas_target'] = t
        if 'GCF_000069185.1' in t.name:
            targets['mabs_target'] = t

    mbol_dist = targets['mbol_target'].get_distance(sample_node[0])
    mmas_dist = targets['mmas_target'].get_distance(sample_node[0])
    mabs_dist = targets['mabs_target'].get_distance(sample_node[0])

    mindist = min([mbol_dist, mmas_dist, mabs_dist])
    if mindist == mbol_dist:
        return 'mbolettii'
    elif mindist == mmas_dist:
        return 'mmassiliense'
    elif mindist == mabs_dist:
        return 'mabscessus'

if __name__ == '__main__':
    in_tree = sys.argv[1]  # "../10-mashtree/MBtree"
    sample = sys.argv[2]  # 1-B
    target = tree_MRCA(in_tree, sample)
    # used in shell to determine MLST
    if target == "mbolettii":
        print('mmassiliense')
    else:
        print(target)
