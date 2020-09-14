#!/usr/bin/env python3
import sys
import ete3

in_tree = sys.argv[1]  # "../10-mashtree/MBtree"
sample = sys.argv[2]  # 1-B

mbtree = ete3.Tree(in_tree)
leaves = mbtree.get_leaves()

# need to assume something about the name
sample_node = [l for l in leaves if sample in l.name][0]

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

mbol_dist = targets['mbol_target'].get_distance(sample_node)
mmas_dist = targets['mmas_target'].get_distance(sample_node)
mabs_dist = targets['mabs_target'].get_distance(sample_node)

mindist = min([mbol_dist, mmas_dist, mabs_dist])
if mindist == mbol_dist:
    print('mmassiliense')
elif mindist == mmas_dist:
    print('mmassiliense')
elif mindist == mabs_dist:
    print('mabscessus')
