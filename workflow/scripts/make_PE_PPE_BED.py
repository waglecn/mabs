#!/usr/bin/env python3
'''
This is a helper script to identify the PE and PPE CDS sequences annoated in
reference sequences from RefSeq and to store these regions in a BED file for
excluding variants from these regions.

Input:
    FILE: A GenBank format file with CDS having db_xref annotations
Output:
    STDOUT: A BED format file regions corresponding to CDS with PE/PPE domains

PE/PPE are highly conserved motifs in the N-term of proteins associated with
secretion systems in Mycobacteria that are thought to be involved in virulence
and/or immune invastion.

In MTB there are dozens of copies, highly homologous. The idea is that they
can interfere with accurate cSNP and hSNP determination.

It is unclear if these sequences would be similarly problematic in
Mycobacteriodes.

The PFAM/InterPro models for these conserved N-terminal domains identified in
annotated CDS need to be used to define regions to be possibly excluded from
variant calling.

Limitations:
- this presumes the GenBank RefSeq files are annotated properly
- that they have CDS sequences and have been scanned by PFAM/InterPro
- this also presumes that no pseudogenes or other features sequence homology

External sites:
<https://pfam.xfam.org/search/keyword?query=Mycobacterium+PE+PPE>
PPE <https://pfam.xfam.org/family/PF00823>
PE <https://pfam.xfam.org/family/PF00934>

References
Gey van Pittius et al (2006) BMC Evol Biol 6:95.
<https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-6-95>

Brennen (2017) Infect Immun 85(6):e00969-16
<https://iai.asm.org/content/85/6/e00969-16>
'''

import sys
from Bio import SeqIO

infile = sys.argv[1]

records = [r for r in SeqIO.parse(infile, format='gb')]
for r in records:
    CDS = [f for f in r.features if f.type == 'CDS']
    for c in CDS:
        label = []
        if 'db_xref' in c.qualifiers:
            for db_xref in c.qualifiers['db_xref']:
                if 'IPR000030' in db_xref:
                    label.append('PPE_IPR000030')
                if 'IPR000084' in db_xref:
                    label.append('PE_IPR000084')
        if len(label) > 0:
            print('{}\t{}\t{}\t{}'.format(
                r.id, c.location.start - 1, c.location.end, ','.join(label)
            ))
