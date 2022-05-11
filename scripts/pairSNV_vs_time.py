#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# defaults
if len(sys.argv[1:]) == 0:
    # pairwise SNV data
    in_csv = '/media/nick/2TB/backup-ignore/MA/results/gubbins/mabs.test.dist.csv'
    # sample data with dates
    in_samples = '/home/nick/Dropbox/MA/sample.csv'
else:
    in_csv = sys.argv[1]
    in_samples = sys.argv[2]

df1 = pd.read_csv(in_csv, index_col=0)
df2 = pd.read_csv(in_samples, index_col=0)
df2["Sample_date"] = pd.to_datetime(df2['Sample_date'])

data = []

all_pts = {}
for i in df1:
    if '-' in i:
        pt = int(i.split('-')[0])
        if pt not in all_pts:
            all_pts[pt] = []
        all_pts[pt].append(i)
    else:
        print(i, file=sys.stderr)

filtered_pts = {pt: hits for (pt, hits) in all_pts.items() if len(hits) > 2}


for pt, samples in filtered_pts.items():
    sample_pairs = []
    included_samples = df1[df1.index.str.startswith('{}-'.format(pt))].index
    pt_samples = df2.loc[included_samples]['Sample_date']
    earliest = pt_samples[pt_samples == pt_samples.min()].index[0]
    other_samples = [i for i in samples if i is  not earliest]
    sample_pairs = [(earliest, i) for i in other_samples]

    if len(sample_pairs) < 3:
        continue
    for pair in sample_pairs:
        data.append(
            (
                str(pt),
                abs(
                    (
                        df2.loc[pair[1]]['Sample_date'] -
                        df2.loc[pair[0]]['Sample_date']
                    ).days
                ),
                df1.loc[pair[1], pair[0]],
                pair[0], df2.loc[pair[0]]['Sample_date'],
                pair[1], df2.loc[pair[1]]['Sample_date']
            )
        )
cols = [
    'patient', 'days', 'dist', 'pt0', 'sample_0_date', 'pt1', 'sample_1_date'
]

data = pd.DataFrame(data, columns=cols)
print(data)
data = data[data['patient'] != '17']

# g = sns.FacetGrid(data, col="patient", hue='patient')
#g.map(sns.lmplot, x="days", y="dist")
g = sns.lmplot(data=data, col='patient', col_wrap=4, x="days", y="dist")
g.set(ylim=(0, 100))
g.set(xlim=(0, 2000))
plt.ylabel('Pairwise SNVs')
plt.xlabel('infection days')
plt.show()

