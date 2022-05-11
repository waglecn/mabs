#!/usr/bin/env python3

import pandas as pd

data = pd.read_csv(
    '~/Dropbox/MA/SRA_metadata/all.csv',
    low_memory=False
)

samples = data
print(samples, "{:.3f} GB total".format(samples['Bytes'].sum() / (1024 ** 3)))
# print(samples['Bytes'].summary())


rule all:
    input:
        expand("{run}.done", run=samples['Run'])


rule download:
    threads: 1
    output:
        "{run}.done",
    shell:
        # the -s flag is necessary to properly split files
        # see ERR4339027, needed to be rerun
        "prefetch -C yes -f yes -p -X 100000000 {wildcards.run} ;"
        "fasterq-dump -p -e 1 -f -S -s {wildcards.run} -O /media/nick/8TB/MA/ ; "
        "rm -Ifr ./{wildcards.run} ;"
        "touch {output} "

