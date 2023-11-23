#!/usr/bin/env python
import sys
import pandas as pd

def main():
    # argv1 ought to be the QC_summary.csv
    inh1 = sys.argv[1]
    # argv2 ought to be the path to the config file relative to snakefile
    inh2 = sys.argv[2]

    data = pd.read_csv(inh1)
    print(data['sample'], file=sys.stderr)

    print("""
mmassiliense:
 - Sample_2B
 - Sample_3B
 - Sample_4B
 - Sample_5B
    """)




if __name__ == '__main__':
    main()
