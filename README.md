# mabs

author:nwaglechner@gmail.com

# Basic Setup

```bash
git clone https://github.com/waglecn/mabs.git
```

Conda and snakemake

Miniconda available from:
https://docs.conda.io/en/latest/miniconda.html

Python 3.8.3 Miniconda
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  
bash Miniconda3-latest-Linux-X86_64.sh
conda env create --name mabs --file environment.yaml
conda activate mabs
```
	- note the version of python installed in the the mabs environment is not necessarily the same as the default miniconda python version
	- asking for ete3 in the default environment will required python 3.6 (200921)

## Required files:
- adapters for trimming - see config.yaml
	- look for adapter files bundled with trimmomatic, ie.
```bash
locate TruSeq3-PE.fa
```
- Kraken database 
(preferred) https://benlangmead.github.io/aws-indexes/k2
(original source) ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20231009.tar.gz
```

# How to run

The workflow proceeds in stages:
- stage1 - QC and processing each sample, generation of stage2 configuration
- stage2 - mapping samples against references in stage2.yaml (see workflow/stage2.smk). Right now, you need to generate this file in the results directory yourself.
- stage3 - checking for recombination, pairwise snv distances, and phlyogenetics

```
snakemake --configfile config/config.yaml --cores 32 --use-conda --conda-frontend mamba -k [ stage1 | stage2 | stage3 ]
```

