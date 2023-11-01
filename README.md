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
- GATK3 jar file
	- available from <https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk>
	- used '''GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2'''
	- https://storage.cloud.google.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
	- see config.yaml
- adapters for trimming - see config.yaml
	- look for adapter files bundled with trimmomatic, ie.
```bash
locate TruSeq3-PE.fa
```
- Kraken database 
ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/
```bash
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
```

# How to run
```
snakemake --configfile config.yaml --cores 8 --use-conda --conda-prefix /path/to/.snakemake/conda
```

