# mabs

author:nwaglechner@gmail.com

# Basic Setup

```bash
git clone https://github.com/waglecn/mabs.git
```

Conda and snakemake

Miniconda available from:
https://docs.conda.io/en/latest/miniconda.html

Python 3.7.7 Miniconda  
```bash
wget https://repo.anaconda.com/miniconda/
```

Miniconda3-latest-Linux-x86_64.sh  
```bash
bash Miniconda3-latest-Linux-X86_64.sh
conda env create --name mabs --file environment.yaml
conda activate mabs
```

## Required files:
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


# TODO 200902

- [ ]download internal project data - deferred
	- [X] configurable data-dir - 200914
- download external project data
	- [X] refseq genomes - done 200904
	- [ ] genomes from Bryant et al, SRA
		- need to know what these are
- [X] download reference assemblies - 200908
	- first used all contig assemblies, changed to 'complete' keyword

- reading in samples somehow, obviously this depends on how/where they are downloaded (see previous TODO item) and the data that is already downloaded
	- need a dummy rule that requires these as input in order to define wildcards

- [X] basic Snakefile - 200905

- [X] build workflow part 1
	- [X] index reference assemblies - deferred 200914
		- moved to resources/alignment_references
	- [X] pre-trim QC - done 200908
	- [X] trim - done 200909
		- specify adapter files, add variable to config
	- [X] post-trim QC done 200909
	- [X] kraken check - done 200910
		- [ ] download kraken db automatically - deferred, added to Required files
	- [X] genome assembly on raw reads - 200914
		- [X] Erm(41) identification on assembly - 200912
	- [X] kraken2 on assembly - 200912
	- [X] mashtree assembly - 200913
	- [X] map everything to ATCC 19977 for basic coverage - 200914

- [ ] build workflow part 2 on available assemblies
	- [ ] TODO big question here - use stats from part 1 to make *new* sample_sheet with QC pass samples?
	- [ ] tree-guided MRCA -
	- [X] MRCA-guided MLST - 200913
	- [ ] MRCA-guided reference mapping
	- [ ] MRCA-guided variant calling
		- read filtering
		- site masking ?
		- called variants
	- [ ] SNP-tree construction

# TODO 200911
- add trimming parameters to config file

# TODO 200914
- sub-species type assemblies are hard-coded in scripts/tree_MRCA.py, it would be useful for these to be configurable but adds layers of complexity to snakefile
