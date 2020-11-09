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

# Notes
200915
- strange bug causing infinite loop in snakemake downloading refseq genomes. I think this is because of the dynamic() output/input in rules. Checking this out, seeing if the bug happens if I run entire pipeline from scratch.

200917
- noticed a bug in running shovill, increased expected memory usage. Shovill version 0.9.0 running from an older miniconda. Removed miniconda, started from scratch, and pinned Shovill 1.1.0 in shovill.yaml
- after fixing, rerunning seems to work with example data, then works after deleting the mashtree and refseq_download directories. 



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
		- [X] download kraken db automatically - deferred, added to Required files
	- [X] genome assembly on raw reads - 200914
		- [X] Erm(41) identification on assembly - 200912
	- [X] kraken2 on assembly - 200912
	- [X] mashtree assembly - 200913
	- [X] map everything to ATCC 19977 for basic coverage - 200914

- [ ] build workflow part 2 on available assemblies
	
	- [X] tree-guided MRCA - 200915
	- [X] MRCA-guided MLST - 200913
	- [X] MRCA-guided reference mapping - 200921
	- [ ] Optional: Mark duplicates with picard
	- [X] read filtering - see Martin et al 2018 and Lee et al 2020
		- [X] filter soft clips - 200922
		- [X] optional GATK realignment, but see for why it was removed in 2015 for gatk4 https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847
			- [X] added 200923, optional 200924
			- intially added gatk4, got errors and followed the rabbit-hole
			- to follow Martin et al, added conda env with gatk3.8, since the resulting bam can be used with any downstream variant caller
	- [ ] annotate regions of interest
		- remove PP/PPE regions (BED file)
			- [X] identify PP/PPE - 200927
		- [X] zero coverage of reference
		- [ ] remove phage, tnp, IS
		- [X] merge ROI BED files
	- [X] MRCA-guided variant calling with bcftools - 200922
		- [X] bcftools mpileup - 200923
		- [X] called variants - 200923
		- [X] variant filtering
			- [X] basic Martin et al - 200925
			- [ ] density filter - see <https://github.com/c2-d2/within-host-diversity/blob/master/fastq_to_vcf_pipeline.py#L122> line 
	- [X] variant annotation with SNPEff
	- [X] SNP-tree construction
		- [X] SNP extraction - custom? merge vcf as per Robyn 201006
		- [X] - merge SNPs - 201013
		- [X] concatenate cSNPSs (exclude hSNPs) 201016
			- snp-sites ? snippy?
		- [X] - vcfmerge 201014

# TODO 200911
- [X] add trimming parameters to config file - 200921

# TODO 200914
- sub-species type assemblies are hard-coded in scripts/tree_MRCA.py, it would be useful for these to be configurable but adds layers of complexity to snakefile

# TODO 200920
- Added GATK info to REQUIREMENTS, and config.yaml

# TODO 200926
- [ ] Tune variant filtering
- [X] TODO big question here - use stats from part 1 to make *new* sample_sheet with QC pass samples? No
	- [X] make list to prune from SNP alignment - not needed 201012
- [X] need separate list of in-complete genomes, as MRCA-guided MLST didn't work as expected, tree has wrong structure (samples from pt 29 should be mmas) - Fixed 201006, need to convert gbff files before mashtree can read

# TODO 201010
- [X] start density filter
- [X] merge completed results without recalculating shovill assemblies for old samples - 201010
- [X] merge 0-coverage bed files and PE_PPE bed files 201013
- [X] filter merged bed from vcf
	- [X] compress vcf with bcftools

# TODO 201013
- [ ] complete density filter

# TODO 201015
- [X] incorporate https://github.com/phac-nml/mab_mabscessus 211021