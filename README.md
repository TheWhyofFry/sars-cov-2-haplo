# SARS-CoV-2 

This repo contains a Snakemake file that will take in a bunch of BAM files,
impute and extract putative hapoltypes and run the result through Virulign. 

Please note that it assumes there is 100% concordance between the positions in your BAM
with those in the [reference SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512).

# Requirements

The following modules are needed:

1. Python:
  - snakemake
  - pysam
  - pandas

2. General:
  - samtools
  - [virulign](https://github.com/rega-cev/virulign) 


# Running example 


Create a folder called "inputdir" and place your BAMs within this folder.  Run the following command

`snakemake -j 4 -p --config inputdir=inputdir` -k



