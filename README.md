# SARS-CoV-2 

This repo contains a Snakemake file that will take in a bunch of BAM files,
impute and extract putative hapoltypes and run the result through Virulign. 



# Requirements

The following modules are needed:

- pysam
- pandas



# Running example 


Create a folder called "inputdir" and place your BAMs within this folder.  Run the following command

`snakemake -j 4 -p --inputdir=inputdir`

*Currently, you will have to run it twice - there is something wrong with the logic*


