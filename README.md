# SARS-CoV-2 haplotype reconstructor 



This repo contains a Snakemake file that will take in a bunch of BAM files,
impute and extract putative hapoltypes with CliqueSNV and run the result through Virulign. 

Please note that it assumes there is 100% concordance between the positions in your BAM
with those in the [reference SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512).

There is nothing to cite and you can use/change it however you like but any critique/concern/bug can raised by opening an issue. 

_*Basic as the code is, I do expect people to not claim it as their own.*_


# Requirements

The following modules are needed:

1. Python:
  - snakemake
  - pysam
  - pandas

2. General:
  - samtools
  - [virulign](https://github.com/rega-cev/virulign) 
  - [cliquesnv](https://github.com/vtsyvina/CliqueSNV/archive/refs/tags/1.5.7.tar.gz)
  - [AGA](https://github.com/emweb/aga)

I've included some tools necessary for the pipeline (under Snakemake/tools), except CliqueSNV.  The version the pipeline was tested on, is CliqueSNV 1.5.7. CliqueSNV is under active development, so if something fails with a new version, please open an issue/try the aforementioned version.


# Running example 


Create a folder called "inputdir" in the path containing the Snakemake file and place your BAMs within this folder.  Run the following command

`snakemake -j 4 -p --config inputdir=inputdir -k`

You can also subset to a gene/list of genes 

`snakemake -j 4 -p --config inputdir=inputdir genes=S,ORF1b`

Will subset the analysis to the listed genes - _case sensitive_.


# Output

Currently, all outputs are contained in a the `output` folder of the Snakemake folder.  

- `aatables` - Contains a list of CSVs.  The format of the name is `{sample}_{gene}.csv`.  This file shows the amino acid residues found for each haplotype as well as the reference AAs.
- `agg_aatables` - Same as `aatables`, but the proportion identical AAs at each position over all haplotypes are aggregated
- `haplotype` - All the reconstructed haplotypes per gene. Separated into folders based on sample name.  
- `virulign` - Virulign AA output for all genes and all samples.  Each subfolder is a gene and csv named according to the sample name
- `mutation_tables` - A list of all the mutations different from the reference set.  Name format is `{sample}_{gene}.csv`



# Workflow overview

This pipeline makes use of various existing tools with intermediate scripts that takes in a set of reference aligned BAM files.

## Pipeline process:

1. Split BAMs into individual gene bams
2. Downsample the BAMs by _depth_.  Rather than using a proportional downsample as one can do with `samtools` and others, there is rather an upper limit of the depth of the BAMs.  This downsampling/depth clipping expidites the haplotype calling process.
3. Haplotype imputation with CliqueSNV.  CliqueSNV is a rather new tool and as of writing still in the process of peer-review.  Nevertheless, it has been independently evaluated. 
4. Alignment and amino acid extraction with Virulign.  Virulign is used using reference XMLs in the `assets/virulign` folder to do a proper alignment against the reference.  Virulign uses a method similar to AGA in that during the sequence alignment, codons are also take into account, with the output containing the miminmal amount of frameshifts and a high quality alignment
	- Sometimes Virulign fails.  In the event this happens, AGA is used to produce the appropriate alignment
5. Mutations are tallied and agglomerated into a single CSV file


## Notes:

For any alignment containing an unknown AA (either by lack of coverage/frameshift error in the alignment), will be masked with X.  This is because virulign performs a global alignment, but strips the gaps beforehand. Newer versions of CliqueSNV marks sites with zero coverage with "N".






