# SARS-CoV-2 haplotype reconstructor 

This repo contains a Snakemake file that will take in a bunch of BAM files,
impute and extract putative hapoltypes with CliqueSNV and run the result through Virulign. 

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

I've included some tools necessary for the pipeline (under Snakemake/tools).  The CliqueSNV version is 1.5.6, BUT I modified a bit of the source - nothing in the algorithm.  Please only use the version on the repo.


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
	- Genes with coverage < 50% will be removed from subsequent analysis.  This "fixes" the CliqueSNV error where the gene coverage is too sparse
2. Downsample the BAMs by _depth_.  Rather than using a proportional downsample as one can do with `samtools` and others, there is rather an upper limit of the depth of the BAMs.  This downsampling/depth clipping expidites the haplotype calling process.
3. Haplotype imputation with CliqueSNV.  CliqueSNV is a rather new tool and as of writing still in the process of peer-review.  Nevertheless, it has been independently evaluated. 
	- Side note: Any gene with a too-small coverage will be excluded from subsequent steps
4. Alignment and amino acid extraction with Virulign.  Virulign is used using reference XMLs in the `assets/virulign` folder to do a proper alignment against the reference.  Virulign uses a method similar to AGA in that during the sequence alignment, codons are also take into account, with the output containing the miminmal amount of frameshifts and a high quality alignment  (Note: Currently, we only export the AA sequence table, but will probably include nucleotide/protein alignments in the near future)
5. Mutations are tallied and agglomerated into a single CSV file


## Notes:

For any alignment containing an unknown AA (either by lack of coverage/frameshift error in the alignment), will be masked with X.  This is because virulign performs a global alignment, but strips the gaps beforehand. For the sake of being conservative, any stretch of deletions longer than 12 (by default) will be masked with Ns and as a consequence, the reference AAs in those regions will be used. This could be changed, of course, to longer stretches.  CliqueSNV reports 0 coverage regions as "gap" characters too.  So stretches of no coverage would explicitly be masked by Ns.  








