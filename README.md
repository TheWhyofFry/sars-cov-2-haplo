# SARS-CoV-2 haplotype reconstructor 

*Major update* - consider this if you cloned the repo prior to 07 July 2021.

This repo contains a Snakemake file that will take in a bunch of BAM files,
impute and extract putative hapoltypes with CliqueSNV and run the result through Virulign. 

Please note that it assumes there is 100% concordance between the positions in your BAM
with those in the [reference SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). Having said that, the pipeline will take a differently named reference in the BAM into consideration.


There is nothing to cite and you can use/change it however you like but any critique/concern/bug can raised by opening an issue. This includes asking for features.

_*Basic as the code is, I do expect people to not claim it as their own.*_

If you want to drop me a line, send and email to werner.smidt@gmail.com

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

3. The pipeline has only been tested on Illumina data

I've included some tools necessary for the pipeline (under Snakemake/tools), except CliqueSNV.  The version the pipeline was tested on, is CliqueSNV 1.5.7. CliqueSNV is under active development, so if something fails with a new version, please open an issue/try the aforementioned version.  So just get the JAR file for CliqueSNV and put it in `tools/`.  

Pretty much a WIP, so changes will be made in the future.

# Workflow overview

This pipeline makes use of various existing tools with intermediate scripts that takes in a set of reference aligned BAM files.

## Pipeline process:

1. Split BAMs into individual gene bams
2. Downsample the BAMs by _depth_.  Rather than using a proportional downsample as one can do with `samtools` and others, there is rather an upper limit of the depth of the BAMs.  This downsampling/depth clipping expidites the haplotype calling process.  The default is a depth of 3000. This can be changed, but at very high depths, CliqueSNV will use a high amount of RAM. 
3. Run the `lofreq` Viterbi algorithm to realign potentially misaligned reads (gives a higher confidence in INDEL calling). In future, this may be replaced with [ABRA2](https://github.com/mozack/abra2).
4. Haplotype imputation with CliqueSNV.  CliqueSNV is a rather new tool and as of writing still in the process of peer-review.  Nevertheless, it has been independently evaluated. 
5. Translation and alignment to reference genes with a script that uses a combination of Virulign/AGA.  
6. Mutations are tallied and agglomerated into a single CSV file


# Limitations

- Currently, the pipeline assumes that the input BAMs are from reads aligned to the mentioned SARS-CoV-2 reference genome. 
- No whole genome haplotypes.  I might include this in the future. 
- By default, I assume the input data is from Illumina origin


# Short term features 

- Adding some quality information to the imputed SNPs in the haplotypes, e.g. depth at the position
- Use Salmon to re-estimate the haplotype frequencies
- Add scripts to extract and align all the output sample gene haplotypes


# Running example 


Create a folder called "inputdir" in the path containing the Snakemake file and place your BAMs within this folder.  Run the following command

`snakemake -j 4 -p --config inputdir=inputdir -k`

By default, the sample names are imputed from the BAM files. I will change this behaviour in the near future so you can include a custom sample sheet to manually name the samples.

You can also subset to a gene/list of genes: 

`snakemake -j 4 -p --config inputdir=inputdir genes=S,ORF1b`

Will subset the analysis to the listed genes - _case sensitive_. You can rerun the pipeline by adding a gene, e.g. `genes=S,ORF1b,ORF1a` and it _should_ just run the bits needed to create the haplotypes for `ORF1a`. *Importantly*, you should delete the folder `output/mutation_tables_agg/` and the file `output/mutation_table_aggregate.csv'.


# Output

Currently, all outputs are contained in a the `output` folder of the Snakemake folder.  

- `aatables` - Contains a list of CSVs.  The format of the name is `{sample}_{gene}.csv`.  This file shows the amino acid residues found for each haplotype as well as the reference AAs.
- `agg_aatables` - Same as `aatables`, but the proportion identical AAs at each position over all haplotypes are aggregated
- `haplotype` - All the reconstructed haplotypes per gene. Separated into folders based on sample name.  
- `virulign` - Virulign AA output for all genes and all samples.  Each subfolder is a gene and csv named according to the sample name
- `mutation_tables` - A list of all the mutations different from the reference set.  Name format is `{sample}_{gene}.csv`









