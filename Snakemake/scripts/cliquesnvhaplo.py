"""

This module assumes that the assembly coordinates are in accordance with the reference


"""


import pandas as pd
import numpy as np
import os
import subprocess

import pysam
import argparse



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", dest="bam", help="BAM file of the assembly")
    parser.add_argument("-g", dest="genelist", help="Genelist file to restrict haplotypes")
    parser.add_argument("-p", dest="output_prefix", default="", help="Output prefix (NOT output folder)")
    parser.add_argument("-o", dest="output_folder", help="Output folder")


    args = parser.parse_args()

    
    genelist = pd.read_csv(args.genelist)

    reads = [0] * len(genelist.gene)
    passed = [False] * len(genelist.gene)

    reads_base = [0.] * len(genelist.gene)
    coords_df = genelist.copy()
    for i, row in genelist.iterrows():
        
        gene = row.gene
        bamfile = pysam.AlignmentFile(args.bam,mode="r")

        reference = bamfile.references[0]
        #outfile = "{prefix}.bam"
        output_filename = os.path.join(args.output_folder, "%s.bam"%(gene))
        command_str =  'samtools view -hb {bamfile} "{reference}:{start}-{end}" -o {outfile} && samtools index {outfile}'.format(reference=reference,
                                                                                                       bamfile=args.bam,
                                                                                                       start=row.start,
                                                                                                       end=row.end,
                                                                                                       outfile=output_filename)

        print("Command: \n%s"%command_str)
        

        # Check bam file - if reads < 100 - fail
        command_out = subprocess.getoutput(command_str)


        command_bamcheck = 'samtools view {outfile} | grep -v "^@" | wc -l'.format(outfile=output_filename)

        num_reads = int(subprocess.getoutput(command_bamcheck))
        print("Num reads:", num_reads)
        reads[i] = num_reads
        if os.path.exists(output_filename):
            passed[i] = True
        reads_base[i] = round(num_reads/float(row.end-row.start+1))

    print("Reads:",reads)
    coords_df["reads"] = reads
    coords_df["reads_base"] = reads_base
    coords_df["pass"] = passed
    coords_df.to_csv(os.path.join(args.output_folder,"{prefix}_haplotype_info.txt".format(prefix=args.output_prefix)))






