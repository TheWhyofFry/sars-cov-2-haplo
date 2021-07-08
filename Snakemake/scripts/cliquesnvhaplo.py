"""

This module assumes that the assembly coordinates are in accordance with the reference


"""


import pandas as pd
import numpy as np
import os
import subprocess

from io import StringIO

import pysam
import argparse

import sys

def get_gene_stats(bam, genelist, ref):

    rows = []
    for i, row in genelist.iterrows():

        reads = bam.count(contig=ref, start=row.start, end=row.end)
        reads_base = int(reads/(row.end-row.start+1))

        rows.append((row.gene, reads, reads_base))

    return pd.DataFrame.from_records(rows, columns=["gene","reads","reads_base"])





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", dest="bam", help="BAM file of the assembly")
    parser.add_argument("-g", dest="genelist", help="Genelist file to restrict haplotypes")
    parser.add_argument("-G", dest="only_gene", default="", help="Restrict to this comma-separated gene list")
    parser.add_argument("-p", dest="output_prefix", default="", help="Output prefix (NOT output folder)")
    parser.add_argument("-o", dest="output_folder", help="Output folder")
    parser.add_argument("-r", dest="ref", default="", help="Reference FASTA or reference name")
    parser.add_argument("-s", dest="stats_only", action="store_true", help="Only generate stats")
    parser.add_argument("-n", dest="no_stats_output", action="store_true", help="Do not output stats")
    args = parser.parse_args()

    
    if args.ref != "":
        if os.path.exists(args.ref):
            with open(args.ref, "r") as fasta_file:
                for line in fasta_file:
                    if line.startswith(">"):
                        line = line.strip()[1:]
                        ref = line
                        break
            ref = None
    else:
        ref = None

    genelist = pd.read_csv(args.genelist)

    if args.only_gene != "":
        genelist = genelist[genelist.gene.isin(args.only_gene.split(","))]

    bamfile = pysam.AlignmentFile(args.bam,mode="rb")




    coords_df = get_gene_stats(bamfile, genelist, ref)
   
    if not args.no_stats_output:

        coords_df.to_csv(os.path.join(args.output_folder,"{prefix}_haplotype_info.txt".format(prefix=args.output_prefix)))

    if args.stats_only:
        sys.exit(0)

    with StringIO(pysam.depth(args.bam)) as sio:
            depth_df = pd.read_table(sio, header=None, names=["ref","pos","depth"], sep="\t")
    
    for i, row in genelist.iterrows():
        
        gene = row.gene


        depth_pass = len(depth_df[(depth_df.pos >= row.start) & (depth_df.pos <= row.end) & (depth_df.depth > 0)])

        gene_length = (row.end - row.start) + 1

        depth_frac = depth_pass / float(gene_length)


        
        reference = bamfile.references[0] if ref is None else ref



        #outfile = "{prefix}.bam"
        output_filename = os.path.join(args.output_folder, "%s.bam"%(gene))
        command_str =  'samtools view -hb {bamfile} "{reference}:{start}-{end}" -o {outfile} && samtools index {outfile}'.format(reference=reference,
                                                                                                       bamfile=args.bam,
                                                                                                       start=row.start,
                                                                                                       end=row.end,
                                                                                                       outfile=output_filename)

        command_out = subprocess.getoutput(command_str)

        






