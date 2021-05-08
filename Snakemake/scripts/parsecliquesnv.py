"""


The purpose of this module is to refine the output produced by cliquesnv:
    * Replace the "fr" with patient id/timepoint data with the following field data (/ separator):
        - patient ID/timepoint/sample name/fraction (sig digits = 5)
    * Generate a table with the subsequent sequence IDs and patient IDs 
    * Add metadata in the form of timepoints

Still, we will allow the basic annotation of the sample based on an input value, i.e. timepoint metadata is not
enforced


"""



import sys,os
import argparse

import coordremap

import pandas as pd

import numpy as np
import argparse
import tempfile
import re


def extract_and_format(cliqueheader):
    
    haplo_num, frac = cliqueheader.split("_fr_")
    frac = float(frac)



    return haplo_num, frac


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i",dest="cliquefasta", type=str, help="Output FASTA from CliqueSNV/Other FASTA (provide the -p switch)")
    parser.add_argument("-n",dest="sample_name", type=str, default="infer", help="Name to replace 'fr' with")
    parser.add_argument("-o", dest="outfile", default="", type=str, help="Output file (default is STDOUT)")
    parser.add_argument("-m", dest="maskgap", default=15, type=int, help="Min length for gaps to be masked as 'N'") 

    args = parser.parse_args()

    fasta_files = coordremap.readfasta(args.cliquefasta)
    haplotypes = coordremap.readfasta(args.cliquefasta)
    

    reg_gap = re.compile("-{%s,}"%args.maskgap)
            


    out_fasta = ""
    for i, (header, seq) in enumerate(fasta_files):

        haplo_num, frac = extract_and_format(header)
        header = "%s/%s/%.4f"%(args.sample_name,haplo_num,frac)
        #print "Header:", header
        
        seq_arr = np.array(list(seq))
        for match in re.finditer(reg_gap, seq):


            start, end = match.start(), match.end()

            seq_arr[start:(end)] = ["N"]*(end-start)
        
        seq = "".join(seq_arr) 
        out_fasta += ">%s\n%s\n"%(header, seq)

    if args.outfile == "":
        print(out_fasta)
    else:
        with open(args.outfile,"w") as f:
            f.write(out_fasta)

    



    
    








