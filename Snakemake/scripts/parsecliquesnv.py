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

def extract_and_format(cliqueheader):
    
    haplo_num, frac = cliqueheader.split("_fr_")
    frac = float(frac)



    return haplo_num, frac

# Gets low/no coverage regions from a referece and replaces all "gaps" with Ns
def refmask(alignment_file, reference_file):
     
    temp_alignment_file = tempfile.mktemp(dir="/tmp")


    command = "mafft --quiet --keeplength --add {reference} {alignment} > {tempfile}".format(**dict(reference=reference_file, alignment=alignment_file,tempfile=temp_alignment_file))
    exit_code = os.system(command)



    alignment = coordremap.readfasta(temp_alignment_file,upper=True)[::-1]
    mat = np.array(list(map(list, [seq for title,seq in alignment])))
    N_instances = mat[0,] == "N"

    mat[:,N_instances] = "N"
    mat = mat[1:,]
    alignment = alignment[1:][::-1]
    
    #mat_cons = np.apply_along_axis(lambda x:"".join(set(x)),0,mat)
    #Delete gap-only
    #mat_filter = mat_cons != "-"


    #mat = mat[:,mat_filter]
    mat = list(np.apply_along_axis(lambda x:"".join(x), 1, mat))

    alignment = [(title,m) for (title,seq),m in zip(alignment,mat)]

    os.unlink(temp_alignment_file)





    return alignment


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i",dest="cliquefasta", type=str, help="Output FASTA from CliqueSNV/Other FASTA (provide the -p switch)")
    parser.add_argument("-n",dest="sample_name", type=str, default="infer", help="Name to replace 'fr' with")
    parser.add_argument("-p", dest="pacbio",action="store_true", default=False, help="Is the sequence pacbio?")
    parser.add_argument("-c",dest="id_column", default="samplename",type=str, help="Column identifying sample")
    parser.add_argument("-r", dest="reference", default="", type=str, help="Reference file for filling in Ns")
    parser.add_argument("-o", dest="outfile", default="", type=str, help="Output file (default is STDOUT)")
    

    args = parser.parse_args()

    haplotypes = coordremap.readfasta(args.cliquefasta)
    

    if args.reference == "":
        fasta_files = coordremap.readfasta(args.cliquefasta)
    else:
        if os.path.exists(args.reference):
            fasta_files = refmask(args.cliquefasta,args.reference)
        else:
            sys.stderr.write("WARNING: Reference %s not found"%(args.reference))
            fasta_files = coordremap.readfasta(args.cliquefasta)
            


    out_fasta = ""
    for i, (header, seq) in enumerate(fasta_files):

        if not args.pacbio:
            haplo_num, frac = extract_and_format(header)
        else:
            haplo_num, frac = i+1, 1
        seqtype = "shortread" if not args.pacbio else "pacbio"
        header = "%s/%s/%s/%.3f"%(args.sample_name,haplo_num,seqtype,frac)
        #print "Header:", header
        out_fasta += ">%s\n%s\n"%(header, seq)

    if args.outfile == "":
        print(out_fasta)
    else:
        with open(args.outfile,"w") as f:
            f.write(out_fasta)

    



    
    








