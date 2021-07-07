"""

The purpose of this script is to run virulign on input sequences and
deal with any pesky FrameShift failures that occur.

Procedure

1. Run Virulign on query sequences using max of 0 FrameShifts
2. If some make it through
   - recreate a reference from those sequences (assuming they do not start/end with gaps)
   - run virulign on these references 
3. Failed sequences are run through AGA
4. May improve this in future - create references from blast hits and realign to HXB2 after the fact


virulign ref.xml/fasta target --exportWithReference yes --exportAlphabet Nucleotides/Aminoacids --exportKind GlobalAlignment


"""

import os
import shutil
import subprocess
import io
import tempfile
import pandas as pd
import math
from runcommand import run_output
import agagen
import re
import pathlib

import psutil

proc = psutil.Process()


def remap(seq):

    seq_nogaps = seq.replace("-","")

    return pd.Series([i for i,s in enumerate(seq) if s != '-'],index=range(len(seq_nogaps)))






def read_fasta(f,upper=True):
    #fasta_dict = {}
    fasta_list = []
    curseq = ""
    currentkey = ""
    fastafile = open(f, "r") if type(f) is str else f
    for line in fastafile:
        if line.startswith(">"):
            if curseq != "":
                #fasta_dict[currentkey] = curseq
                try:
                    fasta_list.append((currentkey,curseq.upper()))
                    currentkey = line[1:].strip()
                    curseq = ""
                except:
                    print(f)
                    print(curseq)
                    raise
            else:
                
                currentkey = line[1:].strip()
        else:
            curseq += line.strip().replace(".","-")

    fastafile.close()


    fasta_list.append((currentkey, curseq))
    return fasta_list


def write_fasta(fasta_list, output_file="string"):

    output = "\n".join([">%s\n%s\n"%(title,seq) for title, seq in fasta_list])

    if output_file == "string":
        return output
    else:
        with open(output_file,"w") as outfile:
            outfile.write(output)
        return True


def parse_virulign(stdout):
    
    with io.StringIO(stdout) as virulign_output:
        fasta_entries = read_fasta(virulign_output)

    return fasta_entries

 
# Fasta is a fasta list produced by read_fasta
def find_missing_entries(query_fasta, virulign_fasta):

    query_entries = [title for title, seq in query_fasta]
    virulign_entries = [title for title, seq in virulign_fasta]

    
    missing_entries =  set(query_entries).difference(virulign_entries)
    return [(title, seq) for title,seq in query_fasta if title in missing_entries]



def initial_alignment(fasta_file, virulign_ref_file, virulign_command="virulignmp",
                      alphabet="AminoAcids", exportkind="GlobalAlignment",
                      exportwithref="yes", maxFrameShifts=3):
    
    virulign_command_template = "{v} {virulign_ref_file} {fasta_file} --exportKind {exportkind} --exportAlphabet {alphabet} --exportReferenceSequence {ref}" + \
                                   " --maxFrameShifts {maxFrameShifts}"

        
    stdout, stderr, returncode = run_output(virulign_command_template.format(v=virulign_command,
                                                                     virulign_ref_file=virulign_ref_file,
                                                                     fasta_file=fasta_file,
                                                                     exportkind=exportkind,
                                                                     alphabet=alphabet,
                                                                     ref=exportwithref,
                                                                     maxFrameShifts=maxFrameShifts))
    # Fix that quirky case when arbitrary characters are written out by virulignmp
    if ">" in stdout:
        stdout = stdout[stdout.index(">"):]

        with io.StringIO(stdout) as fasta_output:
            fasta_list = read_fasta(fasta_output)
    else:
        fasta_list = []
    return fasta_list


        

def mafft_wrapper(input_fasta, options="--globalpair", tmp_dir="/tmp"):
    
    num, tmpfasta = tempfile.mkstemp(prefix="mafft", suffix=".fasta", dir=tmp_dir)

    with open(tmpfasta, "w") as ofile:
        for title, seq in input_fasta:
            ofile.write(">%s\n%s\n"%(title,seq))

    MAFFT_COMMAND = "mafft {options} {input}".format(options=options, input=tmpfasta)
    stdout, stderr, returncode = run_output(MAFFT_COMMAND)


    return stdout


def viruligndf(fasta_entries):

    name, ref = fasta_entries[0]

    name = name.split(" ")[0]    
    current_pos = 0
    within_insertion = False


    columns_seq = []
    for aa in ref:
        if aa != '-':
            current_pos += 1
        
        columns_seq.append("%s_%s"%(name, current_pos))

    for f in re.finditer("-+",ref):
        start, end = f.span()
        columns_seq[start:end] = ["%sins%s"%(c,i+1) for i,c in enumerate(columns_seq[start:end])]

    df_seq = pd.DataFrame(list(map(list,[r[1] for r in fasta_entries])),columns=columns_seq)

    df_name = pd.DataFrame(dict(seqid=[r[0] for r in fasta_entries]))

    df = pd.concat([df_name, df_seq], axis=1)


    return df

if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", type=str, dest="input_fasta", help="Input FASTA file")
    parser.add_argument("-n", type=str, dest="output_fasta_nt", help="Output FASTA file (nuc)")
    parser.add_argument("-a", type=str, dest="output_fasta_aa", help="Output FASTA file (aa)")
    parser.add_argument("-c", type=str, dest="output_csv_aa", help="Output CSV  file (aa - virulign positional table output)")
    parser.add_argument("-r", type=str, dest="virulign_ref_file", help="Initial Virulign reference file")
    parser.add_argument("-s", type=int, dest="num_alignments", help="Number of blast alignments", default=20)
    parser.add_argument("-t", type=str, dest="tmp", default="/tmp", help="Location of temporary directory")
    parser.add_argument("-A", dest="aga_only", action="store_true", help="Skip virulign and align everything with AGA")
    
    args = parser.parse_args()

    # Intial FASTA

    fasta_initial = read_fasta(args.input_fasta)
    n_query = len(fasta_initial)

    # Initial alignment
    fasta_initial_alignment_nt = initial_alignment(args.input_fasta, args.virulign_ref_file, maxFrameShifts=0, alphabet="Nucleotides")
    fasta_initial_alignment_aa = initial_alignment(args.input_fasta, args.virulign_ref_file, maxFrameShifts=0, alphabet="AminoAcids")
    reference_nt = list(fasta_initial_alignment_nt[0])
    reference_aa = list(fasta_initial_alignment_aa[0])
   
    fasta_initial_alignment_nt = fasta_initial_alignment_nt[1:] 
    fasta_initial_alignment_aa = fasta_initial_alignment_aa[1:] 

    reference_nt[0] = reference_nt[0].split(' ')[0]
    reference_aa[0] = reference_aa[0].split(' ')[0]

    reference_nt = tuple(reference_nt)
    reference_aa = tuple(reference_aa)


    # Missing alignemnts

    if args.aga_only:
        missing_entries = fasta_initial
        fasta_initial_alignment_nt = []
        fasta_initial_alignment_aa = []

    else:
        missing_entries = find_missing_entries(fasta_initial, fasta_initial_alignment_nt)
    temp_alignments_aa = []
    temp_alignments_nt = []
    n_missing_entries = len(missing_entries)

    

    if n_missing_entries > 0:
        code, missing_entries_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta",dir=args.tmp)
        write_fasta(missing_entries, missing_entries_file)

            
        gb = agagen.create_aga_gb(reference_nt[1].replace("-","").replace(".",""), prot=reference_nt[0])
        aga_alignments_nt = []
        aga_alignments_aa = []

        for title,seq in missing_entries:

            aga = agagen.run_aga((title,seq,), genbank_str=gb, delete_temp=True)
            aa_align, nt_align, full_align = aga

            aa_align = agagen.extract_aligned_seqs(aa_align,title=title,full=False,replace_char="X")
            nt_align = agagen.extract_aligned_seqs(nt_align,title=title,full=False,replace_char="N")

            full_align = list(full_align[1])
            full_align[0] = title

            aga_alignments_nt.extend(nt_align)
            aga_alignments_aa.extend(aa_align)
            temp_alignments_aa.extend(aa_align)
            temp_alignments_nt.extend(nt_align)

        
        missing_entries_aga = find_missing_entries(missing_entries, aga_alignments_nt)


        missing_entries = missing_entries_aga

    final_alignments_aa = fasta_initial_alignment_aa + temp_alignments_aa
    final_alignments_nt = fasta_initial_alignment_nt + temp_alignments_nt

    if len(final_alignments_aa) > 0:
        with open(args.output_fasta_aa, "w") as ofile_aa:
            s = mafft_wrapper([reference_aa]+final_alignments_aa, options="--localpair --maxiterate 1000 --retree 10")
            ofile_aa.write(s)

        with open(args.output_fasta_nt, "w") as ofile_nt:
            s = mafft_wrapper([reference_nt]+final_alignments_nt)
            ofile_nt.write(s)

    
    try:
        virulign_aa = read_fasta(args.output_fasta_aa)
        virulign_df_aa = viruligndf(virulign_aa)
        virulign_df_aa.to_csv(args.output_csv_aa, index=False)
    except FileNotFoundError:
        pathlib.Path(args.output_csv_aa).touch()
        pathlib.Path(args.output_fasta_aa).touch()
        pathlib.Path(args.output_fasta_nt).touch()


    


            


        

                




            



















    


