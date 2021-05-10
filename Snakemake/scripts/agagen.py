import os
import argparse
from runcommand import run_output

import tempfile 


"""

                     /codon_start="1"
                     /transl_table="1"
                     /protein_id="env"
"""

def create_aga_gb(seq, prot):

    assert len(seq) % 3 == 0

    TEMPLATE =  "FEATURES            Location/Qualifiers\n"\
                " CDS                {RANGE}\n"\
                "                    /codon_start=\"1\"\n"\
                "                    /transl_table=\"1\"\n"\
                "                    /protein_id=\"{PROT}\"\n"\
                "ORIGIN\n"\
                " 1 {SEQ}"


    RANGE = "1..%s"%(len(seq))
    

    return TEMPLATE.format(RANGE=RANGE, SEQ=seq, PROT=prot)


# input_seq = (title, seq) or fasta seq
def run_aga(input_seq, genbank_str=None, genbank_file=None, tmp_dir="/tmp",file_prefix="agawrapper", delete_temp=True):
    
    ALIGNMENT = tempfile.mkstemp(prefix="alignment", suffix=".fasta",dir=tmp_dir)[1]
    NT        = tempfile.mkstemp(prefix="nt", suffix=".fasta", dir=tmp_dir)[1]
    AA        = tempfile.mkstemp(prefix="aa", suffix=".fasta", dir=tmp_dir)[1]
    INPUT     = tempfile.mkstemp(prefix="input", suffix=".fasta", dir=tmp_dir)[1]
    GENBANK   = tempfile.mkstemp(prefix="genbank",suffix=".gb", dir=tmp_dir)[1]
    if type(input_seq) is tuple:
        input_seq = ">%s\n%s\n"%input_seq

    
    if genbank_str is not None:
        with open(GENBANK, "w") as ofile:
            ofile.write(genbank_str)
    if genbank_file is not None:
        GENBANK = genbank_file


    with open(INPUT, "w") as ofile:
        ofile.write(input_seq)


    AGA       = "aga --local --protein-aa-alignments={AA} --protein-nt-alignments={NT} {genbank_file} {INPUT} {ALIGNMENT}"

    aga_command = AGA.format(AA=AA, NT=NT, INPUT=INPUT, ALIGNMENT=ALIGNMENT, genbank_file=GENBANK)

    stdout, stderr, returncode = run_output(aga_command)


    if returncode == 0:
        aa_fasta = read_fasta(AA)
        nt_fasta = read_fasta(NT)
        alignment= read_fasta(ALIGNMENT)
    else:
        aa_fasta, nt_fasta, alignment = False, False, False

        
    if delete_temp:
        os.unlink(ALIGNMENT)
        os.unlink(NT)
        os.unlink(AA)
        os.unlink(INPUT)

        if genbank_file is None:
            pass
            os.unlink(GENBANK)

    


    return aa_fasta, nt_fasta, alignment




def extract_aligned_seqs(aga_alignment, title="", full=False, replace_char="N"):

    ref_alignments = aga_alignment[::2]
    target_alignments = aga_alignment[1::2]

    
    if not full:
        target_alignments = [(title,seq.replace(".",replace_char)) for title_, seq in target_alignments]
    else:
    
        target_alignments = [(reftitle,title,seq.replace(".",replace_char),) for (title_, seq), (reftitle, refseq) in zip(target_alignments,ref_alignments)]

    return target_alignments



    







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
            curseq += line.strip()

    fastafile.close()


    fasta_list.append((currentkey, curseq))
    return fasta_list

