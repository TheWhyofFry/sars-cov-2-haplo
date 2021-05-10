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
    print("MISSING ENTRIES:", missing_entries, len(query_fasta), len(virulign_fasta))
    return [(title, seq) for title,seq in query_fasta if title in missing_entries]


def blast(missing_entries, blast_db, blast_prog="blastn", other_options=""):


    blast_command = "{blast_prog} -query {query} -db {blast_db} -outfmt '6 qseqid sseqid pident length mismatch qstart qend sstart send bitscore qseq sseq' {other_options}"



    stdout, stderr, returncode = run_output(blast_command)



def initial_alignment(fasta_file, virulign_ref_file, virulign_command="virulign",
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
        #stdout = stdout[stdout.index(">"):]

        with io.StringIO(stdout) as fasta_output:
            fasta_list = read_fasta(fasta_output)
    else:
        fasta_list = []
    return fasta_list


        
def do_secondary(missing_entries, ref_entry, virulign_params={}):
    
    code, temp_ref_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta", dir="/tmp")
    code, temp_missing_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta", dir="/tmp")
    print("TEMP MISSING FILE:", temp_missing_file)
    write_fasta(missing_entries, temp_missing_file)

    with open(temp_ref_file, 'w') as temp_ref:
        temp_ref.write(">%s\n%s\n"%(ref_entry))

    
    # We'll just directly assume frameshifts are OK
    secondary_alignment = initial_alignment(temp_missing_file, temp_ref_file,  **virulign_params)

    secondary_missing_entries = find_missing_entries(missing_entries, secondary_alignment)

    if len(secondary_missing_entries) > 0:
        missing_entries = secondary_missing_entries
    

    #os.unlink(temp_ref_file)
    #os.unlink(temp_missing_file)


    return secondary_alignment, secondary_missing_entries

def hamming(s1,s2):
    return sum([1 for a,b in zip(s1,s2) if a != b])

def get_valid_hits(blastx_hits):

        valid_hits = list(blastx_hits.sstart.diff().values[1:] >= 0) + [True]
       

        blastx_hits = blastx_hits[valid_hits] 
        print(blastx_hits[["qstart","qend", "sstart","send"]])
        print(valid_hits)
        valid_hits = [True] + list(blastx_hits.send.diff().values[1:] >= 0) 
        

        blastx_hits = blastx_hits[valid_hits]

        print(blastx_hits[["qstart","qend","sstart","send"]])
        return blastx_hits

def remove_embedded(blastx_hits):
    
    if len(blastx_hits) == 1:
        return blastx_hits


    blastx_hits = blastx_hits.sort_values(["sstart","send"], ascending=[True,False]).copy()
    # Debug counter
    counter = 0 
    l = len(blastx_hits)
    while True:
        if counter > l:
            break
        sstart = blastx_hits.sstart.values
        send   = blastx_hits.send.values


        within_start = sstart[1:] > sstart[:-1]
        within_end   = send[1:] < send[:-1]

        logical = [True]
        logical.extend(~(within_end & within_start))
        
        if len(blastx_hits[logical]) == len(blastx_hits):
            return blastx_hits
        blastx_hits = blastx_hits[logical]
        
        counter += 1



    return blastx_hits


def recreate_seq(blastx_hits):

    # Sort according to query start position
    blastx_hits = blastx_hits.sort_values("qstart")

    # Remove spurious alignments 
    # Check if the qstart/sstart sequences follow ascendingly


    # If multiple hits, we need to check how to concatenate them
    # If there are overlaps, resolve which alignment is the best (hamming distance - with ties, the first is used)
    # If there are not, replace the frameshifted region (non-matchin region) with X
    chain = []
    if len(blastx_hits) > 1:
        blastx_hits = get_valid_hits(blastx_hits)
        # iterate through 
        current_hit = blastx_hits.iloc[0,]
        current_hit.mask = 0
        current_hit.qseq_remap = remap(current_hit.qseq)
        current_hit.sseq_remap = remap(current_hit.sseq)

    
        blastx_hits = blastx_hits.iloc[1:,]


        for idx,hit_ in blastx_hits.iterrows():
            hit = hit_.copy()
            hit.mask = 0
            hit.sseq_remap = remap(hit.sseq)
            hit.qseq_remap = remap(hit.qseq)
            print("Current hit end: {qend}; hit start: {qstart}".format(qend=current_hit.qend, qstart=hit.qstart))
            if hit.qstart < current_hit.qend:
                aa_pos_span = math.ceil((current_hit.qend - hit.qstart + 1)/3)
                
                aa_pos_remap_qseq_current_hit = current_hit.qseq_remap.index.max() - aa_pos_span
                aa_pos_remap_qseq_hit = hit.qseq_remap[aa_pos_span]

                hamming_hit = hamming(hit.qseq[:aa_pos_remap_qseq_hit], hit.sseq[:aa_pos_remap_qseq_hit])
                hamming_current_hit = hamming(current_hit.qseq[:aa_pos_remap_qseq_current_hit], current_hit.sseq[:aa_pos_remap_qseq_current_hit])
                

                if hamming_hit < hamming_current_hit:
                    o = current_hit.qseq_remap.index.max() - aa_pos_span + 1
                    trim_pos = current_hit.qseq_remap.loc[current_hit.qseq_remap.index.max() - aa_pos_span + 1]
                    qseq = current_hit.qseq
                    print("Trim pos:",trim_pos, aa_pos_span, o, len(current_hit.qseq), hit.sstart, hit.send)
                    current_hit.qseq = current_hit.qseq[:trim_pos]
                    print(qseq, len(qseq), current_hit.qseq, len(current_hit.qseq))
                    

                else:
                    print("Hit adjusted")
                    mask = hit.qseq_remap[aa_pos_span-1]
                    hit.qseq = hit.qseq[mask:]
                    hit.sseq = hit.sseq[mask:]
                    hit.qseq_remap = remap(hit.qseq)
                    hit.sseq_remap = remap(hit.sseq)



            else:
                delta = hit.qstart - current_hit.qend + 1
                print("Delta: ",delta)
                aa_pos_span = math.ceil(delta / 3)
                current_hit.qseq += "X" * aa_pos_span

            chain.append(current_hit.qseq)

            current_hit = hit
        chain.append(current_hit.qseq[current_hit.mask:])

    else:
        chain.append(blastx_hits.qseq.value[0])


    return (blastx_hits.qseqid.values[0], "".join(chain))


def blastx_align(fasta_list, blastdb):
    

    temp_blast_query = mkstemp(prefix="virulign-tools-blastx",suffix=".fasta")

    write_fasta(fasta_list, temp_blast_query)
    
    BLASTX_COMMAND = "blastx -db {db} -query {query} -outfmt '6#qseqid#sseqid#bitscore#qstart#qend#qseq'"

    blast_output, stderr, returncode = run_command(BLASTX_COMMAND.format(blastdb, temp_blast_query))

    if returncode == 0:
        with open(io.StringIO(blast_output)) as blast_output_f:
            blastx_df = pd.read_table(blast_output_f, names=["qseqid", "sseqid", "bitscore", "qstart", "qend", "qseq"]).\
                            sort_values(["qseqid","qstart"])

            blastx_df_maxbitscores = blastx_df.group_by(["qseqid","sseqid"]).\
                    aggregate({'bitscore': lambda x:x.bitscore.sum()}).\
                    sort_values("bitscore",ascending=False).reset_index().\
                    groupby("qseqid").first().\
                    reset_index()



            
            



    return

def write_fasta(fasta_list, fasta_file):

    with open(fasta_file, "w") as ofile:
        for title,seq in fasta_list:
            ofile.write(">%s\n%s\n"%(title,seq))


    return


def mafft_wrapper(input_fasta, options="--localpair", tmp_dir="/tmp"):
    
    num, tmpfasta = tempfile.mkstemp(prefix="mafft", suffix=".fasta", dir=tmp_dir)

    with open(tmpfasta, "w") as ofile:
        for title, seq in input_fasta:
            ofile.write(">%s\n%s\n"%(title,seq))

    MAFFT_COMMAND = "mafft {options} {input}".format(options=options, input=tmpfasta)
    stdout, stderr, returncode = run_output(MAFFT_COMMAND)


    return stdout

def reorder_fasta(fasta_list, reorder_dict):

    alignments_ = [True] * len(fasta_list)

    for i, (title,seq) in enumerate(fasta_list):
        alignments_[reorder_dict[title]] = fasta_list[i]

    return alignments_


if __name__ == "__main__":
    import argparse
    import textwrap

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", type=str, dest="input_fasta", help="Input FASTA file")
    parser.add_argument("-n", type=str, dest="output_fasta_nt", help="Output FASTA file (nuc)")
    parser.add_argument("-a", type=str, dest="output_fasta_aa", help="Output FASTA file (aa)")
    parser.add_argument("-r", type=str, dest="virulign_ref_file", help="Initial Virulign reference file")
    #parser.add_argument("-b", type=str, dest="blastn", help="BLASTN DB", default="")
    #parser.add_argument("-B", type=str, dest="blastp", help="BLASTP DB", default="")
    #parser.add_argument("-X", action='store_true', dest="do_blastx", help="If all else fails, align with BLASTX")
    #parser.add_argument("-R", action='store_true', dest="do_selfref", help="Align failed sequences to sequences that passed first")
    #parser.add_argument("-N", action='store_true', dest="do_blastn", help="Repeat alignment for failed sequences using the specified BLASTN DB")
    parser.add_argument("-s", type=int, dest="num_alignments", help="Number of blast alignments", default=20)

    
    args = parser.parse_args()

    # Intial FASTA

    fasta_initial = read_fasta(args.input_fasta)

    fasta_order = dict([(title, i+1) for i, (title,seq) in enumerate(fasta_initial)])


    

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

    missing_entries = find_missing_entries(fasta_initial, fasta_initial_alignment_nt)
    temp_alignments_aa = []
    temp_alignments_nt = []
    n_missing_entries = len(missing_entries)
    
    print("Missing entries:", n_missing_entries)

    if n_missing_entries > 0:
        code, missing_entries_file = tempfile.mkstemp(prefix="virulign-tools", suffix=".fasta",dir="/tmp")
        print("MISSING ENTRY FILE:", missing_entries_file, code)
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
            print("Current recovered:")
            print([title for seq, title in aa_align])
            temp_alignments_aa.extend(aa_align)
            temp_alignments_nt.extend(nt_align)

        
        missing_entries_aga = find_missing_entries(missing_entries, aga_alignments_nt)


        missing_entries = missing_entries_aga

    print("AGA RECOVERED:")
    print([title for title,seq in temp_alignments_aa])
    final_alignments_aa = fasta_initial_alignment_aa + temp_alignments_aa
    final_alignments_nt = reorder_fasta(fasta_initial_alignment_nt + temp_alignments_nt, fasta_order)
    

    if len(final_alignments_aa) > 0:
        with open(args.output_fasta_aa, "w") as ofile_aa:
            
            s = mafft_wrapper([reference_aa]+fasta_alignments_aa)
            ofile_aa.write(s)

        with open(args.output_fasta_nt, "w") as ofile_nt:
            s = mafft_wrapper([reference_nt]+final_alignments_nt)
            ofile_nt.write(s)

        print("Wrote %s AA sequences and %s NT sequences:"%(len(fasta_initial_alignment_nt), len(fasta_initial_alignment_aa)))

    if len(missing_entries) > 0:
        print("Some sequences could not be aligned: "%("\n".join([title for title,seq in missing_entries])))

    


"""
        blastn_fasta_entries = []

        if (args.do_blastn):


            # Get N hits from a blast db to try and use these seqs as reference
            # Note that these are _valid_ references, i.e. seqs where virulign could resolve
            # A fundamental reference sequence (e.g. HXB2 in HIV-1) with a query
            # and WITHOUT any frame shift compensation

            BLASTN = "blastn -db {db} -query {query} -num_alignments {n} -outfmt 6#qseqid#sseqid#bitscore".format(db=args.blastn,
                                                                                                                    query=missing_entries_file,
                                                                                                                    n=args.num_alignments)
            print(BLASTN)
            blastn_output, stderr, returncode = run_output(BLASTN)
            with io.StringIO(blastn_output) as blastn_output_f:
                blastn_df = pd.read_table(blastn_output_f,names=["qseqid", "sseqid", "bitscore"])
                blast_hits = blastn_df.value_counts("sseqid").sort_values(ascending=False).reset_index().sseqid.values
                print("BLAST HITS:")
                blastn_entries, stderr, returncode = run_output('blastdbcmd -entry "{entry}" -db {db}'.format(entry=','.join(blast_hits),db=args.blastn))
                
                with io.StringIO(blastn_entries) as blast_file:
                    blastn_fasta_entries = read_fasta(blast_file)
                    print("BLASTN Entries:", blastn_fasta_entries)

            



        if (args.do_selfref) and (len(fasta_initial_alignment) > 1):
            fasta_initial_alignment_ref = fasta_initial_alignment[:]
        else:
            fasta_initial_alignment_ref = []


            








            


        
        # Placeholder
        rep = False
        if rep:
            print("Finding appropriate references ... ")
            alternative_alignments = fasta_initial_alignment_ref + blastn_fasta_entries
            print([title for title,seq in alternative_alignments])
            for fasta_entry in alternative_alignments:
                title, seq = fasta_entry
                print(title)
                secondary_alignment, missing_entries = do_secondary(missing_entries, fasta_entry, virulign_params=dict(alphabet="AminoAcids", exportwithref="no"))

                if len(secondary_alignment) == 0:
                    continue
                
                temp_alignments.extend(secondary_alignment)

                n_missing_entries -= len(secondary_alignment)

                if n_missing_entries == 0:
                    break
        if (n_missing_entries > 0) and (args.do_blasx):

            blastx = blastx_align(missing_entries, args.blastp)

            if len(blastx) > 1:
                pass                    
            pass
"""
            


        

                




            



















    


