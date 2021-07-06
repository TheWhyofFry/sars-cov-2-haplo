import os, sys

import pysam

import argparse
import coordremap



def getReferenceName(bam, ref_num):
    
    bam = pysam.AlignmentFile(bamfile, "rb")

    # For now, only the first reference

    reference_name = bam.references[ref_num]

    bam.close()

    return reference_name


def getRefCountsLengths(bam):

    return dict([(ref, (bam.count(ref), bam.lengths[i],),) for i, ref in enumerate(bam.references) ])


def str_search(s, r):

    return True in [s.lower() in _r.lower() for _r in r]

def nameSearch(bam, names):

    return [ref for ref in bam.references if str_search(ref, names) ]

def bestMatch(ref_counts_lengths, match_len, names=[]):
    
    best_ref, best_counts = "", 0
    names = names if names != [] else list(ref_counts_lengths.keys())

    for name in names:
        
        count, length = ref_counts_lengths[name]

        if length == match_len:
            if count > best_counts:
                best_ref, best_counts = name, count


    return best_ref




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i",dest="ref_fasta", type=str, help="Input reference FASTA")
    parser.add_argument("-b", dest="bam_file", type=str, help="Input BAM file")

    parser.add_argument("-o", dest="pseudo_fasta", type=str, help="Output pseudo reference FASTA")
    parser.add_argument("-n", dest="nth_alignment", type=int, default=-1, help="Choose the n-th reference in the file (0-based)")
    parser.add_argument("-s", dest="ref_size", default=0, type=int, help="Length of reference")
    parser.add_argument("-k", dest="known_references", type=str, default="", help="Comma separated list of references. A string matching is performed, so the trailing bits can be excluded")
    args = parser.parse_args()

    ref_title, ref_fasta = coordremap.readfasta(args.ref_fasta)[0]
    
    ref_title_search = ref_title.split(".")[0]

    ref_size = args.ref_size if args.ref_size > 0 else len(ref_fasta.strip())

    bam = pysam.AlignmentFile(args.bam_file, "rb")
    
    if (args.nth_alignment == -1) or (args.nth_alignment > (len(bam.references) - 1)):


        names = args.known_references.split(",") + [ref_title]

        name_matches = nameSearch(bam, names)

        read_counts_lengths = getRefCountsLengths(bam)

        pseudo_ref_title =  bestMatch(read_counts_lengths, ref_size, name_matches)
        
    else:
        pseudo_ref_title = bam.references(args.nth_alignment)



    if pseudo_ref_title == "":
        print("Error: no references in the BAM are appropriate")
        sys.exit(1)



    

    print("Pseudo ref:",pseudo_ref_title)
    with open(args.pseudo_fasta,"w") as ofile:
        fasta = ">%s\n%s\n"%(pseudo_ref_title, ref_fasta)
        ofile.write(fasta)






