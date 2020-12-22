import numpy as np
import pandas as pd



class CoordRemapper(object):

    def __init__(self, f,reverse=False, seq_idx=[0,1]):
        
        fastaentries = readfasta(f)
        if reverse:
            ref, target = fastaentries[seq_idx[0]], fastaentries[seq_idx[1]]
        else:
            target, ref = fastaentries[seq_idx[0]], fastaentries[seq_idx[1]]
        self.ref_seq = ref[1]
        self.target_seq = target[1]

        self.ref_remap = selfremap(self.ref_seq)
        self.target_remap = selfremap(self.target_seq)
        self.ref_remap = selfremap(self.ref_seq)
        self.target_remap_rev = pd.Series(self.target_remap.index, index=self.target_remap.values)
        self.ref_remap_rev = pd.Series(self.ref_remap.index, index=self.ref_remap.values)

    def remaptarget(self, start, stop):
        s = slice(start, stop)

        return self.ref_remap_rev[self.target_remap.loc[s].values]
    def gettargetpos(self, start, stop):

        pos_ref = self.ref_remap[start:stop].values

        pos_target_idx  = sorted(set(pos_ref).intersection(self.target_remap_rev.index))

    


        return (self.target_remap_rev[pos_target_idx].values, self.ref_remap_rev[pos_target_idx],)





def coordremapper(s1,s2):

    current_pos1 = 0
    current_pos2 = 0
    pos_list    = []

    #s1 is the aligned seq
    
    while current_pos2 < len(s2):
        if s1[current_pos1] == "-":
            try:
                if s2[current_pos2] == "-":
                    pos_list.append((current_pos1, current_pos2))
                    current_pos1 += 1
                    current_pos2 += 1
                    continue
            except:
                print(current_pos1,current_pos2)
                print(len(s1),len(s2))
                print(pos_list[-30:])
                raise
            current_pos1 += 1
        else:
            pos_list.append((current_pos1, current_pos2))
            current_pos1 += 1
            current_pos2 += 1
    pos_index = [p2+1 for p1,p2 in pos_list]
    pos_ref = [p1+1 for p1,p2 in pos_list]

    return pd.Series(pos_ref, index=pos_index, dtype=np.int32)

#Remap assuming no gaps was in original
def selfremap(s,illegals=".-"):
    t = s
    for i in illegals:
        t = t.replace(i,"")
    
    s = s.replace(".","-")
     
    return coordremapper(s,t)

def readfasta(f,upper=True,replace_dot=True):
    #fasta_dict = {}
    fasta_list = []
    curseq = ""
    
    with open(f, "r") as fastafile:
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
    #fasta_dict[currentkey] = curseq
    if replace_dot:
        curseq = curseq.upper().replace(".","-") if upper else curseq
    else:
        curseq = curseq.upper() if upper else curseq
    curseq = curseq.upper()
    try:
        fasta_list.append((currentkey, curseq))
    except:
        print(f)
        #print currentkey
        print("Curseq:", curseq)
        raise
    return fasta_list


def seqseries(s):
    return pd.Series(list(s), index=list(range(1,len(s)+1)))


def revmapper(target, s):

    target_rev = pd.Series(target.index, index=target.values)

    s_rev = pd.Series(s.index, index=s.values)

    s_target_rev = s_rev[target_rev.index]
    return target_rev, s_rev, s_target_rev


