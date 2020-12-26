import pandas as pd
import argparse
import re
import numpy as np
def fixheadernumbering(header,padding=4):
    header_template = "%%s_%%0%sdins%%0%sd"%(padding,padding)
    
    if "ins" in header:
        gene, position = header.split("_")

        ins_index = position.index("ins")
        pos = int(position[:ins_index])
        trail_num = int(position[ins_index+3:])
        outheader = header_template%(gene,pos,trail_num)
    else:
        outheader = header

    return outheader


class Virualign:

    def __init__(self, positional_table_file,gene="gene"):
        self.protein_pos = []
        self.gene = gene
        self.virulign_df = self._read_pos_table(positional_table_file)
        self.virulign_df_m = self._melt_pos_table(self.virulign_df)
        self.virulign_df_collapsed = self.virulign_df_m.pivot(values="aa",index="seqid",columns="pos")




    def _read_pos_table(self,pos_table_file):

        df = pd.read_csv(pos_table_file)
        df.columns = [fixheadernumbering(col) for col in df.columns]

        return df

    def _melt_pos_table(self, virulign_df):

        m = virulign_df.melt(id_vars=["seqid"]).sort_values(["seqid","variable"])
        m[["gene","pos"]] = m.variable.str.split("_",expand=True)
        m = m.sort_values(["seqid","gene","variable"])
        m["value"] = m["value"].fillna("X").values
        m[m["value"].str.len() > 1]["value"] = "X"
        m["value"][m.value.str.len() > 1] = "X"
        ins_re = re.compile("ins.+")
        m["pos"] = m.pos.apply(lambda x:re.sub(ins_re,"",x)).astype(np.int)

    

        m = m.groupby(["seqid","gene","pos"]).agg({"value":lambda x:"".join(x)}).reset_index()
        
        m = m.rename(columns={"value":"aa"})
        return m



    def _make_seqs(self,gene, start, end,pad_left=0,pad_right=0,annotation="N/A"):
        
        
        v_df = self.virulign_df_m[self.virulign_df_m.gene == gene].pivot(index="seqid",columns="pos",values="aa")
        gene_end = self.virulign_df_m[self.virulign_df_m.gene == gene].pos.max()
        pad_left_str  = "." * (pad_left-start)
        pad_right_str = "." * ((end+pad_right)-gene_end)
        
        pad_left = pad_left if start > 5 else -1
        pad_right = pad_right if (end+pad_right) < gene_end else 1 

        core = np.apply_along_axis(lambda x:"".join(x), 1, v_df[list(range(start,end+1))].values)

        if pad_left > 0:
            pad_left_str = np.apply_along_axis(lambda x:"".join(x)+pad_left_str, 1, v_df[list(range(start-pad_left,start))].values)
        


        if end < gene_end:
            pad_right_str = np.apply_along_axis(lambda x:"".join(x)+pad_right_str, 1, v_df[list(range(end+1,end+pad_right))].values)
        
        
    
        return pd.DataFrame(dict(seqid=v_df.index.values,annotation=annotation, start=start,end=end, core=core, pad_left=pad_left_str, pad_right=pad_right_str))


        

def repairseq(virulign_df_m, ref, excluded="X"):

    v_ref = virulign_df_m[virulign_df_m.seqid == ref]

    v_ref = v_ref.rename(columns={"aa":"aa_ref"})


    v_rest = virulign_df_m[virulign_df_m.seqid != ref].merge(v_ref[["gene","pos","aa_ref"]],on=["gene","pos"])
    v_rest["new_aa"] = v_rest.aa.values[:]
    v_rest.loc[v_rest.aa == excluded, "new_aa"] = v_rest[v_rest.aa == excluded].aa_ref.values
    

    return v_rest

        
def get_muts(v):

     vfilt = v[v.aa_ref != v.new_aa]

     vf_pos = vfilt.drop_duplicates(["gene","pos"])[["gene","pos"]]
     vf = pd.merge(v, vf_pos, on=["gene","pos"])
     if len(vf) == 0:
         return False
     vf["muts"] = vf.apply(lambda x:"".join("%s%s%s"%(x.aa_ref,x.pos,x.new_aa)),axis=1)
     
     return vf

def splitseqid(v):

    v = pd.concat([v, v.seqid.str.split("/",expand=True).rename(columns={0:"samplename",1:"haplonum",2:"prop"})],axis=1)
    
    v["prop"] = v["prop"].astype(np.float)
    return v
        

def aggregate_muts(v):
    if len(v) == 0:
        return None
    vm = get_muts(v).groupby(["samplename","gene","muts","pos","new_aa"]).aggregate({"prop":sum, "haplonum":lambda x: "/".join(x)}).reset_index()
    vm["pos"] = vm["pos"].astype(np.int)
    return vm



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", dest="inputfile", help="Virulign Position Table")
    parser.add_argument("-n", dest="normalfile", help="Normal table output file")
    parser.add_argument("-m", dest="mutaggfile", help="Aggregated mutation file")
    parser.add_argument("-a", dest="aggfile", help="Aggregated residue file")


    args = parser.parse_args()

    v = Virualign(args.inputfile)

    ref_seq = v.virulign_df.iloc[0,0]

    y = repairseq(splitseqid(v.virulign_df_m), ref_seq).sort_values(["gene","pos","prop"],ascending=[True,True,False])

    agg = y.groupby(["samplename","gene","pos","new_aa"]).aggregate({"prop":sum, "haplonum":lambda x: "/".join(x)}).reset_index().sort_values(["samplename","gene","pos","prop"],ascending=[True,True,True,False])

    agg_muts = aggregate_muts(y[y.aa_ref != y.new_aa])
    if agg_muts is not None:
        agg_muts = agg_muts.sort_values(["gene","pos","prop"],ascending=[True,True,False])
        agg_muts.to_csv(args.mutaggfile)
    else:
        with open(args.mutaggfile,"w") as ofile:
            ofile.write("")
    agg.to_csv(args.aggfile)

    y.iloc[:,2:].to_csv(args.normalfile)


