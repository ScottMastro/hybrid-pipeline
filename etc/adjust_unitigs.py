import pandas as pd
import re

def utg_name(tigId):
    return "utg" + re.sub("\D", "", str(tigId)).zfill(8)
def tig_name(tigId, suffix=""):
    return "tig" + re.sub("\D", "", str(tigId)).zfill(8) + suffix

def map_unitig_contig():

    # input 
    print("Provide path of *.contigs.layout.readToTig output file from Canu:")
    contigLayout = input()
    
    cMap = pd.read_csv(contigLayout, sep='\t')

    print("Provide path of *.unitig.layout.readToTig output file from Canu:")
    unitigLayout = input()   
    
    uMap = pd.read_csv(unitigLayout, sep='\t')

    readMap = cMap.merge(uMap, on="#readID", how='right', suffixes=("_ctg", "_utg"))    
    readMap["minpos_ctg"] = readMap[["bgn_ctg","end_ctg"]].min(axis=1)
    readMap["maxpos_ctg"] = readMap[["bgn_ctg","end_ctg"]].max(axis=1)

    tigMapMin = readMap.groupby(["tigID_utg", "tigID_ctg"])["minpos_ctg"].min()
    tigMap = tigMapMin.index.to_frame()
    tigMap = tigMap.reset_index(drop=True)
    tigMap.columns = ["name", "chrom"]
    tigMap["start"] = tigMapMin.reset_index(drop=True).astype("int32")
    tigMapMax = readMap.groupby(["tigID_utg", "tigID_ctg"])["maxpos_ctg"].max()
    tigMapMax.columns = ["name", "chrom", "end"]
    tigMap["end"] = tigMapMax.reset_index(drop=True).astype("int32")
    tigMap = tigMap[["chrom", "start", "end", "name"]]
    
    print("[optional] Provide contig suffix (ex. _pilon_pilon) :")

    suffix = input()   
    
    tigMap["chrom"] = tigMap["chrom"].apply(tig_name, args=(suffix,))
    tigMap["name"] = tigMap["name"].apply(utg_name)
    
    print("Provide output path:")
    outpath = input()   
    if outpath == "":
        outpath = "out.bed"

    tigMap.to_csv(outpath, sep='\t', header=False, index=False)
    
    return tigMap

if __name__== "__main__":
    map_unitig_contig()
    exit()