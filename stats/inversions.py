import sys
import pandas as pd
import regions as rgn
import numpy as np
import matplotlib.pyplot as plt
import math

def parse_paf(pafFile):
    """Loads file containing minimap2 alignments."""

    colnames = ["qid", "qlen", "qstart", "qend",
                "strand", 
                "rid", "rlen", "rstart", "rend", 
                "nmatch", "alnlen", "qual"]
    
    colorder = [c for c in colnames]
    
    typeMap = {"qlen": int, "qstart": int, "qend": int,
                    "rlen": int, "rstart": int, "rend": int,
                    "nmatch": int, "alnlen": int, "qual": int}
    
    lines = []
    with open(pafFile) as paf:
       line = paf.readline()
       while line:
           cols = line.split("\t")
           l = {c : x for c,x in zip(colnames, cols)}
           
           for i in range(len(colnames), len(cols)):
               col = cols[i]
               l[col[0:2]] = col[5:].rstrip("\n\r")
               if col[0:2] not in colorder: colorder.append(col[0:2])
               
           lines.append(l)
           line = paf.readline()
    
    
    df = pd.DataFrame(lines)
    df = df[colorder]
    
    if "dv" in df.columns:
        typeMap["dv"] = float    
    if "NM" in df.columns:
        typeMap["NM"] = int    


    df = df.astype(typeMap)

    return df


def main():
    inversionFile = "/media/scott/Rotom/hybrid2/CF062_19/out/inversions.txt"
    # "/media/scott/HDD/sickkids/CF062_19/inversions.txt"
     #sys.argv[1]
    pafFile = "/media/scott/Rotom/hybrid2/CF062_19/supernova_to_hg38.paf"
    #"/media/scott/HDD/sickkids/CF062_19/supernova_to_hg38.paf" #sys.argv[2]

    logs = pd.read_csv(inversionFile, sep='\t', header=None)
    alignments = parse_paf(pafFile)
    
    logDict, alignmentDict = dict(), dict()
    for i,log in logs.iterrows():
        region, side, info = log
        name = region.split(":")[0]
        if name not in logDict: logDict[name] = []
        logDict[name].append(log)
        
    for i,alignment in alignments.iterrows():
        name = alignment["qid"]
        if name not in alignmentDict: alignmentDict[name] = []
        alignmentDict[name].append(alignment)

    
    BUFFER=10000
    
    strandCounter = dict()
    x = []
    y = []
    for i,log in logs.iterrows():
        
        region, side, info = log
        chrom, pos, strand = region.split(":")
        startPos, endPos = [int(x) for x in pos.split("-")]
        
        if len(alignmentDict[chrom][0]) < 1:
            print("No alignment to human reference...")
            print(log)
            continue
        
        tigLen = alignmentDict[chrom][0].qlen
        
        regionStart = max(0, startPos-BUFFER)
        regionEnd = min(endPos+BUFFER, tigLen)

        flippedRegion = [rgn.SimpleRegion(chrom, startPos, endPos)]
        flankRegions = [rgn.SimpleRegion(chrom, regionStart, startPos),
                        rgn.SimpleRegion(chrom, endPos, regionEnd)]

        flippedMatch, flankMatch = [], []
        chrs = []
        for alignment in alignmentDict[chrom]:
            refRegion = [rgn.SimpleRegion(chrom, alignment.qstart, alignment.qend, 
                                          annotation=alignment.strand)]
            
            chrs.append(alignment.rid)
            flippedMatch.extend(rgn.regions_intersection(refRegion, flippedRegion))
            flippedLeftover = rgn.regions_difference(flippedRegion, refRegion)

            flankMatch.extend(rgn.regions_intersection(refRegion, flankRegions))
            flankLeftover = rgn.regions_difference(flankRegions, refRegion)

            flippedRegion,flankRegions = flippedLeftover,flankLeftover
    
        flip = [len(r) if r.annotation == '+' else 0 for r in flankMatch] > \
            [len(r) if r.annotation == '-' else 0 for r in flankMatch]
        
        flank1 = rgn.SimpleRegion(chrom, regionStart, startPos)
        unused1 = rgn.region_intersection(flank1, flankRegions)
        flank2 = rgn.SimpleRegion(chrom, endPos, regionEnd)
        unused2 = rgn.region_intersection(flank2, flankRegions)

        if sum([len(r) for r in unused1]) > 0.2 * BUFFER or \
             sum([len(r) for r in unused2]) > 0.2 * BUFFER:
                 print("skipping")
                 continue
        
        num = sum([len(r) if r.annotation == '+' else 0 for r in flankMatch])
        dem = sum([len(r) for r in flankMatch])
        if flip: num = dem-num
        x = (num+1)/(dem+1)

        num = sum([len(r) if r.annotation == '+' else 0 for r in flippedMatch])
        dem = sum([len(r) for r in flippedMatch])
        if flip: num = dem-num
        y = (num+1)/(dem+1)
        
        print(round(x,3),round(y,3))
        print(region)
        FLIPSIZE=BUFFER
        
        def strand_to_int(strand, flip=False):
            if strand == "+" and not flip: return 1
            if strand == "-" and flip: return 1
            return 0
        
        for region in flankMatch: 
            for i in range(region.start, region.end+1):
                if i <= startPos:
                    key = i - startPos - 1
                elif i >= endPos:
                    key = i + FLIPSIZE - endPos + 1
                if key not in strandCounter: strandCounter[key] = []
                strandCounter[key].append(strand_to_int(region.annotation, flip))
    
        strandCounterTemp = dict()
        for region in flippedMatch:
            size = endPos-startPos
            for i in range(region.start, region.end+1):
                key = i - startPos
                key = round((key/size) * (FLIPSIZE-1))
                
                if key not in strandCounterTemp: strandCounterTemp[key] = []
                strandCounterTemp[key].append(strand_to_int(region.annotation, flip))

        for key in strandCounterTemp: 
            value = np.mean(strandCounterTemp[key])
            if key not in strandCounter: strandCounter[key] = []
            strandCounter[key].append(value)

    x = []
    y = []

    for key in strandCounter:
        x.append(key)
        y.append(np.mean(strandCounter[key]))
        
    plt.scatter(x, y, alpha=0.5)
    plt.axvline(x=0)
    plt.axvline(x=FLIPSIZE)
    plt.show()        
  

if __name__ == "__main__":
    main()   

