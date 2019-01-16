import pandas as pd
import numpy as np 

def constructMatrix(blockdf, threshold=35):
    
    blockdf = blockdf.sort_values(by=["contig", "chrom", "left"])
    blockdf["% cov"] = np.nan
    blockdf["below_nova"] = blockdf["contig"].shift(-1)
    blockdf["below_canu"] = blockdf["chrom"].shift(-1)

    nova = blockdf.drop_duplicates("contig")
    nova = [ (x[1]["contig"], x[1]["span"]) for x in nova.iterrows()]
    canu = blockdf.drop_duplicates("chrom")
    canu = [ (x[1]["chrom"], x[1]["total_chunk"]*1000) for x in canu.iterrows()]

    novaTigs =  [n[0] for n in nova] 
    canuTigs = [c[0] for c in canu]

    matrix = pd.DataFrame(columns=novaTigs, index=canuTigs)

    count = 1
    indices = []
    edgelist = []

    idxList = list(blockdf.index.values)

    print("Iterating through blocks")

    for index in idxList:
        currentCanu = blockdf.at[index,"chrom"]
        currentNova = blockdf.at[index,"contig"]
        blockcov = round(float((blockdf.at[index,"right"]-blockdf.at[index,"left"])*100)/float(blockdf.at[index, "span"]), 3)
        indices.append(index) 
        
        if blockdf.at[index,"below_nova"] == blockdf.at[index,"contig"] and blockdf.at[index,"below_canu"] == blockdf.at[index,"chrom"]:
            count = count + 1
                    
        elif count == 1: #One nova tig in canu tig
            #if currentNova in novasmol:
            blockdf.at[index, "% cov"] = blockcov
            if blockcov < threshold:
                indices = []
                continue
            
            matrix.at[currentCanu, currentNova] = str(blockcov)+"%", index
            edgelist.append((currentCanu, currentNova))
            
            indices = []
    
        else: #Multiple nova tigs in same canu tig
            tempdf = pd.DataFrame(columns=["left","right"])
            cov = 0.000000000
            for num, i in enumerate(indices):
                tempLeft = blockdf.at[i, "left"]
                tempRight = blockdf.at[i, "right"]
                tempTig = blockdf.at[i, "contig"]
                        
                if num == 0:
                    tempdf.loc[len(tempdf)] = [tempLeft, tempRight]
                    currRight = tempRight
                elif tempLeft > currRight:
                    tempdf.loc[len(tempdf)] = [tempLeft, tempRight]
                    currRight = tempRight
                elif tempLeft < currRight:
                    if tempRight > currRight:
                        tempdf.at[len(tempdf)-1, "right"] = tempRight
                        currRight = tempRight
            for k in range(0, len(tempdf)):
                cov = cov + (tempdf.at[k,"right"]-tempdf.at[k,"left"])
            covpercent = round(cov*100/blockdf.at[index,"span"], 3)
            blockdf.at[index, "% cov"] = covpercent
            
            if covpercent < threshold:
                indices = []
                count = 1
                continue
            
            matrix.at[currentCanu, currentNova] = str(covpercent)+"%", indices
            edgelist.append((currentCanu, currentNova))
    
            indices = []
            count = 1
    
    return matrix

