import pandas as pd
import sys

blockdf = pd.read_csv(sys.argv[1])

chrLenDict = {}
chrLenList = [248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468]
dfList = []
#for i in range(1,2): #for testing
for i in range(1,23):
    dfList.append(blockdf[blockdf["chrom"] == "chr"+str(i)])
    chrLenDict["chr"+str(i)] = chrLenList[i-1]
dfList.append(blockdf[blockdf["chrom"] == "chrX"])
chrLenDict["chrX"] = 156040895
#dfList.append(blockdf[blockdf["chrom"] == "chrY"])
#chrLenDict["chrY"] = 57227415

genomeSpan = 0
chunkLenRed = 0
coveredTotal = 0
for i in range(0, len(dfList)):
    tempdf = pd.DataFrame(columns=["left","right"])
    startEndList = []
    dfList[i] = dfList[i].sort_values(by=["left"])
    idxList = list(dfList[i].index.values)
    chrSpan = chrLenDict[dfList[i].at[idxList[0],"chrom"]]
    
    for count, j in enumerate(idxList):
        tempLeft = dfList[i].at[j, "left"]
        tempRight = dfList[i].at[j, "right"]
        startEndList.append(tempLeft)
        startEndList.append(tempRight)
        chunkLenRed = chunkLenRed + abs(dfList[i].at[j, "right"] - dfList[i].at[j, "left"])
        
        if count == 0:
            tempdf.loc[len(tempdf)] = [tempLeft, tempRight]
            currRight = tempRight
        if tempLeft > currRight:
            tempdf.loc[len(tempdf)] = [tempLeft, tempRight]
            currRight = tempRight
        elif tempLeft < currRight:
            if tempRight > currRight:
                tempdf.at[len(tempdf)-1, "right"] = tempRight
                currRight = tempRight
            
    chrSpan = chrLenDict[dfList[i].at[j,"chrom"]]
    genomeSpan = genomeSpan + chrSpan
    
    for k in range(0, len(tempdf)):
        coveredTotal = coveredTotal + (tempdf.at[k,"right"]-tempdf.at[k,"left"])
    
genomeCovRed = float(chunkLenRed)*100/float(genomeSpan)
genomeCovActual = float(coveredTotal)*100/float(genomeSpan)

print("Percent of human genome covered (With repeats): " + str(genomeCovRed) + "%\n")
print("Percent of human genome covered (No repeats): " + str(genomeCovActual) + "%\n")
