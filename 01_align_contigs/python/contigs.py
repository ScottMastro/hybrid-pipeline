import pandas as pd
import sys

blockdf = pd.read_csv(sys.argv[1])
blockdf = blockdf.sort_values(by=["contig", "left"])
idxList = list(blockdf.index.values)
tigdf = pd.DataFrame(columns=["contig","tig_size","nblocks","avgblocksize","missing_chunks","n_missing_chunks"])

contigList = []
tigSizeList = []
nblocksList = []
avgBlockSizeList = []
allChunksList = []
missingChunksList = []
fullMissingChunksList = []
nFullMissingChunksList = []
nblocks = 0
blockSize = 0
for count, index in enumerate(idxList):
    nblocks = nblocks + 1
    blockSize = blockSize + blockdf.at[index,"size"]
    for j in range(blockdf.at[index,"chunk_start"],blockdf.at[index,"chunk_end"]+1):
        if j not in allChunksList:
            allChunksList.append(j)
        
    if index == idxList[-1]:
        contigList.append(blockdf.at[index,"contig"])
        tigSizeList.append(blockdf.at[index,"nchunk"])
        nblocksList.append(nblocks)
        avgBlockSizeList.append(blockSize/nblocks)
        for k in range(0,len(allChunksList)):
            if k not in allChunksList:
                missingChunksList.append(k)
        fullMissingChunksList.append(missingChunksList)
        nFullMissingChunksList.append(len(missingChunksList))
        break
    
    if blockdf.at[index,"contig"] != blockdf.at[idxList[count+1],"contig"]:
        contigList.append(blockdf.at[index,"contig"])
        tigSizeList.append(blockdf.at[index,"nchunk"])
        nblocksList.append(nblocks)
        avgBlockSizeList.append(blockSize/nblocks)
        for k in range(min(allChunksList), max(allChunksList)+1):
            if k not in allChunksList:
                missingChunksList.append(k)
        fullMissingChunksList.append(missingChunksList)
        nFullMissingChunksList.append(len(missingChunksList))
        nblocks = 0
        blockSize = 0
        allChunksList = []
        missingChunksList = []
        
tigdf["contig"] = contigList
tigdf["tig_size"] = tigSizeList
tigdf["nblocks"] = nblocksList
tigdf["avgblocksize"] = avgBlockSizeList
tigdf["missing_chunks"] = fullMissingChunksList
tigdf["n_missing_chunks"] = nFullMissingChunksList
tigdf["percent_missing_chunks"] = 100*tigdf["n_missing_chunks"]/tigdf["tig_size"]

tigdf.to_csv(sys.argv[2], index=False)
