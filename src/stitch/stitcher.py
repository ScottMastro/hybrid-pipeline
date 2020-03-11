import sys
sys.path.append("..")
import log as logger
import stitch.build_intervals as builder

def stitch_blocks(tigId, alignDict, param):    
    if tigId not in alignDict:
        logger.out("Contig " + str(tigId) + " has no chunk alignments.", 1, param)
        return []
    
    rows = alignDict[tigId]
    blockList = []
    for idx, row in rows.iterrows():
        
        #chunks
        chunks = builder.construct_chunks(row, param)
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        blockList.extend(blocks)

    return blockList


def stitch_megablocks(tigId, blocks, lengthData, param):    
    
    mblockList=[]
    length = lengthData[tigId]
    rids = set([block.rid for block in blocks])
    
    for rid in rids:
        rblocks = [block for block in blocks if block.rid == rid]
        megablocks = builder.construct_megablocks(rblocks, length, param)
        mblockList.extend(megablocks)

    return mblockList


def stitch_contig(tigId, mblocks, lengthData, param):    
    
    length = lengthData[tigId]
    contig = builder.construct_contig(mblocks, tigId, length, param)

    return contig


        