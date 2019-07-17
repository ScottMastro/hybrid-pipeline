import log
import block_builder as blocker
import plot_blocks
import copy

def stitch(tigId, aligndf, lengthData, param, plotDir=None, q=True):    
    
    rows = aligndf[aligndf['qid' if q else 'rid'].astype(str) == tigId]
    
    if len(rows) == 0:
        log.out("Contig " + str(tigId) + " has no chunk alignments.", 1, param)
        return None
    
    clist, blist, mlist = [], [], []

    mblockList=[]
    length = lengthData[tigId]
    
    for idx, row in rows.iterrows():
        #if row[2] == 'tig00007745_pilon_pilon': break
                
        #chunks
        chunks = blocker.construct_chunks(row, param)
        if plotDir is not None: clist.extend(copy.deepcopy(chunks))
        
        #blocks
        blocks, trash = blocker.construct_blocks(chunks, param)
        if plotDir is not None: blist.extend(copy.deepcopy(blocks))

        #megablock
        mblocks = blocker.construct_megablocks(blocks, length, param)
        if plotDir is not None: mlist.extend(copy.deepcopy(mblocks))
        mblockList.extend(mblocks)

    #contig
    contig = blocker.construct_contig(mblockList, tigId, length, param)
    tlist = [] if contig is None else contig.mblocks
    if not q: tlist = []
    
    if len(clist) > 0 and plotDir is not None:
        intervalList=[clist, blist, mlist, tlist]
        plot_blocks.plot_levels(intervalList, length, q=q, outputPath=plotDir + tigId)

    return contig


        