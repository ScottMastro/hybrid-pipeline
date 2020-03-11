import sys

sys.path.append("..")

from stitch.intervals import Block
from stitch.intervals import Megablock
from stitch.intervals import Contig
import regions as rgn

'''
def stitch(tigId, alignDict, lengthData, param, q=True):    
    """Performs all stitching tasks: creating chunks, blocks and megablocks.
    tigId is the ID of the assembled contig to stitch.
    alignDict is the dictionary containing chunk alignments between assemblies.
    lengthData is a dictionary containing the lengths of all assembled contigs.
    param is the Parameters object.
    q is True if tigId is a query contig, False if reference contig.
    
    q=True implies chunks, blocks, etc will be formed along query contig.
    q=False is only intended to plot from the perspective of the reference.

    If param.PLOT_BLOCKS is not None, Contig object will be plotted.

    Returns a Contig object that contains all constructed megablocks.
    Returns None if no chunks alignments exist."""

    plotDir = param.PLOT_BLOCKS
    
    if tigId not in alignDict:
        logger.out("Contig " + str(tigId) + " has no chunk alignments.", 1, param)
        return None
    
    rows = alignDict[tigId]

    #save Intervals for plotting
    cList, bList, mList = [], [], []

    mblockList=[]
    length = lengthData[tigId]
    
    for idx, row in rows.iterrows():
        #if row[2] == 'tig00007745_pilon_pilon': break

     
        #chunks
        chunks = builder.construct_chunks(row, param)
        if plotDir is not None: cList.extend(copy.deepcopy(chunks))
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        if plotDir is not None: bList.extend(copy.deepcopy(blocks))

        #megablock
        megablocks = builder.construct_megablocks(blocks, length, param)
        mblockList.extend(megablocks)
        if plotDir is not None: mList.extend(copy.deepcopy(megablocks))

    #contig
    contig = builder.construct_contig(mblockList, tigId, length, param)
    
    #plot
    if len(cList) > 0 and plotDir is not None:
        
        tList = [] if contig is None else contig.mblocks
        if not q: tList = []
        intervalList=[cList, bList, mList, tList]
        plotter.plot_levels(intervalList, length, q=q, outputPath=plotDir + tigId)

    return contig

'''

def trim_blocks(blocks, param, q=True):
    '''Removes and returns chunks from overlapping blocks.'''
    
    blocks.sort(key=lambda x: x.left(q))
    trimmed = []
    
    i=0
    while ( i < len(blocks)-1 ):
        startPos = blocks[i].right(q)
        endPos =  blocks[i+1].left(q)

        if startPos > endPos:

            blockCov = blocks[i].coverage_between(startPos, endPos, q)
            nextBlockCov = blocks[i+1].coverage_between(startPos, endPos, q)
                        
            #left block is better
            if(blockCov > nextBlockCov):
                tl = blocks[i+1].trim_left(startPos, q)
                trimmed.extend(tl)
                if len(blocks[i+1]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i+1)
                    trimmed.extend(trash.components)
                    continue
                    
            #right block is better
            else:
                tr = blocks[i].trim_right(endPos, q)
                trimmed.extend(tr)
                if len(blocks[i]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i)
                    trimmed.extend(trash.components)
                    continue
            
        i = i + 1 

    return (blocks, trimmed)




def construct_block(alignment, lengthData):
    '''
    Converts a PAF alignment into a block.
    '''
    
    qregion = rgn.SimpleRegion(alignment["qid"], alignment["qstart"], alignment["qend"], lengthData)
    rregion = rgn.SimpleRegion(alignment["rid"], alignment["rstart"], alignment["rend"], lengthData)
    strand = alignment["strand"]
    strand = -1 if strand == "-" else strand
    strand = 1 if strand == "+" else strand
    strand = int(strand)
    
    alignLen = alignment["alnlen"]
    percentIdentity = 1 - alignment["dv"]
    
    if "cg" in alignment:
        cigarString = alignment["cg"]
    else:
        cigarString = None 
        
    return Block(qregion, rregion, strand, percentIdentity, alignLen, cigar=cigarString)


def valid_block_pair(block, nextBlock, param):
    ''' 
    Determines if two blocks are able to to form a megablock together.
    Returns a value of PASS (pairs are valid), SKIP (pairs are not valid) or
    FAIL (pairs are not valid, stop attempting to extend block)
    '''

    if block.rid != nextBlock.rid:
        return False

    maxrDist = param.MBLOCK_RDIST
    maxqDist = param.MBLOCK_QDIST
    percentLen = param.MBLOCK_QLEN_FACTOR

    dist1 = abs(nextBlock.start(q=True) - block.end(q=True))
    dist2 = abs(nextBlock.end(q=True) - block.end(q=True))
    dist3 = abs(nextBlock.start(q=True) - block.start(q=True))
    dist4 = abs(nextBlock.end(q=True) - block.start(q=True))
    
    qdist = min (dist1, dist2, dist3, dist4)
    rdist = nextBlock.start(q=False) - block.end(q=False)
    
    #break megablock formation
    if rdist > abs(maxrDist) or qdist > maxqDist:
        print(block, nextBlock)
        print(rdist, qdist)
        return False

    return True


def construct_megablocks(blocks, length, param):
    ''' 
    Takes a list of blocks and connects them into megablocks.
    Valid megablocks are returned.
    '''

    if len(blocks) < 1: return []

    mblocks=[]
    
    #sort blocks by reference position
    blocks.sort(key=lambda block: (block.qid, block.start(q=False)), reverse=True)
    block = blocks.pop()
    mblock = Megablock(block)
    
    while len(blocks) > 0:
        
        nextBlock = blocks.pop()
        result = valid_block_pair(block, nextBlock, param)
        
        if result:
            mblock.add(nextBlock)
        else:
            mblocks.append(mblock)
            mblock = Megablock(nextBlock)
            
        block = nextBlock
    
    mblocks.append(mblock)

    #sort blocks by query position
    for mblock in megablocks:
        mblock.components.sort(key=lambda block: block.left(q=True))

    return megablocks


def trim_blocks(blocks, param):
    
    blocks.sort(key=lambda x: x.start(q=False))
    trimmed = []
    
    i=0
    while ( i < len(blocks)-1 ):
        startPos = blocks[i].end(q=False)
        endPos =  blocks[i+1].start(q=False)

        if startPos > endPos:

            blockCov = blocks[i].coverage_between(startPos, endPos, q)
            nextBlockCov = blocks[i+1].coverage_between(startPos, endPos, q)
                        
            #left block is better
            if(blockCov > nextBlockCov):
                tl = blocks[i+1].trim_left(startPos, q)
                trimmed.extend(tl)
                if len(blocks[i+1]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i+1)
                    trimmed.extend(trash.components)
                    continue
                    
            #right block is better
            else:
                tr = blocks[i].trim_right(endPos, q)
                trimmed.extend(tr)
                if len(blocks[i]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i)
                    trimmed.extend(trash.components)
                    continue
            
        i = i + 1 

    return (blocks, trimmed)


def stitch2(qids, aligndf, lengthData, param):
    
    minAlignLen = 1000
    
    blockDict = { qid : [] for qid in qids }

    for _,alignment in aligndf.iterrows():
        if alignment["alnlen"] < minAlignLen:
            continue
        
        block = construct_block(alignment, lengthData)
        blockDict[block.qid].append(block)
    
    
    
    mblocks = { qid : [] for qid in qids }
    for qid in qids:
        blocks = blockDict[qid]
        
        ridSet = { block.rid for block in blocks }
        
        for rid in ridSet:
            rblocks = [block for block in blocks if block.rid == rid]
            mblocks = construct_megablocks(rblocks)

            
            
        
        lenDict = { block.rid : 0 for block in blocks }
        for block in blocks: lenDict[block.rid] += len(block)
        
        
        
        mblock = construct_megablocks.add(blocks[qid])
        mblocks[mblock.qid].append(mblock)
    
    
    contigs = []
    for qid in qids:
        contig = Contig()
        for block in blockDict[qid]:
            contig.add(block)
            
        contig.sort()
        contigs.append(contig)
    
    
    
    
        
        
        
        
        rows = alignDict[tigId]

    #save Intervals for plotting
    cList, bList, mList = [], [], []

    mblockList=[]
    length = lengthData[tigId]
    
    for idx, row in rows.iterrows():
        #if row[2] == 'tig00007745_pilon_pilon': break

        '''
        idx=0
        for i in range(idx): row = rows.iterrows().send(None)        
        '''       
        #chunks
        chunks = builder.construct_chunks(row, param)
        if plotDir is not None: cList.extend(copy.deepcopy(chunks))
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        if plotDir is not None: bList.extend(copy.deepcopy(blocks))

        #megablock
        megablocks = builder.construct_megablocks(blocks, length, param)
        mblockList.extend(megablocks)
        if plotDir is not None: mList.extend(copy.deepcopy(megablocks))

    #contig
    contig = builder.construct_contig(mblockList, tigId, length, param)
    
    #plot
    if len(cList) > 0 and plotDir is not None:
        
        tList = [] if contig is None else contig.mblocks
        if not q: tList = []
        intervalList=[cList, bList, mList, tList]
        plotter.plot_levels(intervalList, length, q=q, outputPath=plotDir + tigId)

    return contig

