from stitch.intervals import Chunk
from stitch.intervals import construct_interval
from stitch.intervals import Contig

PASS, SKIP, FAIL = 0, 1, 2

#------- CHUNKS -----------------

def construct_chunks(row, param):
    '''Creates chunks from an aligndf row.'''
    qid=str(row[0])
    rid=str(row[2])
    cid=[int(x) for x in row[4].split(',')]
    rst=[int(x) for x in row[5].split(',')]
    red=[int(x) for x in row[6].split(',')]
    qst=[int(x) for x in row[7].split(',')]
    qed=[int(x) for x in row[8].split(',')]
    length=[int(x) for x in row[11].split(',')]
    identity=[float(x) for x in row[12].split(',')]

    #convert start/end position from chunk to contig
    qst=[(i-1)*param.CHUNK_SIZE + s for i, s in zip(cid, qst)]
    qed=[(i)*param.CHUNK_SIZE - (param.CHUNK_SIZE-e) for i, e in zip(cid, qed)]
    
    return [Chunk(cid[i], rst[i], red[i], qst[i], qed[i], \
                  identity[i], length[i], qid, rid) for i in range(0, len(cid))]


#------- BLOCKS -----------------
    

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


def valid_chunk_pair(chunk, nextChunk, param):
    ''' 
    Determines if two chunks are able to form a valid block together.
    Returns a value of PASS (pairs are valid), SKIP (pairs are not valid) or
    FAIL (pairs are not valid, stop attempting to extend block)
    '''

    idDist = chunk.id - nextChunk.id
    idSkip = idDist-1

    #chunk id improper
    if idDist < 1: return SKIP         
    if idSkip > param.CHUNK_SKIP: return FAIL
    
    bpDist = chunk.mid(q=False) - nextChunk.mid(q=False)
    allowableGap = max(param.CHUNK_MAX_DIST, idSkip*param.CHUNK_SKIP_BP)
    gap = abs(abs(bpDist) - param.CHUNK_SIZE)
    
    #bp gap too big
    if gap > allowableGap: return SKIP

    strand = chunk.get_dir(q=False)
    if strand * bpDist < 0:
        #gap is in wrong direction
        if gap > int(param.CHUNK_SIZE/4): return SKIP

        
    #bpDist = (chunk.mid(q=False) - nextChunk.mid(q=False)) - (param.CHUNK_SIZE * strand)
    #print(str(bpDist*strand) + " < " + str(int(param.CHUNK_SIZE/2)))
    #gap is in wrong direction
    #if (bpDist*strand) < int(param.CHUNK_SIZE/2): return SKIP
        
    #percent identity too low
    if nextChunk.pcid < param.CHUNK_MIN_IDENT: return SKIP

    #chunk oriented incorrectly            
    if not chunk.get_dir(q=False) == nextChunk.get_dir(q=False):
        return SKIP
    
    return PASS


def construct_blocks(chunks, param):
    '''
    Takes a list of chunks and connects them into blocks.
    Valid blocks are returned, leftover chunks are returned in trash.
    '''
    blocks = []
    trash = []

    while len(chunks) > 0:
        i = 0
        chunk = chunks.pop()
        block = construct_interval(chunk)

        while len(chunks) > 0 and i < len(chunks):

            i=i+1
            nextChunk = chunks[-i]
            
            result = valid_chunk_pair(chunk, nextChunk, param)
            if result == SKIP: continue
            if result == FAIL: break
            
            #all good, add chunk to block
            block.add(nextChunk)
            chunk = chunks.pop(-i)
            i=i-1
            
        if(len(block) > param.CHUNK_PER_BLOCK): 
            block.sort(q=True)
            blocks.append(block)
        else:
            trash.append(block)
    
    #trim overlapping regions in block
    blocks, t1 = trim_blocks(blocks, param, q=True)
    blocks, t2 = trim_blocks(blocks, param, q=False)
    
    trash = trash + t1 + t2
    blocks.sort(key=lambda b: b.left())

    return (blocks, trash)

#------- MEGABLOCKS -----------------

def valid_block_pair(block, nextBlock, length, param):
    ''' 
    Determines if two blocks are able to to form a megablock together.
    Returns a value of PASS (pairs are valid), SKIP (pairs are not valid) or
    FAIL (pairs are not valid, stop attempting to extend block)
    '''

    maxrDist = param.MBLOCK_RDIST
    maxqDist = param.MBLOCK_QDIST
    percentLen = param.MBLOCK_QLEN_FACTOR

    dist1 = abs(nextBlock.left(q=True) - block.right(q=True))
    dist2 = abs(nextBlock.right(q=True) - block.right(q=True))
    dist3 = abs(nextBlock.left(q=True) - block.left(q=True))
    dist4 = abs(nextBlock.right(q=True) - block.left(q=True))
    
    qdist = min (dist1, dist2, dist3, dist4)
    rdist = nextBlock.left(q=False) - block.right(q=False)

    #break megablock formation
    if rdist > abs(maxrDist) or qdist > min(maxqDist, percentLen*length):
        return FAIL

    return PASS


def construct_megablocks(blocks, length, param):
    ''' 
    Takes a list of blocks and connects them into megablocks.
    Valid megablocks are returned.
    '''

    if len(blocks) < 1: return []

    megablocks=[]
    
    #sort blocks by reference position
    blocks.sort(key=lambda block: block.left(q=False), reverse=True)
    block = blocks.pop()
    megablock = construct_interval(block)
    
    while len(blocks) > 0:
        
        nextBlock = blocks.pop()
        result = valid_block_pair(block, nextBlock, length, param)
        
        if result == FAIL:
            megablocks.append(megablock)
            megablock = construct_interval(nextBlock)
        elif result == PASS:
            megablock.add(nextBlock)
            
        block = nextBlock
    
    megablocks.append(megablock)

    #sort blocks by query position
    for mblock in megablocks:
        mblock.components.sort(key=lambda block: block.left(q=True))

    return megablocks


#------- CONTIGS -----------------

def trim_megablock_pair(mblockList, param):
    ''' 
    Detects overlapping megablocks and trims overlap if necessary.
    Returns trimmed megablocks.
    '''
    megablocks = []
    
    #sort megablocks by query position
    mblockList.sort(key=lambda megablock: megablock.left(q=True), reverse=True)
    megablock = mblockList.pop()

    overlapTolerance = param.MBLOCK_OVERLAP
    minSpan = param.MIN_MBLOCK_SIZE

    
    while len(mblockList) > 0:
        
        nextMegablock = mblockList.pop()

        pos1 = megablock.right(q=True)
        pos2 = nextMegablock.left(q=True)
            
        if pos2 + overlapTolerance < pos1:
            
            cov1 = megablock.coverage_between(pos1, pos2, q=True)
            cov2 = nextMegablock.coverage_between(pos1, pos2, q=True)
                
            if cov1 > cov2:
                nextMegablock.trim_left(pos1, q=True)
                if nextMegablock.span(q=True) < minSpan:
                    continue
            else:
                megablock.trim_right(pos2, q=True)
                if megablock.span(q=True) < minSpan:
                    if len(megablocks) > 0:
                        megablock = megablocks.pop()
                        mblockList.append(nextMegablock)
                    else:
                        megablock = nextMegablock
                    continue
            
        megablocks.append(megablock)
        megablock = nextMegablock
                            
    megablocks.append(megablock)
    return megablocks
    
def valid_megablock_pair(megablock, nextMegablock, param):
    ''' 
    Determines if two megablocks are able to to form a contig together.
    Returns a value of PASS (pairs are valid), SKIP (pairs are not valid) or
    FAIL (pairs are not valid, stop attempting to extend contig)
    '''

    maxDist = param.MAX_MBLOCK_DIST

    if megablock.rid != nextMegablock.rid:
        return FAIL
    
    rdist = nextMegablock.left(q=False) - megablock.right(q=False)
    if rdist > maxDist:
        return FAIL

    return PASS

def construct_contig(mblockList, tigId, length, param, q=True):
    ''' 
    Takes a list of megablocks and connects them into a Contig.
    A single Contig object is returned.
    '''

    if len(mblockList) < 1: return None

    mblockList = trim_megablock_pair(mblockList, param)
    
    #sort blocks by reference position
    for megablock in mblockList:
        megablock.components.sort(key=lambda block: block.left(q=False))
        if not megablock.is_consistent(q=False):
            megablock.components = megablock.components[::-1]

    contigBlocks = []
    megablock = mblockList.pop()
        
    while len(mblockList) > 0:
        
        nextMegablock = mblockList.pop()
        result = valid_megablock_pair(megablock, nextMegablock, param)
        if result == FAIL:
            contigBlocks.append(megablock)
        elif result == PASS:
            nextMegablock.components = megablock.components + nextMegablock.components

        megablock = nextMegablock
                            
    contigBlocks.append(megablock)
        
    #sort blocks by reference position
    for megablock in contigBlocks:
        megablock.components.sort(key=lambda block: block.left(q=False))
        #if more than 2 blocks, we want to respect the order of the query
        if not megablock.is_consistent(q=False):
            megablock.components = megablock.components[::-1]

    contig = Contig(tigId, length, contigBlocks)

    return contig