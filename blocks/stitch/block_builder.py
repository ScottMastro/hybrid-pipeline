import block_builder_helper as helper
from intervals import Chunk
from intervals import Block
from intervals import Megablock
from intervals import Contig

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

def construct_blocks(chunks, param):
    ''' Takes a list of chunks and connects them into blocks.
    Valid blocks are returned, leftover chunks are returned in trash. '''
    blocks = []
    trash = []
    minPcId = param.CHUNK_MIN_IDENT

    while len(chunks) > 0:
        i = 0
        chunk = chunks.pop()
        block = Block(chunk)

        while i < len(chunks):

            i=i+1
            nextChunk = chunks[-i]
            idDist = chunk.id - nextChunk.id
            idSkip = idDist-1
            allowableGap = max(param.CHUNK_MAX_DIST, idSkip*param.CHUNK_SKIP_BP)
            bpDist = abs(chunk.mid(q=False) - nextChunk.mid(q=False)) - param.CHUNK_SIZE
            pcId = nextChunk.pcid
            
            #chunk id improper
            if idDist < 1: continue         
            if idSkip > param.CHUNK_SKIP: 
                break
            
            #bp gap too big
            if bpDist > allowableGap:
                continue
            #gap is in wrong direction
            strand = block.get_dir(q=False)
            bpDist = (chunk.mid(q=False) - nextChunk.mid(q=False)) - (param.CHUNK_SIZE * strand)
            
            if (bpDist*strand) < int(param.CHUNK_SIZE/2):
                continue
                
            #percent identity too low
            if pcId < minPcId:
                continue

            #chunk oriented incorrectly            
            if not block.rdir == nextChunk.get_dir(q=False):
                continue
            
            #all good, add chunk and remove from list
            block.add(nextChunk)
            chunk = chunks.pop(-i)
            i=i-1

            if len(chunks) == 0:
                break
            
        if(block.nchunk() > param.CHUNK_PER_BLOCK):
            blocks.append(block)
        else:
            trash.append(block)
    
    #trim overlapping regions in block
    blocks, trimmed1 = helper.trim_blocks(blocks, param, q=True)
    #todo: why do I have to trim twice????
    blocks, trimmed2 = helper.trim_blocks(blocks, param, q=True)
    blocks, trimmed3 = helper.trim_blocks(blocks, param, q=False)
    
    trash = trash + trimmed1 + trimmed2 + trimmed3
    
    blocks.sort()
    trash.sort()

    return (blocks, trash)


def construct_megablocks(blocks, length, param):
    
    if len(blocks) < 1: return []

    mblocks=[]
    
    #sort blocks by reference position
    blocks.sort(key=lambda block: block.left(q=False), reverse=True)
    currentBlock = blocks.pop()
    mblock = [currentBlock]
    
    maxrDist = param.MBLOCK_RDIST
    maxqDist = param.MBLOCK_QDIST
    percentLen = param.MBLOCK_QLEN_FACTOR
    
    while len(blocks) > 0:
        
        nextBlock = blocks.pop()

        rdist = nextBlock.left(q=False) - currentBlock.right(q=False)
        qdist = helper.min_dist(currentBlock, nextBlock)
        
        #break mblock formation
        if rdist > abs(maxrDist) or qdist > min(maxqDist, percentLen*length):
            mblocks.append(Megablock(mblock))
            mblock = []

        mblock.append(nextBlock)
        currentBlock = nextBlock
    
    mblocks.append(Megablock(mblock))

    #sort blocks by query position
    for mb in mblocks:
        mb.components.sort(key=lambda block: block.left(q=True))

    return mblocks

def construct_contig(mblockList, tigId, length, param, q=True):
    if len(mblockList) < 1: return None

    mblocks=[]

    #sort megablocks by query position
    mblockList.sort(key=lambda megablock: megablock.left(q=True), reverse=True)
    currentMegablock = mblockList.pop()
            
    overlapTolerance = 100
    minSpan = 5000

    skip=False
    
    #remove overlaps between megablocks
    #---------------------------------
    while len(mblockList) > 0:
        
        if not skip:
            nextMegablock = mblockList.pop()
        else: skip = False
                    
        pos1 = currentMegablock.right(q=True)
        pos2 = nextMegablock.left(q=True)
        
        if pos2 + overlapTolerance < pos1:
            
            cov1 = currentMegablock.coverage_between(pos1, pos2, q=True) #* \
               # math.log(currentMegablock.span_sum(q=True), 100)
            cov2 = nextMegablock.coverage_between(pos1, pos2, q=True) #* \
               # math.log(currentMegablock.span_sum(q=True), 100)
                
            if cov1 > cov2:
                nextMegablock.trim_left(pos1, q=True)
                if nextMegablock.span(q=True) < minSpan:
                    continue
            else:
                currentMegablock.trim_right(pos2, q=True)
                if currentMegablock.span(q=True) < minSpan:
                    if len(mblocks) > 0:
                        currentMegablock = mblocks.pop()
                        skip = True
                    else:
                        currentMegablock = nextMegablock
                    continue
            
        mblocks.append(currentMegablock)
        currentMegablock = nextMegablock
                            
    mblocks.append(currentMegablock)
    #---------------------------------


    mblockList=mblocks
    mblocks = []
    #mblockList.sort(key=lambda megablock: megablock.left(q=True), reverse=True)

    for mblock in mblockList:
        mblock.components.sort(key=lambda block: block.left(q=False))
        #if more than 2 blocks, we want to respect the order of the query
        if not mblock.is_consistent(q=False):
            mblock.components = mblock.components[::-1]
    
    maxDist=10000000

    megablock = mblockList.pop()
        
    while len(mblockList) > 0:
        
        nextMegablock = mblockList.pop()
    
        if megablock.rid != nextMegablock.rid:
            mblocks.append(megablock)
        else:
            rdist = nextMegablock.left(q=False) - megablock.right(q=False)
            if rdist > maxDist:
                mblocks.append(megablock)
            else:
                megablock.components = megablock.components + nextMegablock.components
                continue
            
        megablock = nextMegablock
                            
    mblocks.append(megablock)
        
    
    for mblock in mblocks:
        mblock.components.sort(key=lambda block: block.left(q=False))
        #if more than 2 blocks, we want to respect the order of the query
        if not mblock.is_consistent(q=False):
            mblock.components = mblock.components[::-1]

    contig = Contig(tigId, length, mblocks)

    return contig