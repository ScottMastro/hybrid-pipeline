import blocktree as bt
import math
from intervals import Chunk
from intervals import Block
from intervals import Megablock
from intervals import Contig

def construct_chunks(row, param):
    '''
    Takes a row from the aligndf file, creates chunks.
    '''
    qid=row[0]
    rid=row[2]
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
    '''
    Takes a list of chunks and connects them into blocks.
    Valid blocks are returned, leftover chunks are returned in trash.
    '''
    blocks = []
    trash = []
    minPcId = param.MIN_PERCENT_ID

    while len(chunks) > 0:
         
        i = 0
        chunk = chunks.pop()
        block = Block(chunk)

        while i < len(chunks):

            i=i+1
            nextChunk = chunks[-i]
            idDist = chunk.id - nextChunk.id
            idSkip = idDist-1
            allowableGap = max(param.MAX_DIST_CHUNK, idSkip*param.CHUNK_SKIP_BP)
            bpDist = abs(chunk.rstart - nextChunk.rstart) - param.CHUNK_SIZE
            pcId = nextChunk.pcid
            
            #chunk id improper
            if idDist < 1: continue         
            if idSkip > param.CHUNK_SKIP: 
                break
            
            #bp gap too big
            if bpDist > allowableGap:
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
        
    blocks.sort()
    trash.sort()
    return (blocks, trash)


def get_dist(node1, node2, q=True):
    
    if (node1.left_pos(q) + node1.right_pos(q)) < \
    (node2.left_pos(q) + node2.right_pos(q)):
        leftNode, rightNode = node1, node2  
    else:
        rightNode, leftNode = node1, node2
   
    dist = rightNode.left_pos(q) - leftNode.right_pos(q)
    
    if node1.data.get_dir(q) != node2.data.get_dir(q):
        
        d2 = rightNode.right_pos(q) - leftNode.right_pos(q)
        d3 = rightNode.left_pos(q) - leftNode.left_pos(q)

        if abs(d2) < abs(dist): dist = d2
        if abs(d3) < abs(dist): dist = d3
        
    return dist
    

def can_collapse(node1, node2, distCutoff=1.5, pseudocountFactor=0.6, maxDist=200000, printout=False):
    pseudocount = max(1000, min(node1.span(), node2.span()) * pseudocountFactor)
   
    qdist = get_dist(node1, node2, q=True)

    if qdist > maxDist:
        if printout: print('hit max ' + str(qdist))
        return False
    
    rdist = get_dist(node1, node2, q=False)
        
    if rdist > maxDist: 
        if printout: print('hit max ' + str(rdist))
        return False
    
    if (qdist > 0 and rdist < 0) or (qdist < 0 and rdist > 0):
        num = abs(rdist) + abs(qdist) + pseudocount
        denom = pseudocount
    else:
        num = abs(qdist) + pseudocount
        denom = abs(rdist) + pseudocount
        
    if abs(rdist) < 1001 or abs(qdist) < 1001:
        distCutoff=distCutoff*2
        
    if printout and not abs(math.log(num*1.0/denom*1.0)) < math.log(distCutoff):
        print('=====')
        print("qdist=" + str(qdist))
        print("rdist=" + str(rdist))
        print(str(num) + "/" + str(denom) + "=" + str(num*1.0/denom*1.0))
        print(str(abs(math.log(num*1.0/denom*1.0))) + " vs " + str(math.log(distCutoff)))


    return abs(math.log(num*1.0/denom*1.0)) < abs(math.log(distCutoff))

def get_dist2(leftNode, rightNode, q=True):
    
    leftDir = leftNode.data.get_dir(q)
    rightDir = rightNode.data.get_dir(q)
    
    if leftDir == 1:
        if rightDir == 1:
            return rightNode.left_pos(q) - leftNode.right_pos(q)
        else:
            return rightNode.right_pos(q) - leftNode.right_pos(q)
    else:
        if rightDir == 1:
            return rightNode.left_pos(q) - leftNode.left_pos(q)
        else:
            return rightNode.left_pos(q) - leftNode.right_pos(q)

def get_dist_min(leftNode, rightNode, q=True):
    
    dist1 = abs(rightNode.left_pos(q) - leftNode.right_pos(q))
    dist2 = abs(rightNode.right_pos(q) - leftNode.right_pos(q))
    dist3 = abs(rightNode.left_pos(q) - leftNode.left_pos(q))
    dist4 = abs(rightNode.right_pos(q) - leftNode.left_pos(q))
    return min (dist1, dist2, dist3, dist4)
    

def can_collapse2(node1, node2, distCutoff=2.0, maxQueryDist=2000000, maxRefDist=50000):
   
    rdist = get_dist2(node1, node2, q=False)
    qdist = get_dist_min(node1, node2, q=True)
    
    if abs(qdist) > maxQueryDist or abs(rdist) > maxRefDist:
        return False
    
    ratio = (abs(qdist)+1000.0)/(abs(rdist)+1000)
    if abs(math.log(ratio)) < abs(math.log(distCutoff)):        
        return True
        
    return True

def construct_megablocks(blocks, length, param):

    if len(blocks) < 1: return []

    mblocks=[]
    
    while len(blocks) > 0:    
        #constructs tree, largest blocks at top
        sortedBlocks = sorted(blocks, key=lambda block: block.coverage())
        tree = bt.Node(sortedBlocks.pop())
        while len(sortedBlocks) > 0:
            node = bt.Node(sortedBlocks.pop(), q=False)
            tree.insert(node)
        
        #collapseFn = lambda l, r : can_collapse2(l, r, \
        #    param.BLOCK_DIST_THRESH, 0.6, 200000, printout=False)

        collapseFn = lambda l, r : can_collapse2(l, r)

        mblockNodes, mblockTrash = bt.traverse(tree, None, None, collapseFn)
        
        if len(mblockNodes) < 1:
            continue
        
        mblock = Megablock([node.data for node in mblockNodes ])
        if mblock.get_dir(q=True) == -1:
            mblock.components = mblock.components[::-1]
            
        mblocks.append(mblock)
        blocks = [node.data for node in mblockTrash ]
        
    return mblocks


def get_dist_block(leftBlock, rightBlock, q=True):
    
    leftDir = leftBlock.get_dir(q)
    rightDir = rightBlock.get_dir(q)
    
    return rightBlock.left(q) - leftBlock.right(q)

    
    if leftDir == 1:
        if rightDir == 1:
            return rightBlock.left(q) - leftBlock.right(q)
        else:
            return rightBlock.right(q) - leftBlock.right(q)
    else:
        if rightDir == 1:
            return rightBlock.left(q) - leftBlock.left(q)
        else:
            return rightBlock.left(q) - leftBlock.right(q)


def construct_megablocks2(blocks, length, param):
    if len(blocks) < 1: return []

    mblocks=[]
    
    blocks.sort(key=lambda block: block.left(q=False), reverse=True)
    block = blocks.pop()
    mblock = [block]
    
    maxrDist=20000
    maxqDist=200000
    percentLen = 0.05
    
    def add_mblock(mblock, nextBlock, reason=None):
        mblocks.append(Megablock(mblock))
        mblock = [nextBlock]
        
        #if reason is not None:
         #   print(mblocks[-1].__repr__())
           # print(reason)
        return mblock
    
    while len(blocks) > 0:
        
        nextBlock = blocks.pop()

        #print(block.__repr__() + " vs " + nextBlock.__repr__())

        rdist = get_dist_block(block, nextBlock, q=False)
        if rdist > abs(maxrDist):
            mblock = add_mblock(mblock, nextBlock, reason="rdist: " + str(rdist))
        else:
            dist1 = abs(nextBlock.left(q=True) - block.right(q=True))
            dist2 = abs(nextBlock.right(q=True) - block.right(q=True))
            dist3 = abs(nextBlock.left(q=True) - block.left(q=True))
            dist4 = abs(nextBlock.right(q=True) - block.left(q=True))
            qdist = min (dist1, dist2, dist3, dist4)
            
            if qdist > min(maxqDist, percentLen*length):
                mblock = add_mblock(mblock, nextBlock, reason="qdist: " + str(qdist))
            else:
                mblock.append(nextBlock)
        
        block = nextBlock
    
    mblocks.append(Megablock(mblock))

    for mb in mblocks:
        mb.components.sort(key=lambda block: block.left(q=True))

    return mblocks

def join_megablocks(mblockList, param):
    if len(mblockList) < 1: return []

    mblocks=[]
    
    mblockList.sort(key=lambda megablock: megablock.left(q=True), reverse=True)

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
            rdist = get_dist_block(megablock, nextMegablock, q=False)
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

    return mblocks

def clean_megablocks(mblockList, param):
    if len(mblockList) < 1: return []

    mblocks=[]
    
    mblockList.sort(key=lambda megablock: megablock.left(q=True), reverse=True)
    megablock = mblockList.pop()
            
    overlapTolerance = 100
    minSpan = 5000

    skip=False
    
    while len(mblockList) > 0:
        
        if not skip:
            nextMegablock = mblockList.pop()
        else: skip = False
            
        #print(str(megablock.__repr__()) + " vs " + str(nextMegablock.__repr__()))
        
        pos1 = megablock.right(q=True)
        pos2 = nextMegablock.left(q=True)
        
        if pos2 + overlapTolerance < pos1:
            
            cov1 = megablock.coverage_between(pos1, pos2, q=True) #* \
               # math.log(megablock.span_sum(q=True), 100)
            cov2 = nextMegablock.coverage_between(pos1, pos2, q=True) #* \
               # math.log(megablock.span_sum(q=True), 100)
                
            if cov1 > cov2:
                nextMegablock.trim_left(pos1, q=True)
                if nextMegablock.span(q=True) < minSpan:
                    continue
            else:
                megablock.trim_right(pos2, q=True)
                if megablock.span(q=True) < minSpan:
                    if len(mblocks) > 0:
                        megablock = mblocks.pop()
                        skip = True
                    else:
                        megablock = nextMegablock
                    continue
            
        mblocks.append(megablock)
        megablock = nextMegablock
                            
    mblocks.append(megablock)
            
    return join_megablocks(mblocks, param)   
            

def trim_overlap_query(blocks, param):

    i=0
    while ( i < len(blocks)-1 ):
        block = blocks[i]
        nextBlock=blocks[i+1]
        overlap = block.right_id() - nextBlock.left_id()
        #no overlap
        if overlap < 0:
            i = i + 1 
            continue
        
        #one repeated chunk
        if overlap == 0:
            overlapSize = nextBlock.right(q=True) - block.left(q=True)
            overlapStart = block.right(q=True)
            overlapEnd = nextBlock.left(q=True)
            if abs((overlapEnd - overlapStart) / overlapSize) < param.CHUNK_OVERLAP:
                i = i + 1 
                continue
        
        if overlap > 0:
            
            blockCov = block.coverage()
            nextBlockCov = nextBlock.coverage()
                        
            #left block is better
            if(blockCov > nextBlockCov):
                blocks[i+1].trim_left(block.right(q=True) )
                #blocks[i+1].components = nextBlock[overlap+1:]
                if len(blocks[i+1]) < param.CHUNK_PER_BLOCK:
                    blocks.pop(i+1)
            #right block is better
            else:
                blocks[i].trim_right(nextBlock.left(q=True) )
                #blocks[i].components = block[:-(overlap+1)]
                if len(blocks[i]) < param.CHUNK_PER_BLOCK:
                    blocks.pop(i)
            
    return blocks

def trim_blocks(blocks, param, q=True):

    i=0
    while ( i < len(blocks)-1 ):
        overlap = blocks[i].right(q) - blocks[i+1].left(q)

        if overlap > 0:

            blockCov = blocks[i].coverage()
            nextBlockCov = blocks[i+1].coverage()
                        
            #left block is better
            if(blockCov > nextBlockCov):
                blocks[i+1].trim_left(blocks[i].right(q), q)
                if len(blocks[i+1]) < param.CHUNK_PER_BLOCK:
                    blocks.pop(i+1)
                    continue
                    
            #right block is better
            else:
                blocks[i].trim_right(blocks[i+1].left(q), q)
                #blocks[i].components = block[:-(overlap+1)]
                if len(blocks[i]) < param.CHUNK_PER_BLOCK:
                    blocks.pop(i)
                    continue
            
        i = i + 1 

    return blocks

def remove_overlap(blocks, param):
    blocks.sort(key=lambda x: x.left(q=True))
    blocks = trim_blocks(blocks, param, q=True)
    blocks = trim_blocks(blocks, param, q=True)
    blocks.sort(key=lambda x: x.left(q=False))
    blocks = trim_blocks(blocks, param, q=False)
    blocks.sort()
    return blocks


def megablock_collapse(node1, node2, overlapTolerance=0.2, q=True):
    return True
    pos1 = node1.right_pos(q)
    pos2 = node2.left_pos(q)
    dist = pos2 - pos1 

    #overlap
    if dist < 0:
        smaller = min(node1.span(q), node2.span(q))

        if abs(dist/(smaller+1)) > overlapTolerance:
            return False
            
    return True
                
                
def construct_contig(tigId, size, mblockList, param, q=True):
    
    if len(mblockList) < 1:
        return Contig(tigId, size, [])
    
    sortedMblocks = sorted(mblockList, key=lambda mblock: mblock.coverage())
    tree = bt.Node(sortedMblocks.pop(), q)
    while len(sortedMblocks) > 0:
        node = bt.Node(sortedMblocks.pop(), q)
        tree.insert(node)
    
    collapseFn = lambda l, r : megablock_collapse(l, r, param.MEGABLOCK_OVERLAP, q)
    mblockNodes, trashNodes = bt.traverse(tree, None, None, collapseFn)
    
    if len(mblockNodes) < 1:
        return Contig(id, size, [])
    
    contig = Contig(id, size, [node.data for node in mblockNodes])

    #todo: recover trash???
    return contig

def trim_overlaps_contig(contig, param, q=True):
    i=0
    while i < len(contig.mblocks)-1:
        megablock = contig.mblocks[i]
        nextMegablock = contig.mblocks[i+1]
    
        pos1 = megablock.right(q)
        pos2 = nextMegablock.left(q)
        
        cov1 = megablock.coverage_between(pos1, pos2, q)
        cov2 = nextMegablock.coverage_between(pos1, pos2, q)
    
        if cov1 > cov2:
            nextMegablock.trim_left(pos1)
            if len(nextMegablock) < 1:
                contig.mblocks.pop(i+1)
                continue
        else:
            megablock.trim_right(pos2)
            if len(megablock) < 1:
                contig.mblocks.pop(i)
                continue
            
        i = i + 1
                
        
                
    return contig
