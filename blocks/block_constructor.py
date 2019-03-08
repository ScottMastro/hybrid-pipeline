import blocktree as bt
from math import log
from intervals import Chunk
from intervals import Block
from intervals import Megablock
from intervals import Contig

#todo:should we do this???
def remove_redundancy(blocks, chunkPerBlock, repeatOverlap):

    i=0
    while ( i < len(blocks)-1 ):
        
        x = blocks[i+1].left_id() - blocks[i].right_id()
        
        #no overlap
        if x > 0:
            i = i + 1 
            continue
        
        #overlapping
        if x <= 0:
            #one repeated chunk
            if x == 0:
                overlapSize = blocks[i+1].right(q=True) - blocks[i].left(q=True)
                overlapStart = blocks[i].right(q=True)
                overlapEnd = blocks[i+1].left(q=True)
                if abs((overlapEnd - overlapStart) / overlapSize) < repeatOverlap:
                    i = i + 1 
                    continue
            
            #left block is bigger
            if(blocks[i].nchunk() > blocks[i+1].nchunk()):
                blocks[i+1].trim_left(-x + 1)
                if(blocks[i+1].nchunk() < chunkPerBlock):
                    blocks.pop(i+1)
            #right block is bigger
            else:
                blocks[i].trim_right(-x + 1)
                if(blocks[i].nchunk() < chunkPerBlock):
                    blocks.pop(i)
            
    return blocks

def construct_chunks(row, chunkSize=1000):
    qid=row[0]
    rid=row[2]
    cid=[int(x) for x in row[4].split(',')]
    rst=[int(x) for x in row[5].split(',')]
    red=[int(x) for x in row[6].split(',')]
    qst=[int(x) for x in row[7].split(',')]
    qed=[int(x) for x in row[8].split(',')]
    
    #convert start/end position from chunk to contig
    qst=[(i-1)*chunkSize + s for i, s in zip(cid, qst)]
    qed=[(i)*chunkSize - (chunkSize-e) for i, e in zip(cid, qed)]

    return [Chunk(cid[i], rst[i], red[i], qst[i], qed[i], qid, rid) for i in range(0, len(cid))]


def construct_blocks(chunks, chunkSize=1000, chunkPerBlock=3, maxDist=5000, chunkSkip=2, bpSkip=2000):
    '''
    Takes a row from the aligndf file, creates chunks and connects them into blocks.
    Valid blocks are returned, leftover chunks are returned in trash.
    '''
    blocks = []
    trash = []

    while len(chunks) > 0:
         
        i = 0
        chunk = chunks.pop()
        block = Block(chunk)
                
        while i < len(chunks):

            i=i+1
            nextChunk = chunks[-i]
            idDist = chunk.id - nextChunk.id
            idSkip = idDist-1
            allowableGap = max(maxDist, idSkip*bpSkip)
            bpDist = abs(chunk.rstart - nextChunk.rstart) - chunkSize

            #chunk id improper
            if idDist < 1: continue         
            if idSkip > chunkSkip: break
            
            #bp gap too big
            if bpDist > allowableGap:
                continue
            
            #chunk oriented incorrectly
            if not block.verify_direction(nextChunk):
                continue
            
            #all good, add chunk and remove from list
            block.add(nextChunk)
            chunk = chunks.pop(-i)
            i=i-1

            if len(chunks) == 0:
                break
            
            
        if(block.nchunk() > chunkPerBlock):
            blocks.append(block)
        else:
            trash.append(block)
        
    blocks.sort()
    trash.sort()
    return (blocks, trash)

def can_collapse(node1, node2, distCutoff=1.5, pseudocount=0, printout=False):
    
    if(printout):
        print("=====================\nleft=" + node1.__repr__())
        print("right=" + node2.__repr__())

    q=True
    if (node1.left_pos(q) + node1.right_pos(q)) < \
        (node2.left_pos(q) + node2.right_pos(q)):
            leftNode, rightNode = node1, node2
    else:
        leftNode, rightNode = node2, node1
        
    qdist = rightNode.left_pos(q) - leftNode.right_pos(q)

    q=False

    if (node1.left_pos(q) + node1.right_pos(q)) < \
        (node2.left_pos(q) + node2.right_pos(q)):
            leftNode, rightNode = node1, node2
    else:
        leftNode, rightNode = node2, node1
        
    rdist = rightNode.left_pos(q) - leftNode.right_pos(q)
    
    if qdist >= 0:
        if rdist >= 0:
            #no overlap
            num = abs(qdist) + pseudocount
            denom = abs(rdist) + pseudocount
        else:
            #ref overlaps
            num = abs(qdist) + pseudocount
            denom = abs(rdist) + abs(qdist) + pseudocount
    elif rdist >= 0:
        #query overlaps
        num = abs(qdist) + abs(rdist) + pseudocount
        denom = abs(rdist) + pseudocount
    else:
        #both overlap
        num = abs(qdist) + pseudocount
        denom = abs(qdist) + pseudocount

    if printout:
        print("qdist=" + str(qdist))
        print("rdist=" + str(rdist))
        print(str(num) + "/" + str(denom) + "=" + str(num*1.0/denom*1.0))
        print(str(abs(log(num*1.0/denom*1.0))) + " vs " + str(log(distCutoff)))

    return abs(log(num*1.0/denom*1.0)) < abs(log(distCutoff))
        
def construct_megablocks(row, blocks, trash, length):

    if len(blocks) < 1:
        return []

    q=True
    mblocks=[]
    
    while len(blocks) > 0:    
        #constructs tree, largest blocks at top
        sortedBlocks = sorted(blocks, key=lambda block: block.sum(q))
        tree = bt.Node(sortedBlocks.pop())
        while len(sortedBlocks) > 0:
            node = bt.Node(sortedBlocks.pop())
            tree.insert(node)
        
        collapseFn = lambda l, r : can_collapse(l, r, 1.5, length*0.05, printout=False)
        
        mblockNodes, mblockTrash = bt.traverse(tree, None, None, collapseFn)
        
        if len(mblockNodes) < 1:
            continue
        
        mblock = Megablock([node.data for node in mblockNodes ])    
        mblocks.append(mblock)
        blocks = [node.data for node in mblockTrash ]
        
    return mblocks


def megablock_collapse(node1, node2, overlapTolerance=0.01):
   
    dist = node2.left_pos() - node1.right_pos() 

    #overlap
    if dist < 0:
        smaller= min(node1.span(), node2.span())
        
        if abs(dist/smaller) > overlapTolerance:
            return False
    
    return True
                
                
def construct_contig(id, size, mblockList, q=True):
    
    mblocks = []
    for mblockSublist in mblockList:
        mblocks.extend(mblockSublist)
    
    if len(mblocks) < 1:
        return Contig(id, size, [])
    
    sortedMblocks = sorted(mblocks, key=lambda mblock: mblock.sum(q))
    tree = bt.Node(sortedMblocks.pop())
    while len(sortedMblocks) > 0:
        node = bt.Node(sortedMblocks.pop())
        tree.insert(node)
    
    collapseFn = lambda l, r : megablock_collapse(l, r, 0.01)
    mblockNodes, trashNodes = bt.traverse(tree, None, None, collapseFn)
    
    if len(mblockNodes) < 1:
        return Contig(id, size, [])
    
    contig = Contig(id, size, [node.data for node in mblockNodes])
    #todo: recover trash???
    
    return contig