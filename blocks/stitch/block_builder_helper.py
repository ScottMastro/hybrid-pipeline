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
                blocks[i+1].trim_left(startPos, q)
                if len(blocks[i+1]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i+1)
                    trimmed.extend(trash.components)
                    continue
                    
            #right block is better
            else:
                blocks[i].trim_right(endPos, q)
                if len(blocks[i]) < param.CHUNK_PER_BLOCK:
                    trash = blocks.pop(i)
                    trimmed.extend(trash.components)
                    continue
            
        i = i + 1 

    return (blocks, trimmed)


def min_dist(A, B, q=True):
    dist1 = abs(B.left(q) - A.right(q))
    dist2 = abs(B.right(q) - A.right(q))
    dist3 = abs(B.left(q) - A.left(q))
    dist4 = abs(B.right(q) - A.left(q))
    return  min (dist1, dist2, dist3, dist4)

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