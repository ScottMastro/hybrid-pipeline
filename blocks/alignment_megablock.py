import uuid
import blocktree as bt
from math import log

class Megablock:
    blocks=[]
    rid=""
    qid=""
  
    def __init__(self, blocks):
        self.index=uuid.uuid4()
        if len(blocks) > 0:
            self.rid=blocks[0].rid
            self.qid=blocks[0].qid
        self.blocks = blocks
        
    def left_q(self): return self.blocks[0].left_q()
    def right_q(self): return self.blocks[-1].right_q()
    def left_r(self): return self.blocks[0].right_q()
    def right_r(self): return self.blocks[0].left_q()

    def cov_size_q(self): return sum([block.size_q() for block in self.blocks])
    def cov_size_r(self): return sum([block.size_r() for block in self.blocks])    
    def size_q(self): return abs(self.left_q() - self.right_q())
    def size_r(self): return abs(self.left_r() - self.right_r())

    def nblocks(self): return len(self.blocks)
    def __repr__(self):
        return str(self.blocks) 

    def __str__(self):
        return "qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                " query=" + str(self.left_q()) + "-" + str(self.right_q()) + \
                " reference=" + str(self.left_r()) + "-" + str(self.right_r()) + \
                " nblocks=" + str(self.nblocks())


def can_collapse(left, right, distCutoff=1.5, pseudocount=0):
    qdist = right.data[0] - left.data[1]
    rdist = right.data[2] - left.data[3]
    
    qoverlap = min(0, qdist)
    if qoverlap != 0: qdist = 0
    
    #overlap = total size - size of prevNode - size of node
    size = max(left.data[2], left.data[3], right.data[2], right.data[3]) - \
           min(left.data[2], left.data[3], right.data[2], right.data[3])
    roverlap = size - abs(left.data[3] - left.data[2]) - \
           abs(right.data[3] - right.data[2])       

    roverlap = min(0, roverlap)
    if roverlap != 0: rdist = 0

    if qoverlap == 0:
        if roverlap == 0:
            num = abs(qdist) + pseudocount
            denom = abs(rdist) + pseudocount
        else:
            num = abs(qdist) + pseudocount
            denom = abs(roverlap) + abs(qdist) + pseudocount
    elif roverlap == 0:
        num = abs(qoverlap) + abs(rdist) + pseudocount
        denom = abs(rdist) + pseudocount
    else:
        num = abs(qoverlap) + pseudocount
        denom = abs(roverlap) + pseudocount

    if abs(log(num/denom)) < abs(log(distCutoff)):
        return True
    else:
        return False    

def construct_megablocks(row, blocks, trash, chunkSize=1000, distCutoff=1.5):

    if len(blocks) < 1:
        return []

    #constructs tree
    intervals=[(block.left_q(), block.right_q(), \
            block.left_r(), block.right_r(), \
            block.size_q(), block) for block in blocks] 
    intervals.sort(key=lambda tup: tup[4])  # sorts by size of query

    tree = bt.Node(intervals.pop())
    while len(intervals) > 0:
        node = bt.Node(intervals.pop())
        tree.insert(node)
    
    
    length = row[1]*chunkSize
    collapseFn = lambda l, r : can_collapse(l, r, distCutoff, length*0.01)
    
    blocksx, trash = bt.traverse(tree, None, None, collapseFn)
    
    if len(blocks) < 1:
        return []
    
    mblock = Megablock(blocks)
    #todo: recover trash???
    return [mblock]
