import blocktree as bt

class Contig:
    mblocks=[]
    id=""
    size=-1
    def __init__(self, id, size, mblocks):
        self.id=id
        self.size=size

        self.mblocks = mblocks

    def is_empty(self): return len(self.mblocks) == 0
    def nblocks(self): return len(self.mblocks)
    def __repr__(self):
        return str(self.mblocks) 

    def __str__(self):
        return "qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                " query=" + str(self.left_q()) + "-" + str(self.right_q()) + \
                " reference=" + str(self.left_r()) + "-" + str(self.right_r()) + \
                " nblocks=" + str(self.nblocks())
    
def megablock_collapse(left, right, overlapTolerance=0.01):
   
    dist = right.data[0] - left.data[1]

    #overlap
    if dist < 0:
        smaller= min(left.data[3], right.data[3])
        
        if abs(dist/smaller) > overlapTolerance:
            return False
    
    return True
                
                
def construct_contig(id, size, mblocks, isQuery=True):
    
    if len(mblocks) < 1:
        return Contig(id, size, [])
    
    #constructs tree
    if isQuery:
        intervals=[(mblock.left_q(), mblock.right_q(), \
                mblock.cov_size_q(), mblock.size_q(), mblock) for mblock in mblocks]
    else:
        intervals=[(mblock.left_r(), mblock.right_r(), \
                mblock.cov_size_r(), mblock.size_r(), mblock) for mblock in mblocks]
    
    intervals.sort(key=lambda tup: tup[2])  # sorts by size

    tree = bt.Node(intervals.pop())
    while len(intervals) > 0:
        node = bt.Node(intervals.pop())
        tree.insert(node)
    
    collapseFn = lambda l, r : megablock_collapse(l, r, 0.01)

    mblocks, trash = bt.traverse(tree, None, None, collapseFn)
    
    
    if len(mblocks) < 1:
        return Contig(id, size, [])
    
    contig = Contig(id, size, mblocks)
    
    #todo: recover trash???
    
    return contig
