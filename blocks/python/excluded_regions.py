class SimpleRegion:
    def __init__(self, chrom, start, end):     
        self.chrom = chrom
        self.start = start
        self.end = end
                            
    def __sub__(self, other):
        if other.chrom != self.chrom: return [self]
     
        if self.start < other.start:
            if self.end < other.start:
                return [self]
            elif self.end > other.end:
                sr1 = SimpleRegion(self.chrom, self.start, other.start-1)
                sr2 = SimpleRegion(self.chrom, other.end+1, self.end)
                return [sr1, sr2]
            elif self.end >= other.start:
                return [SimpleRegion(self.chrom, self.start, other.start-1)]
            
        if self.start >= other.start:
            if other.end < self.start:
                return [self]
            if self.end <= other.end:
                return []
            if other.end >= self.start:
                return [SimpleRegion(self.chrom, other.end+1, self.end)]
            
    def __len__(self):
        return self.end - self.start + 1
    def __repr__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp)"
    def __str__(self):
        return self.__repr__()
           
            
        
def region_difference(region, otherRegions, minSize=None):
    
    if minSize is not None and len(region) < minSize:
        return []
    
    if minSize is None: minSize = 2

    split = [region]
    
    for splitter in otherRegions:
        newSplit = []
        for s in split: newSplit.extend(s - splitter)
        split = newSplit
        
        if minSize is not None:
            split = [s for s in split if len(s) >= minSize]

    return split


def get_unused_contig(contig, lengthData, minSize=None, megablockLevel=True):

    tigId = contig.id
    fullRegion = SimpleRegion(tigId, 0, lengthData[tigId]-1) 
    
    usedRegions = []
    for mblock in contig.mblocks:
        if not megablockLevel:
                for block in mblock:
                    usedRegion = SimpleRegion(tigId, block.left(q=True), block.right(q=True)) 
                    usedRegions.append(usedRegion)            
        else:
            usedRegion = SimpleRegion(tigId, mblock.left(q=True), mblock.right(q=True)) 
            usedRegions.append(usedRegion)
            
    unusedRegions = region_difference(fullRegion, usedRegions, minSize)


     

def get_unused_regions(paths, tigIds, lengthData, minSize=10000):
    
    fullRegions = [ SimpleRegion(tigId, 0, lengthData[tigId]-1) for tigId in tigIds ]
    usedRegions = { tigId : [] for tigId in tigIds }
    
    def normalized_pos(fork, tigId):
        pos = fork.get_pos_by_id(tigId)
        if pos is None: return None
        if fork.get_strand_by_id(tigId) == -1:
            pos = lengthData[str(tigId)] - pos
        return pos

    
    for path in paths:
        if len(path) < 2: continue
    
        prevFork = path[0]
        for fork in path[1:]:
            if not fork.is_Nfork() and not prevFork.is_Nfork():
                
                tigId = prevFork.before_id()
                if (tigId in usedRegions) and fork.after_id() == tigId:
                    start = normalized_pos(prevFork, tigId)
                    end = normalized_pos(fork, tigId)
                    used = SimpleRegion(tigId, min(start, end), max(start, end))
                    usedRegions[tigId].append(used)
                    
                tigId = prevFork.after_id()
                if (tigId in usedRegions) and fork.before_id() == tigId:
                    start = normalized_pos(prevFork, tigId)
                    end = normalized_pos(fork, tigId)
                    used = SimpleRegion(tigId, min(start, end), max(start, end))
                    usedRegions[tigId].append(used)
                    
            prevFork = fork

    unusedRegions = { tigId : [] for tigId in tigIds }
    
    for region in fullRegions:        
        unusedRegions[region.chrom] = region_difference(region, usedRegions[region.chrom], minSize)
    
    return unusedRegions
    
    
    
    
    
    
    
    