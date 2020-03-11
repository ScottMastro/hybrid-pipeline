class SimpleRegion:
    def __init__(self, chrom, start, end, annotation=""):     
        self.chrom = chrom
        self.start = start
        self.end = end
        self.annotation = annotation
        
    def contains(self, chrom, pos):
        if chrom != self.chrom: return False
        if pos < self.start: return False
        if pos > self.end: return False
        return True
        
    def __sub__(self, other):
        if other.chrom != self.chrom: return [self]
     
        if self.start < other.start:
            if self.end < other.start:
                return [self]
            elif self.end > other.end:
                sr1 = SimpleRegion(self.chrom, self.start, other.start-1, self.annotation)
                sr2 = SimpleRegion(self.chrom, other.end+1, self.end, self.annotation)
                return [sr1, sr2]
            elif self.end >= other.start:
                return [SimpleRegion(self.chrom, self.start, other.start-1, self.annotation)]
            
        if self.start >= other.start:
            if other.end < self.start:
                return [self]
            if self.end <= other.end:
                return []
            if other.end >= self.start:
                return [SimpleRegion(self.chrom, other.end+1, self.end, self.annotation)]
    
    def __and__(self, other):
        if other.chrom != self.chrom: return []
        
        if (other.start > self.end) or (self.start > other.end):
            return []
        else:
            s = max(self.start, other.start)
            e = min(self.end, other.end)
            return [SimpleRegion(self.chrom, s, e, self.annotation)]

    def __eq__(self, other):
        return self.chrom == other.chrom and \
               self.start == other.start and \
               self.end == other.end
        
    def __len__(self):
        return self.end - self.start + 1
    def __repr__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp) " + self.annotation
    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + " "  + self.annotation

def region_from_string(string, lengthData=None):
    chrom = string.split(":")[0]
    start = int(string.split(":")[1].split("-")[0])
    end = int(string.split(":")[1].split("-")[1])
    return SimpleRegion(chrom, start, end, lengthData)             
        
def region_difference(region, otherRegions, minSize=None):
    
    if minSize is not None and len(region) < minSize:
        return []

    split = [region]
    
    for splitter in otherRegions:
        newSplit = []
        for s in split: newSplit.extend(s - splitter)
        split = newSplit
        
        if minSize is not None:
            split = [s for s in split if len(s) >= minSize]

    return split

def regions_difference(regions, otherRegions, minSize=None):
    
    result = []
    for region in regions:
        result.extend(region_difference(region, otherRegions, minSize))
    return result

def region_overlap(region, otherRegions, minOverlap=0):

    overlaps = []
    
    for otherRegion in otherRegions:
        intersection = region & otherRegion
        
        if len(intersection) > 0:
            if len(intersection[0]) > minOverlap:
                overlaps.append(otherRegion)
            
    return overlaps

def region_intersection(region, otherRegions, minOverlap=0):

    intersections = []
    
    for otherRegion in otherRegions:
        intersection = region & otherRegion
        
        if len(intersection) > 0:
            for inters in intersection:
                if len(inters) > minOverlap: intersections.append(inters)
            
    return intersections


def regions_intersection(regions, otherRegions, minOverlap=0):
    
    result = []
    for region in regions:
        result.extend(region_intersection(region, otherRegions, minOverlap))
    return result
