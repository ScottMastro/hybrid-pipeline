import helper
import external_tools as tools

class SimpleRegion:
    def __init__(self, chrom, start, end, lengthData=None):     
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = None
        if lengthData is not None:
            self.length = lengthData[chrom]
            
    def contains(self, chrom, pos):
        if chrom != self.chrom: return False
        if pos < self.start: return False
        if pos > self.end: return False
        return True
    
    def contains_pos(self, pos):
        if pos < self.start: return False
        if pos > self.end: return False
        return True

    def overlaps(self, otherRegion):
        if self.chrom != otherRegion.chrom: return False
        return max(self.start, otherRegion.start) <= min(self.end, otherRegion.end)
        
    def str_base1(self):
        self.start += 1
        self.end += 1
        string = str(self)
        self.start -= 1
        self.end -= 1
        return string

    def extend_left(self, size):
        newStart = max(0, self.start - size)
        return SimpleRegion(self.chrom, newStart, self.end)
        
    def extend_right(self, size):
        if self.length is None:
            newEnd = self.end + size
        else:
            newEnd = min(self.length-1, self.end + size)
        return SimpleRegion(self.chrom, self.start, newEnd)

    def extend(self, size):
        newStart = max(0, self.start - size)
        if self.length is None:
            newEnd = self.end + size
        else:
            newEnd = min(self.length-1, self.end + size)
        return SimpleRegion(self.chrom, newStart, newEnd)
        
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
    
    def __and__(self, other):
        if other.chrom != self.chrom: return []
        
        if (other.start > self.end) or (self.start > other.end):
            return []
        else:
            s = max(self.start, other.start)
            e = min(self.end, other.end)
            return [SimpleRegion(self.chrom, s, e)]

    def __eq__(self, other):
        return self.chrom == other.chrom and \
               self.start == other.start and \
               self.end == other.end
        
    def __len__(self):
        return self.end - self.start + 1
    def __repr__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp)"
    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)

def region_from_string(string, lengthData=None):
    chrom = string.split(":")[0]
    start = int(string.split(":")[1].split("-")[0])
    end = int(string.split(":")[1].split("-")[1])
    return SimpleRegion(chrom, start, end, lengthData)     


class PhasedRegion(SimpleRegion):
    def __init__(self, simpleRegion):   
        self.chrom = simpleRegion.chrom
        self.start = simpleRegion.start
        self.end = simpleRegion.end
        self.length = simpleRegion.length
        self.readsA, self.readsB, self.readsUnphased= [],[],[]
        self.variantsA, self.variantsB = None, None
        self.collapsedRegions = []
        
        self.fastaA, self.fastaB = None, None

    def add_fasta(self, fastaA=None, fastaB=None):
        self.fastaA = fastaA
        self.fastaB = fastaB
    
    def add_reads(self, readIdsA=None, readIdsB=None, readIdsUnphased=None):
        if readIdsA is not None:
            self.readsA = readIdsA
        if readIdsB is not None:
            self.readsB = readIdsB
        if readIdsUnphased is not None:
            self.readsUnphased = readIdsUnphased
        
    def add_collapsed_region(self, primaryRegion, collapsedRegion):
        self.collapsedRegions.extend(primaryRegion, collapsedRegion)

    def add_variants(self, variantSetA=None, variantSetB=None):
        if variantSetA is not None:
            self.variantsA = variantSetA
        if variantSetB is not None:
            self.variantsB = variantSetB
            
    def __repr__(self):
        info = []
        info.append(self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp)")
        
        info.append("readsA: " + str(len(self.readsA)) + \
                    ", readsB: " + str(len(self.readsB)) + \
                    ", readsUnphased: " + str(len(self.readsUnphased)))
        
        if self.variantsA is not None:
            info.append("A: " + self.variantsA.__repr__())
        if self.variantsB is not None:
            info.append("B: " + self.variantsB.__repr__())
   
        if self.fastaA is not None:
            info.append("fastaA: " + self.fastaA)
        if self.fastaB is not None:
            info.append("fastaB: " + self.fastaB)
        return "\n".join(info)


class GraphSegment(SimpleRegion):
    def __init__(self, regionA, regionB=None):
        
        self.regionA = regionA
        self.regionB = regionB
        
        self.chrom = regionA.chrom
        self.start = regionA.start
        self.end = regionA.end
        self.length = regionA.length

        self.seqA, self.seqB = None, None
        self.bamA, self.bamB, self.bamU = None, None, None
        self.variantsA, self.variantsB = None, None


    def write_seq(self, prefix, A=True, B=True):
        faA, faB = None, None
        
        if A and self.seqA is not None:
            faDictA = { str(self.regionA) + "_A": self.seqA }
            faA = tools.dict2fasta(faDictA, prefix + "_A")
            
        if B and self.seqB is not None:
            faDictB = { str(self.regionB) + "_B": self.seqB }
            faB = tools.dict2fasta(faDictB, prefix + "_B")
        
        return(faA, faB)
        
    def set_sequence(self, fastaA=None, fastaB=None):
        if fastaA is not None:
            self.seqA = helper.get_fasta_seq(fastaA)
        if fastaB is not None:
            self.seqB = helper.get_fasta_seq(fastaB)
        
        self.variantsA, self.variantsB = None, None
        self.collapsedRegions = []
        
        self.fastaA, self.fastaB = None, None
        
    def haplo_flip(self):
        self.seqA, self.seqB = self.seqB, self.seqA
        self.readsA, self.readsB = self.readsB, self.readsA
        self.regionA, self.regionB = self.regionB, self.regionA

    def set_reads(self, bamA=None, bamB=None, bamU=None):
        if bamA is not None:
            self.bamA = bamA
        if bamB is not None:
            self.bamB = bamB
        if bamU is not None:
            self.bamU = bamU
        
    def set_variants(self, variantSetA=None, variantSetB=None):
        if variantSetA is not None:
            self.variantsA = variantSetA
        if variantSetB is not None:
            self.variantsB = variantSetB
                
            
    def __repr__(self):
        info = []
        info.append(self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp)")
        
        info.append("readsA: " + ("" if self.bamA is None else self.bamA) + \
                    ", readsB: " + ("" if self.bamB is None else self.bamB) + \
                    ", readsUnphased: " + ("" if self.bamU is None else self.bamU))
        
        if self.variantsA is not None:
            info.append("A: " + self.variantsA.__repr__())
        if self.variantsB is not None:
            info.append("B: " + self.variantsB.__repr__())
   
        if self.seqA is not None:
            info.append("seqA: " + str(len(self.seqA)) + " bp")
        if self.seqB is not None:
            info.append("seqB: " + str(len(self.seqB)) + " bp")
        return "\n".join(info)


class AssemblyGraph():
    def __init__(self, tigName):   
        self.segments = []

    def add_segment(self, segment):
        self.regions.append(segment)
        
        
        
        
        
        
def region_difference(region, otherRegions, minSize=2):
    
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



def region_overlap(region, otherRegions, minOverlap=0):

    ovelaps = []
    
    for otherRegion in otherRegions:
        intersection = region & otherRegion
        
        if len(intersection) > 0:
            if len(intersection[0]) > minOverlap:
                ovelaps.append(otherRegion)
            
    return ovelaps



    '''

    #todo: region to sequence letters (Haplotype specific)      

    def summarize(self, printLines=True):
        lines = ["Summarizing region, " + self.__repr__()]
        
        A = self.alnA["qid"] + ":" + str(self.alnA["qstart"]) + "-" + str(self.alnA["qend"])
        
        B = None
        if self.alnB is not None:
            B = self.alnB["qid"] + ":" + str(self.alnB["qstart"]) + "-" + str(self.alnB["qend"])
        lines.append(A + " / " + B)

        for vsKey in self.variants:
            lines.append(vsKey + " variantset:")
            vsLines = self.variants[vsKey].summarize(printLines=False)
            vsLines = ['\t' + l for l in vsLines]
            lines.extend(vsLines)

        summary = '\n'.join(lines)
        if printLines: print(summary)
        
        return summary

    '''
     
