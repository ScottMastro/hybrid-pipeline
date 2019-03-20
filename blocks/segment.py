reverseMap = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N', '-':'-'} 

def reverse_complement(seq):
    return ''.join([reverseMap[x] for x in seq[::-1]])


class Segment:
    def __init__(self, tigId, startPos, startStrand, endPos, endStrand, length):
        self.id=tigId
        self.startPos=startPos
        self.endPos=endPos
        self.startStrand=startStrand
        self.endStrand=endStrand
        self.length=length

    #these may give nonsense values if strands are not valid
    def left(self):
        return self.startPos if self.startPos < self.endPos else self.endPos
    def right(self):
        return self.endPos if self.startPos < self.endPos else self.startPos
    def start(self):
        return self.startPos
    def end(self):
        return self.endPos

    def __len__(self):
        return abs(self.endPos - self.startPos)
    
    def start_strand(self, string=False):
        if self.startStrand == 1:
            return "+" if string else 1
        if self.startStrand == -1:
            return "-" if string else -1
        else:
            return "?" if string else 0
        
    def end_strand(self, string=False):
        if self.endStrand == 1:
            return "+" if string else 1
        if self.endStrand == -1:
            return "-" if string else -1
        else:
            return "?" if string else 0
        
    def __repr__(self):
        return str(self.id) + ":" + \
            str(self.startPos) +  "(" + self.start_strand(string=True)  + ")" +\
            "-" + str(self.endPos) +  "(" + self.end_strand(string=True)  + ")"


    def __str__(self):
        return str(self.id) + ":" + \
            str(self.startPos) +  "(" + self.start_strand(string=True)  + ")" +\
            "-" + str(self.endPos) +  "(" + self.end_strand(string=True)  + ")"
            
    def get_sequence(self, seq):
        if not self.is_valid():
            return ""
        
        if self.startStrand == -1:
            reverseSelf = self.make_reverse_complement(getCopy=True)
            s = reverseSelf.get_sequence(seq)
            s = reverse_complement(s)
            return s
        
        return seq[self.left():self.right()]
    
    def get_validation_sequence(self, seq, side, nbases=50):
        if not self.is_valid():
            return ""

        if self.startStrand == -1:
            reverseSelf = self.make_reverse_complement(getCopy=True)
            flipSide = "l" if side == "r" else "r"
            s = reverseSelf.get_validation_sequence(seq, flipSide, nbases=50)
            s = reverse_complement(s)
            return s
        
        if side == 'l':
            return seq[self.left() - int(nbases/2) : self.left() + int(nbases/2)]
        if side == 'r':
            return seq[self.right() - int(nbases/2) : self.right() + int(nbases/2) ] 

    def valid_strands(self):
        return self.startStrand == self.endStrand
    
    def flip_end(self):
        self.endPos = self.length - self.endPos
        self.endStrand = self.endStrand*(-1)


    def is_valid(self):
        if self.startStrand != self.endStrand:
            return False
        if self.startPos == self.endPos:
            return False
        if self.startPos > self.endPos:
            return False
        return True
    
    def try_merge(self, other):
        '''assumes self starts before other'''
        
        if not self.is_valid() or not other.is_valid():
            return None
        
        if not self.id == other.id or not self.startStrand == other.startStrand:
            return None
        
        isOverlapping = self.right() >= other.left() and self.left() <= other.right()
        if not isOverlapping:
            return None
        
        minPos = min(self.startPos, self.endPos, other.startPos, other.endPos)
        maxPos = max(self.startPos, self.endPos, other.startPos, other.endPos)

        if self.startStrand == 1:
            return Segment(self.id, minPos, 1, maxPos, 1, self.length)
        elif self.startStrand == -1:
            return Segment(self.id, maxPos, -1, minPos, -1, self.length)

        return None     
    
    def make_reverse_complement(self, getCopy=False):
        #todo: might be off-by-one error here????
        startPosrc = self.length - self.endPos
        endPosrc = self.length - self.startPos
        startStrandrc = self.endStrand*(-1)
        endStrandrc = self.startStrand*(-1)

        if getCopy:
            return Segment(self.id, startPosrc, startStrandrc, endPosrc, endStrandrc, self.length)
        else:
            self.startPos = startPosrc
            self.endPos = endPosrc
            self.startStrand = startStrandrc
            self.endStrand = endStrandrc

            return self
        
        '''    
        def get_segment_between(self, other):
        #assumes self is left of other
        
        if self.dir != other.dir:
            print("improper direction")
            return None
        
        if self.id != other.id:
            print("improper id")
            return None
        
        if self.length != other.length:
            print("improper length")
            return None
        
        newSeg = Segment(self.id, self.endPos, other.startPos, self.dir, self.length)
        
        return newSeg        

    def get_sequence_between(self, other, seq):
        #assumes self is left of other
        
        if self.dir != other.dir:
            print("improper direction")
            return "___"
        
        if self.dir == -1:
            reverseSelf = self.make_reverse_complement(getCopy=True)
            reverseOther = other.make_reverse_complement(getCopy=True)

            s = reverseSelf.get_sequence_between(reverseOther, seq)
            s = reverse_complement(s)
            return s
        
        l = min(self.endPos, other.startPos)
        r = max(self.endPos, other.startPos)

        return seq[l:r]
        '''