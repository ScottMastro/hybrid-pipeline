class Segment:
    id=""
    left=None
    right=None
    length=-1
    dir=1

    def __init__(self, tigId, left, right, dir, length):
        self.id=tigId
        self.left=left
        self.right=right
        self.dir=dir
        self.length=length

    def to_string(self):
        return self.id + " " + str(self.left) + " " + str(self.right) + " " + str(self.dir) 
    
    def size(self):
        return self.right-self.left
        
    def equals(self, otherSegment):
        
        if self.left != otherSegment.left:
            return False
        if self.right != otherSegment.right:
            return False
        if self.dir != otherSegment.dir:
            return False
        if self.id != otherSegment.id:
            return False
        return True
    
    def reverse_complement(self, getCopy=False):
        
        #todo: might be off-by-one error here????
        leftrc = self.length - self.left
        rightrc = self.length - self.right
        dirrc = self.dir*(-1)
                
        if getCopy:
            return Segment(self.id, rightrc, leftrc, dirrc, self.length)
        else:
            self.left = rightrc
            self.right = leftrc
            self.dir = dirrc
            return self

class Path:
    
    segments=[]
    ids=set()
    
    def add_segments(self, segs):
        for seg in segs:
            self.segments.append(seg)
            self.ids.add(seg.id)
        
    def collapse_segments(self):
        collapsed=[]
        prevSeg = self.segments[0]
        for i in range(1, len(self.segments)):
            seg = self.segments[i]
            
            if seg.id == prevSeg.id:       #same contig
                if seg.dir != prevSeg.dir:   #different strands
                    print("WARNING: unexpected segment join (reverse complements)")
                else:
                    prevSeg = Segment(prevSeg.id, prevSeg.left, seg.right, prevSeg.dir, prevSeg.length)
                    continue
            
            collapsed.append(prevSeg)
            prevSeg = seg
            
        collapsed.append(prevSeg)
        self.segments = collapsed
        
    def __init__(self, segs=None):
        self.segments=[]
        self.ids=set()
        if segs is not None:
            self.add_segments(segs)
    
    def get_left_id(self):
        if len(self.segments) < 1: return None
        return self.segments[0].id

    def get_right_id(self):
        if len(self.segments) < 1: return None
        return self.segments[-1].id
        
    def print_path(self):
        string = ""
        for seg in self.segments:
            string = string + "(" + seg.to_string() + ") "
        print(string)
        
    def is_empty(self):
        return len(self.segments) == 0
    
    def equals(self, otherPath):
        
        if len(self.segments) != len(otherPath.segments):
            return False
        if self.is_empty() and otherPath.is_empty():
            return True
        for i in range(len(self.segments)):
            if not self.segments[i].equals(otherPath.segments[i]):
                return False
            
        return True
        
    def reverse_complement(self):
        
        reverseSegs = []
        
        nsegs = len(self.segments)
        for i in range(nsegs):
            reverseSeg = self.segments[nsegs-(i+1)].reverse_complement(getCopy=True)
            reverseSegs.append(reverseSeg)  
        
        return Path(reverseSegs)
    

def join_paths(left, right):
    if left.is_empty():
        return right
    
    if right.is_empty():
        return left
    
    segments=left.segments
    segments.extend(right.segments)
    joined = Path(segments)
    joined.collapse_segments()
    return joined

def can_join(left, right):
    if left.is_empty() or right.is_empty():
        return True

    leftSeg = left.segments[-1]
    rightSeg = right.segments[0]
    
    if not leftSeg.id == rightSeg.id:
        return False
    
    if not leftSeg.dir == rightSeg.dir:
        return False
    
    if not rightSeg.right >= leftSeg.left:
        return False
    
    if not rightSeg.right <= leftSeg.right:
        return False

    return True
    
def try_join(seg1, seg2, getResult=True):
    if seg1.is_empty() or seg2.is_empty():
        return True if not getResult else join_paths(seg1, seg2)
    
    if can_join(seg1, seg2):
        return True if not getResult else join_paths(seg1, seg2)
    
    if can_join(seg2, seg1):
        return True if not getResult else join_paths(seg2, seg1)
    
    seg1rc = seg1.reverse_complement()
    
    if can_join(seg1rc, seg2):
        return True if not getResult else join_paths(seg1rc, seg2)
    
    if can_join(seg2, seg1rc):
        return True if not getResult else join_paths(seg2, seg1rc)

    return False if not getResult else None

def pp(paths):
    
    try:
        iterator = iter(paths)
    except TypeError:
        paths.print_path()  
    else:
        for path in paths:
            path.print_path()
