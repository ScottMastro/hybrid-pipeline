LEFT, RIGHT = 'l', 'r'

class Interval:
    '''
    A region of homology between a query and reference contig.
    Each Interval has a qid and rid and a set of component Intervals
    which define smaller and higher resolution regions of homology.
    
    ex. Interval may be a single block, components would be Chunks
    
    The most basic level of the Interval class is the Chunk class,
    which represent exact alignments between query and reference contigs.
    
    Most functions have parameter q where q=True means to evaluate the function
    with respect to the query data and q=False is with respect to the reference.
    '''

    def __init__(self, qid="?", rid="?", subclass=None):
        self.iid=None
        self.rid=str(rid)
        self.qid=str(qid)
        self.components = []
        self.subclass = subclass
    
    def set_interval_id(self, iid):
        self.iid=str(iid)

    def add(self, interval): self.components.append(interval)  
    def sort(self, q=True): self.components.sort(key=lambda c: c.left(q))

    def left(self, q=True): 
        '''Gets the smallest (leftmost) position in this Interval.'''
        return min([x.left(q) for x in self.components])
    def right(self, q=True):
        '''Gets the largest (rightmost) position in this Interval.'''
        return max([x.right(q) for x in self.components])
    def start(self, q=True): 
        '''Gets the start position of the first component in this Interval.'''
        return self.components[0].start(q)
    def end(self, q=True): 
        '''Gets the end position of the last component in this Interval.'''
        return self.components[-1].end(q)

    def similarity(self): 
        '''(length * percent identity) across entire Interval '''
        return sum([component.similarity() for component in self.components])
    def similarity_between(self, pos1, pos2, q=True): 
        '''(length * percent identity) across a given range of positions'''
        return sum([component.similarity_between(pos1, pos2, q) for component in self.components])

    def should_trim(self, position, side, q=True):
        '''
        Determines if this Interval should be trimmed.
        position determines where the Interval will be trimed
        side (LEFT or RIGHT) determines which side of position is to trimmed.
        Returns True if entire Interval is inside the region to be trimmed,
        False otherwise
        '''
        if side == LEFT and self.right(q) < position: return True
        if side == RIGHT and self.left(q) > position: return True
        return False

    def trim(self, position, side, q=True):
        trimmed, i = [], 0
        while i < len(self):
            if self[i].should_trim(position, side, q):
                x = self.components.pop(i)
                trimmed.append(x)
            else: 
                self[i].trim(position, side, q)
                if len(self[i]) < 1: self.components.pop(i)
                else: i = i+1
        return trimmed

    def trim_left(self, position, q=True):
        '''
        Components to the left of position will be trimmed.
        Returns components that have been fully trimmed.
        '''
        return self.trim(position, LEFT, q)
                    
    def trim_right(self, position, q=True):
        '''
        Components to the right of position will be trimmed.
        Returns components that have been fully trimmed.
        '''
        return self.trim(position, RIGHT, q)

    def span(self, q=True): 
        '''Distance between the first position and last position'''
        if len(self) < 1: return 0
        return abs(self.start(q) - self.end(q))    
    def span_sum(self, q=True):
        '''Span of all component Intervals, summed together'''
        if len(self) < 1: return 0
        return sum([component.span(q) for component in self.components])
    
    def total_coverage(self):
        '''coverage across all components'''
        return sum([x.total_coverage() for x in self.components])
    def coverages(self):
        '''coverage across all components'''
        return [x.total_coverage() for x in self.components]

    def percent_identity(self):
        '''Average percent identity of all components'''
        covs = self.coverages()
        pcid = sum([x.percent_identity()*cov for x,cov in zip(self.components, covs)])
        return pcid / sum(covs)
    
    def percent_identities(self): 
        '''Gets the percent identity of all components, returns as list'''
        return [x.percent_identity() for x in self.components]

    def is_consistent(self, q=True):
        '''Returns True if the sum of the distance between consecutive 
        components is smaller than the sum of the distance between consecutive
        reverse complement components. Returns False otherwise'''

        dist=0
        revDist=0
        for i in range(len(self.components)-1):
            dist = dist + abs(self.components[i].end(q) - self.components[i+1].start(q))
            revDist = revDist + abs(self.components[i].start(q) - self.components[i+1].end(q))

        if revDist < dist:
            return False
        
        return True
    
    
    
    def get_dir(self, q=True, string=False):
        if(self.start(q) < self.end(q)):    return "+" if string else 1
        elif(self.end(q) < self.start(q)):  return "-" if string else -1
        else:                               return "?" if string else 0

    def left_of(self, other, q=True):
        return self.left(q) < other.left(q)
    def right_of(self, other, q=True):
        return self.right(q) < other.right(q)
    def overlapping(self, other, q=True):
        return self.right(q) >= other.left(q) and self.left(q) <= other.right(q)
    
    def contains_position(self, basePosition, q=True): 
        return basePosition >= self.left(q) and basePosition <= self.right(q)

    def closest_corresponding_position(self, basePosition, q=True, side=None): 
        
        if self.contains_position(basePosition, q):
            
            distances = [min(abs(basePosition-component.start(q)), abs(basePosition-component.end(q))) \
                         for component in self.components]
            
            index_min = min(range(len(distances)), key=distances.__getitem__)
            
            if side is None or self.components[index_min].contains_position(basePosition, q) :
                return self.components[index_min].closest_corresponding_position(basePosition, q, side)
            #before index_min
            elif abs(basePosition-self.components[index_min].start(q)) < \
                abs(basePosition-self.components[index_min].end(q)):
                    if side == RIGHT:
                        return self.components[index_min].start(not q)
                    elif side == LEFT:
                        #min_index-1 should exist because self contains the point
                        return self.components[index_min-1].end(not q)
            #after index_min
            else:
                    if side == RIGHT:
                        #min_index+1 should exist because self contains the point
                        return self.components[index_min+1].start(not q)
                    elif side == LEFT:
                        return self.components[index_min].end(not q)
                    
        else:
            if abs(basePosition-self.start(q)) < abs(basePosition-self.end(q)):
                return self.start(not q)
            else:
                return self.end(not q)       
        

    def __eq__(self, other): 
        if not isinstance(other, Interval):
            return NotImplemented
        
        if not other.rid == self.rid: return False
        if not other.qid == self.qid: return False
        if not len(other.components) == len(self.components): return False

        for i, cpnt in enumerate(self.components):
            if not cpnt == other.components[i]:return False

        return True

    def __len__(self):
        return len(self.components)
            
    def __getitem__(self, key):
         return self.components[key]

    def __repr__(self):
        string = "" if self.subclass is None else str(self.subclass + ":\t")
        string = string + "qid=" + str(self.qid) + "\trid=" + str(self.rid) + "\t"
        string = string + "len=" + str(len(self))
        return string

    def __str__(self):
        return self.__repr__() + "\n".join([c.__str__() for c in self.components])

def construct_interval(component):
    newInterval = Interval(component.qid, component.rid)
    newInterval.add(component)
    return newInterval


class Chunk(Interval):

    def __init__(self, cid, rst, red, qst, qed, pcid, alen, qid="?", rid="?"):
        
        super().__init__(qid, rid)
        self.id=cid
        self.rstart,self.rend = rst,red
        self.qstart,self.qend = qst,qed
        self.pcid=pcid
        self.alen=alen
            
    def start(self, q=True): return self.qstart if q else self.rstart
    def end(self, q=True): return self.qend if q else self.rend
    
    def left(self, q=True): 
        return min(self.qstart, self.qend) if q else min(self.rstart, self.rend)
    def right(self, q=True): 
        return max(self.qstart, self.qend) if q else max(self.rstart, self.rend)
    def span(self, q=True): 
        return abs(self.qstart - self.qend) if q else abs(self.rstart - self.rend)
    def mid(self, q=True):
        return int((self.qstart + self.qend)/2) if q else int((self.rstart + self.rend)/2)
    
    def similarity(self): return (self.alen*self.pcid)/100.0
    def similarity_between(self, pos1, pos2, q=True):
        minPos, maxPos = min(pos1, pos2), max(pos1, pos2)
           
        o1 = max(0, minPos - self.left(q) )
        o2 = max(0, self.right(q) - maxPos)
        o = o1 + o2
        
        if o > self.span(q): return 0
        else: return (self.span(q) - o)*(self.alen*self.pcid)/100.0

    def trim(self, position, side, q=True): return
    def should_trim(self, position, side, q=True):
        '''
        Determines if this Interval should be trimmed.
        position determines where the Interval will be trimed
        side (LEFT or RIGHT) determines which side of position is to trimmed.
        Returns True if Chunk overlaps the region to be trimmed, False otherwise.
        
        Note: this implementation is more strict than the Interval method 
        as it doesn't require complete overlaps.
        '''
        if side == LEFT and self.left(q) < position: return True
        if side == RIGHT and self.right(q) > position: return True
        return False
    
    def percent_identity(self): return self.pcid
    def total_coverage(self): return self.alen

    def closest_corresponding_position(self, basePosition, q=True, side=None): 
        displacement = self.start(q) - basePosition
        direction = self.get_dir(not q) * self.get_dir(q)
        return self.start(not q) - displacement*direction

    def __len__(self): return 1

    def __eq__(self, other): 
        if not isinstance(other, Chunk): return False
        if not other.qid == self.qid: return False
        if not other.rid == self.rid: return False
        if not other.rstart == self.rstart: return False
        if not other.rend == self.rend: return False
        if not other.qstart == self.qstart: return False
        if not other.qend == self.qend: return False

        return True

    def __repr__(self):
        return str(self.id) +  " (" + str(self.qstart) + "-" + str(self.qend) + ")"

    def __str__(self):
        return "CHUNK " + str(self.id) + " q=" + str(self.qid) + \
    " (" + str(self.qstart) + "-" + str(self.qend) + ")" + \
    " r=" + str(self.rid) + " (" + str(self.rstart) + "-" + str(self.rend) + ")"

