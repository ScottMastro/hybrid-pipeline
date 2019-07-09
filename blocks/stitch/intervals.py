class Interval:
        
    def __init__(self, qid="?", rid="?"):
        self.rid=str(rid)
        self.qid=str(qid)
        self.components = []
        self.sorted=True


    '''Get position x percent identity for a range of positions'''
    def coverage(self): 
        return sum([component.coverage() for component in self.components])
    def coverage_between(self, pos1, pos2, q=True): 
        return sum([component.coverage_between(pos1, pos2, q) for component in self.components])

    def add(self, interval):
        self.components.append(interval)  
        self.sorted=False
        
    def ncomponents(self): return len(self.components)
        
    
    def span(self, q=True): 
        if len(self) < 1: return 0
        return abs(self.start(q) - self.end(q))    
    def span_sum(self, q=True):
        if len(self) < 1: return 0
        return sum([component.span(q) for component in self.components])
    
    def percent_identity(self):
        return sum([x.percent_identity() for x in self.components])/len(self.components)
    def percent_identities(self): 
        return [x.percent_identity() for x in self.components]

    def should_trim(self, position, side, q=True):
        if side == 'l':
            if self.right(q) < position:
                return True
        elif side == 'r':
            if self.left(q) > position:
                return True
        return False

    def trim_left(self, position, q=True):
        trimmed = []
        i=0
        while i < len(self):
            if self[i].should_trim(position, 'l', q):
                x = self.components.pop(i)
                trimmed.append(x)
            else: 
                self[i].trim_left(position, q)
                if len(self[i]) < 1:
                    self.components.pop(i)
                else:
                    i = i+1
        return trimmed
                    
    def trim_right(self, position, q=True):
        trimmed = []
        i=0
        while i < len(self):
            if self[i].should_trim(position, 'r', q):
                x = self.components.pop(i)
                trimmed.append(x)
            else: 
                self[i].trim_right(position, q)
                if len(self[i]) < 1:
                    self.components.pop(i)
                else:
                    i = i+1
        return trimmed

    def get_dir(self, q=True, string=False):
        if(self.start(q) < self.end(q)):
            return "+" if string else 1
        elif(self.end(q) < self.start(q)):
            return "-" if string else -1
        else:
            return "?" if string else 0
        
    def is_consistent(self, q=True):
        
        dist=0
        revDist=0
        for i in range(len(self.components)-1):
            dist = dist + abs(self.components[i].end(q) - self.components[i+1].start(q))
            revDist = revDist + abs(self.components[i].start(q) - self.components[i+1].end(q))

        if revDist < dist:
            return False
        
        return True

        
    def sort(self):
        self.components.sort()
        self.sorted=True
    
    def left_id(self, q=True):
        if not self.sorted: self.sort()
        return self.components[0].left_id()
    def right_id(self, q=True): 
        if not self.sorted: self.sort()
        return self.components[-1].right_id()
    
    def left(self, q=True):
        if not self.sorted: self.sort()
        return min(self.components[0].left(q), self.components[-1].left(q))
    def right(self, q=True): 
        if not self.sorted: self.sort()
        return max(self.components[0].right(q), self.components[-1].right(q))

    def start(self, q=True): 
        if not self.sorted: self.sort()
        return self.components[0].start(q)
    def end(self, q=True): 
        if not self.sorted: self.sort()
        return self.components[-1].end(q)

    def left_of(self, other, q=True):
        return self.left(q) < other.left(q)
    def right_of(self, other, q=True):
        return self.right(q) < other.right(q)
    def overlapping(self, other, q=True):
        return self.right(q) >= other.left(q) and self.left(q) <= other.right(q)
    
    def contains_position(self, basePosition, q=True): 
        return basePosition >= self.left(q) and basePosition <= self.right(q)

    def closest_corresponding_position(self, basePosition, q=True, side=None): 
        if not self.sorted: self.sort()
        
        if self.contains_position(basePosition, q):
            
            distances = [min(abs(basePosition-component.start(q)), abs(basePosition-component.end(q))) \
                         for component in self.components]
            
            index_min = min(range(len(distances)), key=distances.__getitem__)
            
            if side is None or self.components[index_min].contains_position(basePosition, q) :
                return self.components[index_min].closest_corresponding_position(basePosition, q, side)
            #before index_min
            elif abs(basePosition-self.components[index_min].start(q)) < \
                abs(basePosition-self.components[index_min].end(q)):
                    if side == 'r':
                        return self.components[index_min].start(not q)
                    elif side == 'l':
                        #min_index-1 should exist because self contains the point
                        return self.components[index_min-1].end(not q)
            #after index_min
            else:
                    if side == 'r':
                        #min_index+1 should exist because self contains the point
                        return self.components[index_min+1].start(not q)
                    elif side == 'l':
                        return self.components[index_min].end(not q)
                    
        else:
            if abs(basePosition-self.start(q)) < abs(basePosition-self.end(q)):
                return self.start(not q)
            else:
                return self.end(not q)       
        

    def __len__(self):
        return len(self.components)
    
    #for sorting, sort by leftmost id
    def __lt__(self, other):
        value = self.left_id() - other.left_id()
        if value < 0: return True
        elif value > 0: return False
        else: return self.right_id() < other.right_id()
        
    def __getitem__(self, key):
         return self.components[key]


class Chunk(Interval):

    def __init__(self, cid, rst, red, qst, qed, pcid, alen, qid="?", rid="?"):
        
        super().__init__(qid, rid)
        self.id=cid
        self.rstart=rst
        self.rend=red
        self.qstart=qst
        self.qend=qed
        self.pcid=pcid
        self.alen=alen
        
    def get_dir(self, q=True, string=False):
        #query is by default + direction
        if q:
            return "+" if string else 1

        if(self.rstart < self.rend):
            return "+" if string else 1
        elif(self.rend < self.rstart):
            return "-" if string else -1
        else:
            return "?" if string else 0
    
    def left_id(self): return self.id 
    def right_id(self): return self.id 
    
    def start(self, q=True): 
        if q: return self.qstart
        return self.rstart
    def end(self, q=True):         
        if q: return self.qend
        return self.rend
    
    def left(self, q=True): 
        if q: return min(self.qstart, self.qend)
        return min(self.rstart, self.rend)
    def right(self, q=True): 
        if q: return max(self.qstart, self.qend)
        return max(self.rstart, self.rend)
    
    def span(self, q=True): 
        if q: return abs(self.qstart - self.qend)
        return abs(self.rstart - self.rend)
    
    def coverage(self): 
        return (self.alen*self.pcid)/100.0
    def coverage_between(self, pos1, pos2, q=True):
        if pos1 < pos2:
            minPos, maxPos = pos1, pos2
        else:
            minPos, maxPos = pos2, pos1
            
        o1 = max(0, minPos - self.left(q) )
        o2 = max(0, self.right(q) - maxPos)
        o = o1+o2
        
        if o > self.span(q):
            return 0
        
        return (self.span(q) - o)*(self.alen*self.pcid)/100.0

    def trim_left(self, position, q=True):
        return
    def trim_right(self, position, q=True):
        return
        
    def should_trim(self, position, side, q=True):
        if side == 'l':
            if self.left(q) < position:
                return True
        elif side == 'r':
            if self.right(q) > position:
                return True
        return False
    
    def percent_identity(self): 
        return self.pcid
    
    def closest_corresponding_position(self, basePosition, q=True, side=None): 
        displacement = self.start(q) - basePosition
        
        direction = self.get_dir(not q) * self.get_dir(q)
        
        return self.start(not q) - displacement*direction

    def __len__(self):
        return 1

    def __repr__(self):
        return str(self.id) +  " (" + str(self.qstart) + "-" + str(self.qend) + ")"

    def __str__(self):
        return "CHUNK " + str(self.id) + " q=" + str(self.qid) + \
    " (" + str(self.qstart) + "-" + str(self.qend) + ")" + \
    " r=" + str(self.rid) + " (" + str(self.rstart) + "-" + str(self.rend) + ")"


class Block(Interval):
  
    def __init__(self, chunk):
        super().__init__(chunk.qid, chunk.rid)        
        self.rdir=0
        self.add(chunk)
        self.rdir=chunk.get_dir(q=False)
    
    def nchunk(self): return len(self)
    def nchunk_from(self, start):
        n=0
        for chunk in self.components:
            if chunk.id > start: n = n+1
        return n         
                
    def __repr__(self):
        if len(self.components) == 0: return "Empty"
        return str(self.left_id()) + "-" + str(self.right_id())

    def __str__(self):
        if len(self.components) == 0:
            return "Empty BLOCK"
        return "BLOCK qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                "\nchunks=\n" + "\n".join([chunk.__str__() for chunk in self.components])
                
class Megablock(Interval):
  
    def __init__(self, blocks):
        if len(blocks) > 0:
            super().__init__(blocks[0].qid, blocks[0].rid)        
        else:
            super().__init__()        
        self.components = blocks


    def __repr__(self):
        return str(self.components) 

    def __str__(self):
        return "MEGABLOCK qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                " query=" + str(self.start(q=True)) + "-" + str(self.end(q=True)) + \
                " reference=" + str(self.start(q=False)) + "-" + str(self.end(q=False)) + \
                " nblocks=" + str(len(self.components))
                
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
        return "id=" + str(self.qid)