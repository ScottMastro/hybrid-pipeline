class Interval:
        
    def __init__(self, qid="?", rid="?"):
        self.rid=rid
        self.qid=qid
        self.components = []
        self.sorted=True
        self.fill=[0,0,0]

    def add(self, interval):
        self.components.append(interval)  
        self.sorted=False
        
    def ncomponents(self): return len(self.components)
        
    def span(self, q=True): 
        return abs(self.start(q) - self.end(q))    
    def sum(self, q=True): 
        return sum([component.span(q) for component in self.components])

    def get_dir(self, q=True, string=False):
        if(self.start(q) < self.end(q)):
            return "+" if string else 1
        elif(self.end(q) < self.start(q)):
            return "-" if string else -1
        else:
            return "?" if string else 0
    
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
        return self.components[-1].end(q)

    def left_of(self, other, q=True):
        return self.left(q) < other.left(q)
    def right_of(self, other, q=True):
        return self.right(q) < other.right(q)
    def overlapping(self, other, q=True):
        return self.right(q) >= other.left(q) and self.left(q) <= other.right(q)
    
    '''
    def __eq__(self, other):
        if not self.sorted: self.sort()
        if not other.sorted: other.sort()

        if len(self.components) != len(other.components):
            return False
        if len(self.rid) != len(other.rid):
            return False
        if len(self.qid) != len(other.qid):
            return False

        for c1, c2 in zip(self.components, other.components):
            if not c1.__eq__(c2):
                return False
            
        return True
    '''
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

    def __init__(self, cid, rst, red, qst, qed, qid="?", rid="?"):
        
        super().__init__(qid, rid)
        self.id=cid
        self.rstart=rst
        self.rend=red
        self.qstart=qst
        self.qend=qed
    
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
    
    '''
    def __eq__(self, other):
        if len(self.cid) != len(self.cid): return False
        if len(self.rid) != len(other.rid): return False
        if len(self.qid) != len(other.qid): return False
        if len(self.qstart) != len(other.qstart): return False
        if len(self.qend) != len(other.qend): return False
        if len(self.rstart) != len(other.rstart): return False
        if len(self.rend) != len(other.rend): return False
    '''
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
        self.rdir=chunk.get_dir(q="False")
        
    def verify_direction(self, chunk): 
        return self.rdir == chunk.get_dir(q="False")
    
    def trim_left(self, n):
        if not self.sorted: self.sort()
        
        if n >= len(self.components):
            self.components = []
        elif self.rdir == 1:
            self.components = self.components[n:]
        else:
            self.components = self.components[:-n]
    
    def trim_right(self, n):
        if not self.sorted: self.sort()

        if n >= len(self.components):
            self.components = []
        elif self.rdir == 1:
            self.components = self.components[:-n]
        else:
            self.components = self.components[n:]   
    
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

    '''
    def left_block(self): 
        if(self.blocks[0].left_q() < self.blocks[-1].left_q()):
            return self.blocks[0]
        else:
            return self.blocks[-1]
    def right_block(self): 
        if(self.blocks[0].left_q() < self.blocks[-1].left_q()):
            return self.blocks[-1]
        else:
            return self.blocks[0]
    '''
    
    '''
    def cov_size(self, q=True): return sum([block.size_q() for block in self.blocks])
    def cov_size_r(self): return sum([block.size_r() for block in self.blocks])    
    def size_q(self): return abs(self.left_q() - self.right_q())
    def size_r(self): return abs(self.left_r() - self.right_r())
    '''    
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
        return "qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                " query=" + str(self.left_q()) + "-" + str(self.right_q()) + \
                " reference=" + str(self.left_r()) + "-" + str(self.right_r()) + \
                " nblocks=" + str(self.nblocks())
