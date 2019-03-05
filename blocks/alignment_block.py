class Block:
    dir=0
    chunks=[]
    rid=""
    qid=""
  
    def __init__(self, chunk, rid, qid):
        self.dir=0
        self.chunks=[]
        self.rid=rid
        self.qid=qid

        self.why=""
        self.chunks.append(chunk)

    def get_dir(self, string=False): 
        if self.dir == 0:
            return "*" if string else 0
        if self.dir == 1:
            return "+" if string else 1
        if self.dir == -1:
            return "-" if string else -1
        
    def verify_direction(self, chunk): 
        d=1
        if(chunk.left_r() > self.chunks[-1].left_r()):
            d=-1
        
        if self.dir == 0:
            self.dir = d
            return True
        
        return d == self.dir
    
    def add(self, chunk): self.chunks.append(chunk)

    #left and right relative to query
    def left_chunk(self): 
        return self.chunks[0] if self.dir == 1 else self.chunks[-1]
    def right_chunk(self): 
        return self.chunks[0] if self.dir == -1 else self.chunks[-1]
    
    def trim_left(self, n):
        if n >= len(self.chunks):
            self.chunks = []
        elif self.dir == 1:
            self.chunks = self.chunks[n:]
        else:
            self.chunks = self.chunks[:-n]
    
    def trim_right(self, n):
        if n >= len(self.chunks):
            self.chunks = []
        elif self.dir == 1:
            self.chunks = self.chunks[:-n]
        else:
            self.chunks = self.chunks[n:]   
    
    def left_q(self): return self.chunks[0].left_q() if self.dir == 1 else self.chunks[-1].left_q()
    def right_q(self): return self.chunks[0].right_q() if self.dir == -1 else self.chunks[-1].right_q()
    def left_r(self): return self.chunks[0].left_r() if self.dir == 1 else self.chunks[-1].left_r()
    def right_r(self): return self.chunks[0].right_r() if self.dir == -1 else self.chunks[-1].right_r()

    def size_q(self): return sum([chunk.size_q() for chunk in self.chunks])
    def size_r(self): return sum([chunk.size_r() for chunk in self.chunks])

    def nchunk(self): return len(self.chunks)
    def nchunk_from(self, start):
        n=0
        for chunk in self.chunks:
            if chunk.id > start: n = n+1
        return n         
            
    def __repr__(self):
        if len(self.chunks) == 0:
            return "Empty"
        return str(self.left_chunk().id) + "-" + str(self.right_chunk().id) + " (" + str(self.nchunk()) + ")"

    def __str__(self):
        if len(self.chunks) == 0:
            return "Empty"
        return "qid=" + str(self.qid) + " rid=" + str(self.rid) + \
                " chunks=" + str(self.left_chunk().id) + "-" + \
                str(self.right_chunk().id) + " (" + str(self.nchunk()) + ")" + \
                " dir=" + self.get_dir(string=True)

    def __lt__(self, other):
        return self.left_chunk().id < other.left_chunk().id

def chunk_split(lchunk, rchunk, repeatOverlap):
    overlapSize = rchunk.right_q() - lchunk.left_q()
    overlapStart = lchunk.right_q()
    overlapEnd = rchunk.left_q()
    return (overlapEnd - overlapStart) / overlapSize < repeatOverlap

def remove_redundancy(blocks, chunkPerBlock, repeatOverlap):
    blocks.sort()

    i=0
    while ( i < len(blocks)-1 ):
        
        x = blocks[i+1].left_chunk().id - blocks[i].right_chunk().id
        
        #no overlap
        if x > 0:
            i = i + 1 
            continue
        
        #overlapping
        if x <= 0:
            #one repeated chunk
            if x ==0:
                if(chunk_split(blocks[i].right_chunk(), blocks[i+1].left_chunk(), repeatOverlap)):
                    i = i + 1 
                    continue
            
            #left block is bigger
            if(blocks[i].nchunk() > blocks[i+1].nchunk()):
                blocks[i+1].trim_left(-x + 1)
                if(blocks[i+1].nchunk() < chunkPerBlock):
                    blocks.pop(i+1)
            #right block is bigger
            else:
                blocks[i].trim_right(-x + 1)
                if(blocks[i].nchunk() < chunkPerBlock):
                    blocks.pop(i)
            
    return blocks

def construct_blocks(rid, qid, chunks, chunkSize=1000, chunkPerBlock=3, maxDist=5000, chunkSkip=2, bpSkip=2000):
    '''
    Takes a row from the aligndf file, creates chunks and connects them into blocks.
    Valid blocks are returned, leftover chunks are returned in trash.
    '''
    blocks = []
    trash = []

    while len(chunks) > 0:
         
        i = 0
        chunk = chunks.pop()
        block = Block(chunk, rid, qid)
                
        while i < len(chunks):

            i=i+1
            nextChunk = chunks[-i]
            idDist = chunk.id - nextChunk.id
            idSkip = idDist-1
            allowableGap = max(maxDist, idSkip*bpSkip)
            bpDist = abs(chunk.rstart - nextChunk.rstart) - chunkSize

            #chunk id improper
            if idDist < 1: continue         
            if idSkip > chunkSkip: break
            
            #bp gap too big
            if bpDist > allowableGap:
                continue
            
            #chunk oriented incorrectly
            if not block.verify_direction(nextChunk):
                continue
            
            #all good, add chunk and remove from list
            block.add(nextChunk)
            chunk = chunks.pop(-i)
            i=i-1

            if len(chunks) == 0:
                break
            
            
        if(block.nchunk() > chunkPerBlock):
            blocks.append(block)
            #block.print_info()
        else:
            trash.append(block)
        
    return (blocks, trash)