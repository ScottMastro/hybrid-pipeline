#from reverse_complement import reverse_complement

class Path:
    def __init__(self): 
        self.path = []

    def id_on_end(self, id):
        if len(self.path) == 0:
            return None        
        return self.path[0].has_id(id) or self.path[-1].has_id(id)

    def add_fork(self, fork): self.path.append(fork)
    def add_fork_front(self, fork): self.path.insert(0, fork)
    def extend(self, forks):  self.path.extend(forks)
    def add_path(self, other): self.path.extend(other.path)

    def pop(self, index=None):
        if len(self.path) < 1: return None        
        if index is None:      return self.path.pop()        
        return self.path.pop(index)

    def flip_strands(self, lengthData):
        self.path = self.path[::-1]
        for fork in self.path: fork.flip_strands(lengthData)
        
    def has_id(self, tigId):
        for fork in self.path:
            if fork.has_id(tigId): return True
        return False
    
    def last_fork(self):
        if len(self.path) == 0: return None
        else:                   return self.path[-1]
    def first_fork(self):
        if len(self.path) == 0: return None
        else:                   return self.path[0]

    def head(self, index=0):
        for fork in self.path:
            if not fork.is_Nfork(): 
                index = index -1
                if index < 0: return fork
        return None

    def tail(self, index=0):
        for fork in reversed(self.path):
            if not fork.is_Nfork(): 
                index = index -1
                if index < 0: return fork
        return None

    def __getitem__(self, key):
         return self.path[key]
    def __len__(self):
         return len(self.path)
     
    def __repr__(self):
        return str(self.first_fork()) + "...\t..." + str(self.last_fork()) + "\t(" + str(len(self.path)) + ")"
    def __str__(self):
        return "\n".join([f.__repr__() for f in self.path])
    
    
    '''
    def get_interval(self, tigId, lengthData):
        tigSize = lengthData[str(tigId)]
        
        minPos = tigSize
        maxPos = 0
        for fork in self.path:
            pos = fork.get_pos_by_id(tigId)
            if pos is None: continue
            pos = pos if fork.get_strand_by_id(tigId) == 1 else tigSize - pos
            
            minPos = min(minPos, pos)
            maxPos = max(maxPos, pos)

        return(minPos, maxPos)
        
    def normalize(self, tigId, lengthData):
        tigSize = lengthData[str(tigId)]
        if not self.path[-1].has_id(tigId):
            if not self.path[0].has_id(tigId):
                return
            else:
                if self.path[0].get_strand_by_id(tigId) == -1:
                    self.flip_strands(lengthData)
            return
        
        elif not self.path[0].has_id(tigId):
            if self.path[-1].get_strand_by_id(tigId) == -1:
                self.flip_strands(lengthData)
            return
        else:
            firstPos = self.path[0].get_pos_by_id(tigId)
            firstPos = firstPos if self.path[0].get_strand_by_id(tigId) == 1 else tigSize - firstPos
            lastPos = self.path[-1].get_pos_by_id(tigId)
            lastPos = lastPos if self.path[-1].get_strand_by_id(tigId) == 1 else tigSize - lastPos
            if lastPos < firstPos:
                self.flip_strands(lengthData)
                
            return

    

    
    def head(self):
        if len(self.path) == 0:
            return None
        
        id = self[0].before_id()
        strand = self[0].before_strand()
        
        if self[0].before_switch() == 'q':
            h = Fork(id, None, strand, None, None, None)
            h.switch_query()
            return h
        if self[0].before_switch() == 'r':
            h = Fork(None, None, None, id, None, strand)
            h.switch_reference()
            return h
        return None
    

    def tail(self):
        if len(self.path) == 0:
            return None
        
        id = self[-1].after_id()
        strand = self[-1].after_strand()
        
        if self[-1].after_switch() == 'q':
            t = Fork(id, None, strand, None, None, None)
            t.switch_reference()
            return t
        if self[-1].after_switch() == 'r':
            t = Fork(None, None, None, id, None, strand)
            t.switch_query()
            return t
        return None    
        
            def get_fork_sequence(self, seqData, nbases=1000, q=True):
    
        forks = dict()
        for fork in self:
            
            if (q and fork.switch == 'r') or (not q and fork.switch == 'q'):
                id = fork.before_id()
                pos = fork.before_pos()
                seq = str(seqData[str(id)])
                if fork.before_strand() == -1:
                    pos = len(seq) - pos
    
                forkSeq = seq[pos-nbases:pos+nbases]
    
                if fork.before_strand() == -1:
                    forkSeq = reverse_complement(forkSeq)
            
                forks[fork] = forkSeq
                
            if (not q and fork.switch == 'r') or (q and fork.switch == 'q'):
                id = fork.after_id()
                pos = fork.after_pos()
                seq = str(seqData[str(id)])
                if fork.after_strand() == -1:
                    pos = len(seq) - pos
    
                forkSeq = seq[pos-nbases:pos+nbases]
    
                if fork.after_strand() == -1:
                    forkSeq = reverse_complement(forkSeq)
            
                forks[fork] = forkSeq
    
        return forks        
    '''
    



class Fork:
    
    def __init__(self, qid, qpos, qstrand, rid, rpos, rstrand, Nfork=False):
        self.qid, self.rid = qid, rid
        self.qpos, self.rpos = qpos, rpos
        self.qstrand, self.rstrand = qstrand, rstrand
        
        self.switch = "q"
        self.Nfork = Nfork
             
    def is_Nfork(self):            return self.Nfork
    def is_switch_query(self):     return self.switch == "q"
    def is_switch_reference(self): return self.switch ==  "r"

    def before_id(self):     return self.rid     if self.switch == "q" else self.qid
    def after_id(self):      return self.qid     if self.switch == "q" else self.rid
    def before_pos(self):    return self.rpos    if self.switch == "q" else self.qpos
    def after_pos(self):     return self.qpos    if self.switch == "q" else self.rpos
    def before_strand(self): return self.rstrand if self.switch == "q" else self.qstrand
    def after_strand(self):  return self.qstrand if self.switch == "q" else self.rstrand
    def before_switch(self): return "r"          if self.switch == "q" else "q"
    def after_switch(self):  return self.switch

    def has_id(self, tigId): return self.rid == tigId or self.qid == tigId

    def switch_query(self):     self.switch = "q"
    def switch_reference(self): self.switch = "r"
    def flip_switch(self):      self.switch = "r" if self.switch == "q" else "q"
    
    def flip_strands(self, lengthData, makeCopy=False):
        
        if self.Nfork:
            if makeCopy: return get_Nfork()
            else:        return
        
        if makeCopy:
            flipped = Fork(self.qid, self.qpos, self.qstrand, \
                           self.rid, self.rpos, self.rstrand)
            flipped.switch = self.switch
            flipped.flip_strands(lengthData)
            return flipped
            
        self.qpos = lengthData[str(self.qid)] - self.qpos
        self.qstrand = self.qstrand*(-1)
        
        self.rpos = lengthData[str(self.rid)] - self.rpos
        self.rstrand = self.rstrand*(-1)
        self.flip_switch()


    def get_pos(self, q=True):
        if q: return self.qpos
        else: return self.rpos
        
    def get_pos_by_id(self, tigId):
        if self.qid == tigId:    return self.qpos
        elif self.rid == tigId:  return self.rpos
        else:                    return None
                
    def get_strand_by_id(self, tigId):
        if self.qid == tigId:    return self.qstrand
        elif self.rid == tigId:  return self.rstrand
        else:                    return None
        
    def __repr__(self):
        sq = str(self.qid) + " " + str(self.qpos) + " (" + str(self.qstrand) + ")"
        sr = str(self.rid) + " " + str(self.rpos) + " (" + str(self.rstrand) + ")"
        
        if self.switch == 'q': return sr + " " + sq
        else:                  return sq + " " + sr
        
    def __str__(self):
        return self.__repr__()
    
    
def get_Nfork(): return Fork("NNN", -1, 0, "NNN", -1, 0, Nfork=True)

def get_Npath():
    p = Path()
    p.add_fork(get_Nfork())
    return p
    