class Path:
    '''
    A ordered set of Forks used to specify a traversal though both
    query and reference sequence that defines a hybrid sequence.
        
    Most functions have parameter q where q=True means to evaluate the function
    with respect to the query data and q=False is with respect to the reference.
    '''
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

    def extract_subpath(self, idx1=None, idx2=None): 
        subpath = Path()
        subpath.extend(self[idx1:idx2])
        return subpath
    
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
    
    def path_length(self, tigId=None, strict=True):
        '''
        Calculates the total length of the path. If tigId is provided,
        only counts segments corresponding to tigId.
        If strict is False, any segments containing tigId will be counted otherwise
        segments will only be counted if after_id and before_id match tigId
        '''
    
        if len(self.path) < 2: return 0
        pathSum = 0
    
        prevFork = self.path[0]
        for fork in self.path[1:]:
            if not fork.is_Nfork() and not prevFork.is_Nfork():
                                
                if (tigId is None) or \
                    (strict and prevFork.after_id() == tigId and fork.before_id() == tigId) or \
                    (not strict and prevFork.has_id(tigId) and fork.has_id(tigId)):
                        pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
                    
            prevFork = fork
            
        return pathSum

    def __eq__(self, other): 
        if not isinstance(other, Path): return False
        if not len(other.path) == len(self.path): return False

        for i, fork in enumerate(self.path):
            if not fork == other[i]: return False

        return True

    def __getitem__(self, key):
         return self.path[key]
    def __len__(self):
         return len(self.path)
     
    def __repr__(self):
        return str(self.first_fork()) + "...\t..." + str(self.last_fork()) + "\t(" + str(len(self.path)) + ")"
    def __str__(self):
        return "\n".join([f.__repr__() for f in self.path])
    

class Fork:
    '''
    A dual position within both query and reference sequence that specifies
    a transition from one sequence to the other; defines a hybrid sequence.
    
    Contig id, position, and strand are explicitly defined for both query and ref.
    Most functions have parameter q where q=True means to evaluate the function
    with respect to the query data and q=False is with respect to the reference.
    
    A special case of this class ("Nfork") exists to define an unknown sequence
    which will be represented by sequential NNN characters.
    '''
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
        return self.qpos if q else self.rpos
    def get_id(self, q=True):
        return self.qid if q else self.rid
    def get_strand(self, q=True):
        return self.qstrand if q else self.rstrand

    def get_pos_by_id(self, tigId):
        if self.qid == tigId:    return self.qpos
        elif self.rid == tigId:  return self.rpos
        else:                    return None

    def get_strand_by_id(self, tigId):
        if self.qid == tigId:    return self.qstrand
        elif self.rid == tigId:  return self.rstrand
        else:                    return None
    
    def normalize_pos(self, pos, strand, tigId, lengthData):
        if pos is None or self.is_Nfork(): return pos
        if strand == -1:
            return lengthData[tigId] - pos
        return pos
    
    def before_pos_norm(self, lengthData):   
        return self.normalize_pos(self.before_pos(), self.before_strand(), 
                                  self.before_id(), lengthData)
    def after_pos_norm(self, lengthData):   
        return self.normalize_pos(self.after_pos(), self.after_strand(), 
                                  self.after_id(), lengthData)
    def get_pos_norm(self, lengthData, q=True):
        if q: return self.normalize_pos(self.qpos, self.qstrand, 
                                        self.qid, lengthData)
        return self.normalize_pos(self.rpos, self.rstrand, 
                                        self.rid, lengthData)
    def get_pos_by_id_norm(self, tigId, lengthData):
        if self.qid == tigId:    return self.get_pos_norm(lengthData, q=True)
        elif self.rid == tigId:  return self.get_pos_norm(lengthData, q=False)
        else:                    return None

    def __eq__(self, other): 
        if not isinstance(other, Fork): return False        
        if not other.qid == self.qid: return False
        if not other.rid == self.rid: return False
        if not other.qpos == self.qpos: return False
        if not other.rpos == self.rpos: return False
        if not other.qstrand == self.qstrand: return False
        if not other.rstrand == self.rstrand: return False
        if not other.switch == self.switch: return False

        return True

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
    