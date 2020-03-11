from . import fork as FORK

class Path:
    '''
    A ordered set of Forks used to specify a traversal though both
    query and reference sequence that defines a hybrid sequence.
        
    Most functions have parameter q where q=True means to evaluate the function
    with respect to the query data and q=False is with respect to the reference.
    '''
    def __init__(self): 
        self.path = []
        self.pid = None
        
    def set_path_id(self, pid): 
        self.pid = str(pid)

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
        

def get_Npath():
    p = Path()
    p.add_fork(FORK.get_Nfork())
    return p
    
