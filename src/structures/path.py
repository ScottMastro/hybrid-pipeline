#==================================================
# Fork class
#==================================================

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

#==================================================
# Path class
#==================================================

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

    def get_tig_ids(self, q=None):
       #True for query, False for ref, None for both                
        ids = set()
        for fork in self.path:
            if q is None or not q:
                ids.add(fork.rid)
            if q is None or q:
                ids.add(fork.qid)
    
        return ids        

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
    p.add_fork(get_Nfork())
    return p

#==================================================
# Helper functions
#==================================================

def consistent_forks(beforeFork, afterFork):
    '''
    Checks if an adjacent pair of forks in a path are valid. Checks contig ID,
    strand and position. Returns True if valid pair and False otherwise.
    '''
    if beforeFork.is_Nfork() or afterFork.is_Nfork():
        return True
    
    if beforeFork.after_id() != afterFork.before_id(): return False
    if beforeFork.after_strand() != afterFork.before_strand(): return False
    if beforeFork.after_pos() >= afterFork.before_pos(): return False

    return True

def interleave_paths(path1, path2, q=True, preferRef=True):
    '''
    Takes two Paths. Sorts all Forks by query postion (if q=True) or
    ref position (if q=False). Assumes total ordering of positions.
    
    Removes a Fork from non-consistent Fork pair. Fork is removed such that
    more reference sequence is preferred (preferRef=True) or more query
    is preferred (preferRef=False).
    
    Returns an interleaved Path.    
    '''

    forks = path1.path + path2.path
    if q:
        forks = sorted(forks, key=lambda fork: fork.qpos)
    else:
        forks = sorted(forks, key=lambda fork: fork.rpos)

    i=0
    while i < len(forks)-1:
        if forks[i].is_switch_reference() and forks[i+1].is_switch_reference():
            if preferRef:  forks.pop(i+1)
            else:          forks.pop(i)
            continue
        if forks[i].is_switch_query() and forks[i+1].is_switch_query():
            if preferRef:  forks.pop(i)
            else:          forks.pop(i+1)

            continue
        i += 1
    
    interleave = Path()
    interleave.extend(forks)
    return interleave

def filter_paths(paths, tigId, endsOnly=True):
    '''
    Takes a list of Paths and returns the following two lists:
    1) Paths that contain tigId, sorted from smallest to largest
    2) Paths that do not start or end with tigId
    
    If endsOnly is True, only start/end Fork will be checked, otherwise
    every Fork in the Path will be checked.    
    '''
    
    validPaths,invalidPaths = [],[]

    for path in paths:
    
        if endsOnly:        
            if not path[0].has_id(tigId) and not path[-1].has_id(tigId):
                invalidPaths.append(path)
            else:        
                validPaths.append(path)      
        else:
            found = False
            for fork in path:
                if fork.has_id(tigId):
                    validPaths.append(path)
                    found = True
                    break
            if not found:
                invalidPaths.append(path)      
           
    validPaths = sorted(validPaths, key=lambda path: path.path_length(tigId))
    return (validPaths, invalidPaths)

def check_path_consistency(path, printInfo=True, waitForInput=False):
    '''
    Validates path. Validates ID, strand an position between every fork. Prints
    warning and returns False if path is inconsistent, otherwise returns True.
    Used for debugging and validation.
    '''
    if len(path) < 1: return True
    
    startFork = path[0]

    for endFork in path[1:]:
        if startFork.is_Nfork() or endFork.is_Nfork():
            startFork = endFork
            continue
                                
        if not startFork.after_id() == endFork.before_id():  
            print("Contig IDs do not match:")
            print(startFork)
            print(endFork)
            if waitForInput:
                input()
            return False
        
        if not startFork.after_strand() == endFork.before_strand():
            print("Strands do not match:")
            print(startFork)
            print(endFork)
            if waitForInput:
                input()
            return False

        if startFork.after_pos() > endFork.before_pos():
            print("Positional issue (going backwards):")
            print(startFork)
            print(endFork)
            if waitForInput:
                input()
            return False
        
        if startFork.before_id() == endFork.after_id() and \
        startFork.before_strand() == endFork.after_strand() and \
        startFork.before_pos() > endFork.after_pos():
            print("Positional issue (ending backwards):")
            print(startFork)
            print(endFork)
            
            if startFork.is_switch_reference():
                print("Might be intentional...")
            else:
                if waitForInput:
                    input()
                return False

        startFork = endFork

    return True

def clean_overlapping_forks(path, param, q=True):
    '''
    Identifies pairs of forks in path that represent overlapping segments and 
    combines into a single pair of forks. Returns a cleaned path.
    '''
    if len(path) < 1: return path
    
    cleanPath = Path()
    startFork = path[0]
    skip = False
    
    for i in range(1, len(path)):
        if skip:
            skip=False
            continue
        
        endFork = path[i]
        
        # check if pairs overlap,
        if (startFork.after_pos() > endFork.before_pos() or \
            startFork.before_pos() > endFork.after_pos()):
            
            if (q and startFork.is_switch_query()) or \
               (not q and startFork.is_switch_reference()):
                   
                    # combine gap filling region by removing two forks if it does overlap
                    if len(cleanPath) > 0:
                        startFork = cleanPath.pop()
                    elif len(path) > i+1:
                        startFork = path[i+1]
                        skip=True
                    
                    continue

        cleanPath.add_fork(startFork)
        startFork = endFork
        
    cleanPath.add_fork(endFork)
    return cleanPath

def clean_Nforks_if_consistent(path, q=False):
        
    #remove NNN that are not necessary
    toPop = []
    for i,fork in enumerate(path):
        if not fork.is_Nfork(): continue
        if i-1 < 0 or i+1 >= len(path): continue
        if not consistent_forks(path[i-1], path[i+1]): continue

        if q and not (path[i-1].is_switch_query() and path[i+1].is_switch_reference()):
            continue
        if not q and not (path[i-1].is_switch_reference() and path[i+1].is_switch_query()):
            continue
        
        toPop.append(i)

    for i in reversed(toPop): path.pop(i)
    return path

def clean_Nforks(path):
    '''
    Cleans up unnecessary NNNs from a path. Returns the cleaned path.
    '''   
    if len(path) < 1: return path

    #remove duplicated NNN
    toPop = []
    for i,fork in enumerate(path[:-1]):
        if fork.is_Nfork() and path[i+1].is_Nfork():
            toPop.append(i)

    for i in reversed(toPop): path.pop(i)
    
    #remove NNN from start or end
    if path[0].is_Nfork():  path.pop(0)
    if path[-1].is_Nfork(): path.pop()      
    
    return path

def make_consistent(path):
    '''
    Adds Nfork in path between invalid pairs of forks.
    Returns path with Nforks inserted.
    '''

    newPath = Path()
    for i in range(len(path)-1):
        newPath.add_fork(path[i])
        if not consistent_forks(path[i], path[i+1]):
            newPath.add_fork(get_Nfork())
            
    newPath.add_fork(path[-1])
    return newPath

def make_ends_consistent(path, lengthData, q=True):
    '''
    Adds Nfork at the start or end of the path if the start/end does not 
    correspond to the most extreme query (q=True) or reference (q=False)
    position within the path.
    Returns path with Nforks inserted.
    '''

    startNFlag, endNFlag = False, False
    if len(path) > 0:
        f1 = path[0]
        f2 = path[-1]
        
        firstPos = f1.get_pos_norm(lengthData, q)
        lastPos = f2.get_pos_norm(lengthData)
            
        for fork in path:
            if fork.is_Nfork(): continue
            pos = fork.get_pos_norm(lengthData, q) 
            if pos > lastPos:
                endNFlag = True
                break
            if firstPos is not None and pos < firstPos:
                startNFlag = True
                break
        
    if startNFlag: 
        path.add_fork_front(get_Nfork())
    if endNFlag:
        path.add_fork(get_Nfork())

    return path