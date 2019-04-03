from reversecomp import reverse_complement

class Path:
    def __init__(self):
        self.path = []
        self.supplementaryPath = []

    def id_on_end(self, id):
        if len(self.path) == 0:
            return None        
        return self.path[0].has_id(id) or self.path[-1].has_id(id)

    def add_fork(self, fork):
        self.path.append(fork)
  
    def pop(self, index=None):
        if len(self.path) < 1:
            return None        
        if index is None:
            return self.path.pop()
        
        return self.path.pop(index)

    def add_fork_front(self, fork):
        self.path.insert(0, fork)

    def add_path(self, other): 
        self.path.extend(other.path)
    
    def extend(self, forks): 
        self.path.extend(forks)

    def add_supplementary(self, path, startIndex, endIndex):
        self.supplementaryPath.append((path, startIndex, endIndex))
        
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
    
    def has_id(self, tigId):
        for fork in self.path:
            if fork.has_id(tigId):
                return True
        return False
    
    def last_fork(self):
        if len(self.path) == 0:
            return None
        return self.path[-1]
    
    def first_fork(self):
        if len(self.path) == 0:
            return None
        return self.path[0]
    
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
    
    def flip_strands(self, lengthData):
        self.path = self.path[::-1]
        for fork in self.path:
            fork.flip_strands(lengthData)
        
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
    
    def __getitem__(self, key):
         return self.path[key]
    def __len__(self):
         return len(self.path)
     
    def __repr__(self):
        return str(self.first_fork()) + "...\t..." + str(self.last_fork()) + "\t(" + str(len(self.path)) + ")"
    def __str__(self):
        return "\n".join([f.__repr__() for f in self.path])

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

class Fork:
    def __init__(self, qid, qpos, qstrand, rid, rpos, rstrand):
        self.qid = qid
        self.qpos = qpos
        self.qstrand = qstrand
        
        self.rid = rid
        self.rpos = rpos
        self.rstrand = rstrand
        self.switch = "q"
        
    def switch_query(self):
        self.switch = "q"

    def switch_reference(self):
        self.switch = "r"
        
    def before_id(self):
        return self.rid if self.switch == "q" else self.qid
    def after_id(self):
        return self.qid if self.switch == "q" else self.rid
    def before_pos(self):
        return self.rpos if self.switch == "q" else self.qpos
    def after_pos(self):
        return self.qpos if self.switch == "q" else self.rpos
    def before_strand(self):
        return self.rstrand if self.switch == "q" else self.qstrand
    def after_strand(self):
        return self.qstrand if self.switch == "q" else self.rstrand
    def before_switch(self):
        return "r" if self.switch == "q" else "q"
    def after_switch(self):
        return self.switch

    def flip_switch(self):
        self.switch = "r" if self.switch == "q" else "q"

    def flip_strands(self, lengthData, makeCopy=False):
        
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

    def has_id(self, id):
        return self.rid == id or self.qid == id

    def get_pos(self, q=True):
        if q: return self.qpos
        else: return self.rpos
        
    def get_pos_by_id(self, tigId):
        if self.qid == tigId:
            return self.qpos
        elif self.rid == tigId:
            return self.rpos
        else: return None
                
    def get_strand_by_id(self, tigId):
        if self.qid == tigId:
            return self.qstrand
        elif self.rid == tigId:
            return self.rstrand
        else: return None
        
    def __repr__(self):
        sq = str(self.qid) + " " + str(self.qpos) + " (" + str(self.qstrand) + ")"
        sr = str(self.rid) + " " + str(self.rpos) + " (" + str(self.rstrand) + ")"
        
        if self.switch == 'q': 
            return sr + " " + sq
        else:
            return sq + " " + sr
        
    def __str__(self):
        return self.__repr__()
    
def can_join_forks(beforeFork, afterFork):
    
    if not beforeFork.after_id() == afterFork.before_id():
        return False
    if not beforeFork.after_strand() == afterFork.before_strand():
        return False
    if beforeFork.after_strand() == 1:
        if not beforeFork.after_pos() < afterFork.before_pos():
            return False
    else:
        if not beforeFork.after_pos() > afterFork.before_pos():
            return False
    return True
    

def clean_strand(path, lengthData, verbose=False):

    strandedPaths = []
    toFlip = []

    flipStrand = False
    currentPath = Path()
    
    #if len(path) < 1: return currentPath
    
    startFork = path[0]
    flipStrand = False
    
    for i in range(1, len(path)):

        endFork = path[i]
    
        if verbose:
            print("New pair:\t" + str(startFork) + "\t" + str(endFork))

        if startFork.after_id() != endFork.before_id():           
            if verbose:
                print("Excluding fork because IDs do not match:")
                print(str(startFork) + " (try to match)")
                print(str(endFork) + " (excluded)")
            continue


        if not startFork.after_strand() == endFork.before_strand():
            
            currentPath.add_fork(startFork)
            strandedPaths.append(currentPath)
            currentPath = Path()
            toFlip.append(flipStrand)
            flipStrand = not flipStrand
            if verbose:
                print("Flipping at:\t" + str(endFork))

            startFork = endFork
            continue

        currentPath.add_fork(startFork)
        startFork = endFork
        
        if i == len(path)-1:
            currentPath.add_fork(endFork)

    strandedPaths.append(currentPath)
    toFlip.append(flipStrand)
    singleStrandPath = Path()

    for flip, forks in zip(toFlip, strandedPaths):
        if flip:
            forks.flip_strands(lengthData)

        singleStrandPath.add_path(forks)

    return singleStrandPath



def clean_path(path, lengthData, verbose=False):

    cleanPath = Path()

    if len(path) < 1: return cleanPath
    
    singleStrandPath = clean_strand(path, lengthData, verbose=False)
    startFork = singleStrandPath[0]
        
    for i in range(1, len(singleStrandPath)):

        endFork = singleStrandPath[i]
    
        if verbose:
            print("New pair:\t" + str(startFork) + "\t" + str(endFork))
            
        if startFork.after_id() != endFork.before_id():
            if verbose:
                print("Excluding fork because IDs do not match:")
                print(str(startFork) + " (try to match)")
                print(str(endFork) + " (excluded)")
            continue
        
        if not startFork.after_strand() == endFork.before_strand():
            print("Strand issue with forks!!!!")

        if startFork.after_pos() > endFork.before_pos():
            #this could happen if sets of NNNs are close

            if verbose:
                print("Excluding forks because positional issue")
                print(str(startFork) + " (start fork, excluding)")
                print(str(endFork) + " (end fork, excluding)")
                print("Popping last fork:")
                print(cleanPath[-1])

            startFork = cleanPath.pop()
            continue
            
        cleanPath.add_fork(startFork)
        #print("Keeping fork:\t" + str(startFork))
        
        startFork = endFork
        
        if i == len(singleStrandPath)-1:
            cleanPath.add_fork(endFork)

    return cleanPath