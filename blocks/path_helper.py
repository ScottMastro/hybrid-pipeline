from reversecomp import reverse_complement
import log
import matplotlib.pyplot as plt

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
    def __init__(self, qid, qpos, qstrand, rid, rpos, rstrand, Nfork=False):
        self.qid = qid
        self.qpos = qpos
        self.qstrand = qstrand
        
        self.rid = rid
        self.rpos = rpos
        self.rstrand = rstrand
        self.switch = "q"
        self.Nfork = Nfork
             
    def is_Nfork(self):
        return self.Nfork
        
    def switch_query(self):
        self.switch = "q"
    def is_switch_query(self):
        return self.switch == "q"

    def switch_reference(self):
        self.switch = "r"
    def is_switch_reference(self):
        return self.switch ==  "r"

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
        
        if self.Nfork:
            if makeCopy: return get_Nfork()
            return
        
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
    
    
def get_Nfork():
    return Fork("NNN", -1, 0, "NNN", -1, 0, Nfork=True)
def get_Npath():
    p = Path()
    p.add_fork(get_Nfork())
    return p
    
def can_join_forks(beforeFork, afterFork):
    
    if not beforeFork.after_id() == afterFork.before_id():
        return False
    if not beforeFork.after_strand() == afterFork.before_strand():
        return False
    if beforeFork.after_strand() == 1:
        if not beforeFork.after_pos() < afterFork.before_pos():
            return False
    else:
        if not beforeFork.after_pos() < afterFork.before_pos():
            return False
    return True

def get_tig_ids(path, source=None):
                       #'r' for ref, 'q' for query, None = both                
    tigIds = set()
    for fork in path:
        if source is None or source == 'r':
            tigIds.add(fork.rid)
        if source is None or source == 'q':
            tigIds.add(fork.qid)

    return tigIds


def path_length(path):

    pathSum = 0

    if len(path) >= 2: 
    
        prevFork = path[0]
        
        for fork in path[1:]:
            pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
            prevFork = fork
        
    return pathSum


def plot_length(path, lengthData, printout=True):

    pathSum = 0
    x=[]
    y=[]
    if len(path) >= 2: 
        prevFork = path[0]
        
        for fork in path[1:]:
            print(fork)
            pos = fork.qpos
            if fork.qstrand == -1:
                pos = lengthData[str(fork.qid)] - pos
            pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
            
            x.append(pos)
            y.append(pathSum)
            if printout:
                print(str(pathSum) + " / " + str(pos) + "  " + str(round(pathSum/pos,2)))
            
            prevFork = fork
    
    plt.scatter(x,y)


def path_overlap(path1, path2, lengthData, source=None):
                                           #'r' for ref, 'q' for query, None = both
    def normalized_pos(fork, tigId):
        pos = fork.get_pos_by_id(tigId)
        if pos is None: return None
        if fork.get_strand_by_id(tigId) == -1:
            pos = lengthData[str(tigId)] - pos
        return pos
    
    def fill_dict(path):
        starts = dict()
        ends = dict()
        for fork in path:
            if source is None :
                tigIds = [fork.before_id(), fork.after_id()]
            else:
                tigIds = ([fork.rid] if source == 'r' else []) + \
                 ([fork.qid] if source == 'q' else [])
                 
            for tigId in tigIds:
                pos = normalized_pos(fork, tigId)
                tigId = str(tigId)
                if tigId not in starts:
                    starts[tigId] = pos
                    ends[tigId] = pos
                else:
                    if pos < starts[tigId]: starts[tigId] = pos
                    if pos > ends[tigId]: ends[tigId] = pos
        return (starts,ends)
               
    starts1, ends1 = fill_dict(path1)
    starts2, ends2 = fill_dict(path2)


    overlap = 0.0
    total1 = 0.0
    total2 = 0.0
    for tigId in set(starts1.keys()).union(set(starts2.keys())):
        if tigId in starts1 and tigId in starts2:
            start = max(starts1[tigId], starts2[tigId])
            end = min(ends1[tigId], ends2[tigId])
            overlap = overlap + max(0, end - start)
            #print(tigId)
            #print(str(start) + " - " + str(end))
            #print (overlap)
            
        if tigId in starts1:
            total1 = total1 + abs(ends1[tigId] - starts1[tigId])
        if tigId in starts2:
            total2 = total2 + abs(ends2[tigId] - starts2[tigId])
        
    return (overlap/(total1 +1), overlap/(total2 +1))
        
def clean_strand(path, lengthData, param):

    strandedPaths = []
    toFlip = []

    flipStrand = False
    currentPath = Path()
    
    if len(path) < 1: return currentPath
    
    log.out("Fixing strand orientation.", 2, param)

    startFork = path[0]
    flipStrand = False
    
    for i in range(1, len(path)):

        endFork = path[i]
        log.out("New pair:\t" + str(startFork) + "\t" + str(endFork), 3, param)

        if startFork.after_id() != endFork.before_id():           
            #this should not happen if all alignments were successful !
            log.out("Excluding fork because IDs do not match:", 2, param)
            log.out(str(startFork) + " (start fork, try to match)", 2, param)
            log.out(str(endFork) + " (end fork, excluded)", 2, param)
            continue


        if not startFork.after_strand() == endFork.before_strand():
            
            currentPath.add_fork(startFork)
            strandedPaths.append(currentPath)
            currentPath = Path()
            toFlip.append(flipStrand)
            flipStrand = not flipStrand
            log.out("Flipping strand at:" + str(endFork), 2, param)
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
        
    log.out("Strand orientation fixing complete.", 2, param)

    return singleStrandPath



def clean_path(path, lengthData, param):

    if len(path) < 1: return Path()

    log.out("Cleaning path.", 1, param)
    log.out("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", 1, param)
    
    cleanPath = Path()

    report = log.ReportSet(log.PATH_CLEAN_ATTEMPT)
    param.add_report(report)
    
    singleStrandPath = clean_strand(path, lengthData, param)
    startFork = singleStrandPath[0]
    skip = False
    
    for i in range(1, len(singleStrandPath)):
        if skip:
            skip=False
            continue
        
        endFork = singleStrandPath[i]
    
        log.out("New pair:\t" + str(startFork) + "\t" + str(endFork), 3, param)
            
        if startFork.after_id() != endFork.before_id():      
            #this should not happen after cleanStrand() !
            log.out("Excluding fork because IDs do not match:", 2, param)
            log.out(str(startFork) + " (start fork, try to match)", 2, param)
            log.out(str(endFork) + " (end fork, excluded)", 2, param)
            continue
        
        if not startFork.after_strand() == endFork.before_strand():
            #this should not happen after cleanStrand() !
            log.out("Excluding fork because strands do not match:", 2, param)
            log.out(str(startFork) + " (start fork, try to match)", 2, param)
            log.out(str(endFork) + " (end fork, excluded)", 2, param)
            continue

        if startFork.after_pos() > endFork.before_pos():
            log.out("Excluding forks because positional issue:", 2, param)

            #this could happen if sets of NNNs are close
            #if startFork.is_switch_query():
            if startFork.switch == 'q':

                log.out(str(startFork) + " (start fork, excluding)", 2, param)
                log.out(str(endFork) + " (end fork, excluding)", 2, param)
                
                if len(cleanPath) > 0:
                    log.out("Popping last fork: " + str(cleanPath[-1]), 2, param)
                    startFork = cleanPath.pop()
                elif len(singleStrandPath) > i+1:
                    startFork = singleStrandPath[i+1]
                    skip=True
                continue
            
            #this could happen if blocks overlap
            #if startFork.is_switch_reference():
            if startFork.switch == 'r':

                if len(cleanPath) > 0:
                    log.out("Popping and removing last fork: " + str(cleanPath[-1]), 2, param)
                    cleanPath.pop()

                log.out(str(startFork) + " (start fork, excluding)", 2, param)
                log.out(str(endFork) + " (end fork, try to match)", 2, param)
                
                startFork = endFork
                continue

            
            
        cleanPath.add_fork(startFork)
        log.out("Keeping fork:" + str(startFork), 3, param)
        
        startFork = endFork
        
        if i == len(singleStrandPath)-1:
            cleanPath.add_fork(endFork)
            
    log.out("Path cleaning complete.", 1, param)
    log.out("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", 1, param)

    return cleanPath