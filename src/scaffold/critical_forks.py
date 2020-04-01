from pybedtools import Interval

LEFT, RIGHT ="l", "r"
class CriticalForks:
    
    def __init__(self, minFork, maxFork, path, tigId, lengthData, q=False): 
        self.minFork = minFork
        self.maxFork = maxFork
        self.path = path
        self.tigId = tigId
        self.tigLength = lengthData[tigId]
        self.minPos = minFork.get_pos_by_id_norm(tigId, lengthData)
        self.maxPos = maxFork.get_pos_by_id_norm(tigId, lengthData)
        self.minTigId = minFork.get_id(not q)
        self.maxTigId = maxFork.get_id(not q)
        self.span = abs(self.maxPos - self.minPos)

    def merge(self, other):
        self.maxFork = other.maxFork
        self.maxPos = other.maxPos
        self.maxTigId = other.maxTigId
        self.span = abs(self.maxPos - self.minPos)

    def overlap(self, other):
        a1, a2 = self.minPos, self.maxPos
        b1, b2 = other.minPos,other.maxPos
        
        if b1 > a2: return 0
        
        os, oe = max(a1,b1), min(a2,b2)
        return oe - os
    
    def create_interval(self, other):
        return Interval(self.tigId, self.maxPos, other.minPos)
    
    def dist(self, other):
        a1, a2 = self.minPos,  self.maxPos
        b1, b2 = other.minPos, other.maxPos
        
        return max(0, b1 - a2)

    def _peel_forward(self, targetFork):
        length = 0
        prevFork = None
        for i,fork in enumerate(self.path):
            if prevFork is not None \
                and not fork.is_Nfork() and not prevFork.is_Nfork():
                length += abs(prevFork.after_pos() - fork.before_pos())    
    
            if fork == targetFork: return i, length
            prevFork = fork
    
    def _peel_reverse(self, targetFork):
        length = 0
        prevFork = None
        for i,fork in enumerate(reversed(self.path)):
            if prevFork is not None \
                and not fork.is_Nfork() and not prevFork.is_Nfork():
                length += abs(prevFork.before_pos() - fork.after_pos())    
            
            if fork == targetFork: return len(self.path)-1 - i, length
            prevFork = fork

    def _peel(self, targetFork, side):
        minIdx, maxIdx = None, None
        for i,fork in enumerate(self.path):
            if self.minFork == fork: minIdx = i
            if self.maxFork == fork: maxIdx = i
            if minIdx is not None and maxIdx is not None:
                break

        if minIdx is None or maxIdx is None: return (None, None, None)
        orientation = 1
        if side == LEFT:
            if minIdx < maxIdx: peelResult = self._peel_forward(targetFork)
            else:
                peelResult = self._peel_reverse(targetFork)
                orientation = -1
        if side == RIGHT:
            if minIdx < maxIdx: peelResult = self._peel_reverse(targetFork)
            else: 
                peelResult = self._peel_forward(targetFork)
                orientation = -1
            
        return (peelResult[0], peelResult[1], orientation)
            
    def peel_to_min(self): return self._peel(self.minFork, LEFT)
    def peel_to_max(self): return self._peel(self.maxFork, RIGHT)      
    def peel_left(self, targetFork): return self._peel(targetFork, LEFT)      
    def peel_right(self, targetFork): return self._peel(targetFork, LEFT)      

        
    def peel_max_reverse(self): return self.peel_reverse(self.maxFork)

    def __repr__(self):
        return self.tigId + " " + str(self.minPos) + " - " + str(self.maxPos) 
        
    def __str__(self):
        return "Critical forks for " + self.tigId + "\n" + \
        str(self.minPos) + " - " + str(self.maxPos) + "\n" + \
        "(" + self.minTigId + " / " + self.maxTigId + ")\n" + \
        "Span = " + str(self.span)


def build(path, tigId, lengthData):
    forkSets = []
    forkSet = []
    
    for fork in path:
        if fork.is_Nfork():
            continue
        elif fork.has_id(tigId):
            forkSet.append(fork)
        elif len(forkSet) > 0:
            forkSets.append(forkSet)
            forkSet=[]
            
    if len(forkSet) > 0: forkSets.append(forkSet)

    criticalForks = []
    
    def pos(fork): return fork.get_pos_by_id_norm(tigId, lengthData)
    
    for forkSet in forkSets:
        maxFork, minFork = forkSet[0], forkSet[0]
        
        for fork in forkSet:
            if pos(fork) < pos(minFork): minFork = fork
            if pos(fork) > pos(maxFork): maxFork = fork
        
        criticalForks.append(CriticalForks(minFork, maxFork, path, tigId, lengthData))     
        
    return criticalForks