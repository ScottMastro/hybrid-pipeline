import sys
sys.path.append("..")
import log.log as logger

import weld.path_helper as path_helper
from scaffold.scaffold_codes import fail_enum as f

def scaffold_pair(path1, path2, fork1, fork2, lengthData, param):
    '''
    Takes a pair of Paths and attempts to scaffold them.
    Returns a single scaffolded Path or None if scaffold attempt fails.
    '''

    fork1Flip = fork1.flip_strands(lengthData, makeCopy=True)
    fork2Flip = fork2.flip_strands(lengthData, makeCopy=True)

    idx1 = None
    for i,fork in enumerate(path1):
        if fork == fork1: 
            idx1 = i
            break
        if fork == fork1Flip:
            path1.flip_strands(lengthData)
            idx1 = len(path1)-1 - i
            
    idx2 = None
    for i,fork in enumerate(path2):
        if fork == fork2: 
            idx2 = i
            break
        if fork == fork2Flip:
            path2.flip_strands(lengthData)
            idx2 = len(path2)-1 - i

    if idx1 is None or idx2 is None: return None
    
    if path_helper.valid_fork_pair(path1[idx1], path2[idx2]):
        scaffold = path1.extract_subpath(None, idx1+1)
        scaffold.add_path(path2.extract_subpath(idx2, None))
        return scaffold
    
    return None


def keep_peeling(leftCf, rightCf, leftOrientation, rightOrientation, \
                 leftIdx, rightIdx, lengthData):
    nullResult = (None,None)
    leftPath, rightPath = leftCf.path, rightCf.path
    leftFork, rightFork = leftPath[leftIdx], rightPath[rightIdx]
    
    def flip(fork, orient):
        if orient == -1: 
            return fork.flip_strands(lengthData, makeCopy=True)
        else: return fork
    
    LEFT,RIGHT=-1,1
    def next_fork(fork, path, idx, orient, side):
        idx += orient * side
        if not 0 <= idx < len(path): return nullResult
        fork = flip(path[idx], orient)
        return (fork, idx)

    leftFork = flip(leftFork, leftOrientation)
    rightFork = flip(rightFork, rightOrientation)
    
    while leftFork.after_id() != leftCf.tigId:
        leftFork, leftIdx = next_fork(leftFork, leftPath, leftIdx, leftOrientation, LEFT)
        if leftFork is None: return nullResult
    
    while rightFork.before_id() != rightCf.tigId:
        rightFork, rightIdx = next_fork(rightFork, rightPath, rightIdx, rightOrientation, RIGHT)
        if rightFork is None: return nullResult

    preferLeft = (leftCf.span > rightCf.span)
    
    while True:
        if path_helper.valid_fork_pair(leftFork, rightFork):
            return (leftIdx, rightIdx)
        
        if not preferLeft:
            leftFork, leftIdx = next_fork(leftFork, leftPath, leftIdx, leftOrientation, LEFT)
        else:
            rightFork, rightIdx = next_fork(rightFork, rightPath, rightIdx, rightOrientation, RIGHT)
            
        if leftFork is None or rightFork is None:
            return nullResult


def can_scaffold(criticalForkList, i, unitigs, lengthData, param):
    if len(criticalForkList) < i+2:
        return (f.UNSURE, None)
    
    leftCf, rightCf = criticalForkList[i], criticalForkList[i+1]
    
    # check if both CriticalForks are from the same path (ie. already scaffolded)
    # (assumes each min/maxTigId can only be found in a single path)
    if leftCf.minTigId == rightCf.minTigId or leftCf.minTigId == rightCf.maxTigId:
        return (f.SAME_PATH, None)
    if leftCf.maxTigId == rightCf.minTigId or leftCf.maxTigId == rightCf.maxTigId:
        return (f.SAME_PATH, None)

    # check if one critical region completely overlaps the other
    # can be indication of a SV or a small contig that was supposed to fill a gap
    overlap = leftCf.overlap(rightCf)
    overlapPc = overlap / min(leftCf.span, rightCf.span)
    
    if overlapPc > 0.99:
        
        leftPathLength = leftCf.path.path_length()
        rightPathLength = rightCf.path.path_length()
        
        if leftPathLength <= rightPathLength and 1.0*overlap / leftPathLength > 0.97:
               return (f.FULL_OVERLAP, i)
        elif rightPathLength <= leftPathLength and 1.0* overlap / rightPathLength  > 0.97:
               return (f.FULL_OVERLAP, i+1)
        else:
            return (f.PARTIAL_OVERLAP, None)
                
    
    # check if a significant amount of data is lost if a scaffold happens
    # can be indication of a missassembly between the critical forks
    gainThresh = -10000
    dist = leftCf.dist(rightCf) + rightCf.span
    leftIdx, leftPeel, leftOrientation = leftCf.peel_to_max()
    rightIdx, rightPeel, rightOrientation = rightCf.peel_to_min()
    if leftIdx is None or rightIdx is None: 
        return (f.PEEL_FAIL, (leftIdx, rightIdx))
        
    gain = dist - (leftPeel + rightPeel)    
    if gain < gainThresh:
        #print(str(dist) + " vs " + str(gain))
        return (f.EXCESSIVE_TRIM, ((leftPeel + rightPeel), dist))
    
    if unitigs is not None:
        inv = leftCf.create_interval(rightCf)
        #print(inv)
        hits = unitigs.all_hits(inv, overlap=1.0)
        #print(hits)
        
        if len(hits) < 1:
            return (f.NO_UNITIG, None)
    
    leftFork = leftCf.maxFork
    rightFork = rightCf.minFork

    if leftOrientation == -1:
        leftFork = leftFork.flip_strands(lengthData, makeCopy=True)
    if rightOrientation == -1:
        rightFork = rightFork.flip_strands(lengthData, makeCopy=True)    
    
    if path_helper.valid_fork_pair(leftFork, rightFork):
        return (f.PASS, (leftFork, rightFork))
    
    #invalid fork pair, keep peeling!
    i1,i2 = keep_peeling(leftCf, rightCf, leftOrientation, rightOrientation, \
                 leftIdx, rightIdx, lengthData)
    
    if i1 is None or i2 is None:
        return (f.PEEL_FAIL, None)
        
    leftFork = leftCf.path[i1]
    rightFork = rightCf.path[i2]
        
    leftIdx, leftPeel, leftOrientation = leftCf.peel_right(leftFork)
    rightIdx, rightPeel, rightOrientation = rightCf.peel_left(rightFork)

    if leftOrientation == -1:
        leftFork = leftFork.flip_strands(lengthData, makeCopy=True)
    if rightOrientation == -1:
        rightFork = rightFork.flip_strands(lengthData, makeCopy=True)    

    if leftIdx is None or rightIdx is None: 
        return (f.PEEL_FAIL, (leftIdx, rightIdx))
        
    if path_helper.valid_fork_pair(leftFork, rightFork):
        print("overlap=",overlap)
        return (f.PASS, (leftFork, rightFork))

    return (f.UNSURE, None)
