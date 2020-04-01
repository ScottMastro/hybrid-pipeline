import sys
sys.path.append("..")

from . import critical_forks
import structures.path_operations as pathops

class Result:
    #failure enums
    FULL_OVERLAP      = 1
    PARTIAL_OVERLAP   = 2
    EXCESSIVE_TRIM    = 3
    SAME_PATH         = 4
    PEEL_FAIL         = 5
    NO_UNITIG         = 6
    PASS              = 0
    UNSURE            = -1
    
    failDict = {UNSURE:"????", PASS:"pass", FULL_OVERLAP:"overlap (full)",
                PARTIAL_OVERLAP:"overlap (partial)", 
                EXCESSIVE_TRIM:"excessive trimming", SAME_PATH:"redundant",
                NO_UNITIG:"no unitig support", PEEL_FAIL:"could not find critical forks"}
    
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
    
    if pathops.consistent_forks(path1[idx1], path2[idx2]):
        scaffold = path1.extract_subpath(None, idx1+1)
        scaffold.add_path(path2.extract_subpath(idx2, None))
        scaffold.set_path_id(path1.pid + "," + path2.pid)
        print(scaffold.pid)
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
        if pathops.consistent_forks(leftFork, rightFork):
            return (leftIdx, rightIdx)
        
        if not preferLeft:
            leftFork, leftIdx = next_fork(leftFork, leftPath, leftIdx, leftOrientation, LEFT)
        else:
            rightFork, rightIdx = next_fork(rightFork, rightPath, rightIdx, rightOrientation, RIGHT)
            
        if leftFork is None or rightFork is None:
            return nullResult


def can_scaffold(criticalForkList, i, unitigs, lengthData, param):
    if len(criticalForkList) < i+2:
        return (Result.UNSURE, None)
    
    leftCf, rightCf = criticalForkList[i], criticalForkList[i+1]
    
    # check if both CriticalForks are from the same path (ie. already scaffolded)
    # (assumes each min/maxTigId can only be found in a single path)
    if leftCf.minTigId == rightCf.minTigId or leftCf.minTigId == rightCf.maxTigId:
        return (Result.SAME_PATH, None)
    if leftCf.maxTigId == rightCf.minTigId or leftCf.maxTigId == rightCf.maxTigId:
        return (Result.SAME_PATH, None)

    # check if one critical region completely overlaps the other
    # can be indication of a SV or a small contig that was supposed to fill a gap
    overlap = leftCf.overlap(rightCf)
    overlapPc = overlap / min(leftCf.span, rightCf.span)
    
    if overlapPc > 0.99:
        
        leftPathLength = leftCf.path.path_length()
        rightPathLength = rightCf.path.path_length()
        
        if leftPathLength <= rightPathLength and 1.0*overlap / leftPathLength > 0.97:
               return (Result.FULL_OVERLAP, i)
        elif rightPathLength <= leftPathLength and 1.0* overlap / rightPathLength  > 0.97:
               return (Result.FULL_OVERLAP, i+1)
        else:
            return (Result.PARTIAL_OVERLAP, None)
                
    
    # check if a significant amount of data is lost if a scaffold happens
    # can be indication of a missassembly between the critical forks
    gainThresh = -10000
    dist = leftCf.dist(rightCf) + rightCf.span
    leftIdx, leftPeel, leftOrientation = leftCf.peel_to_max()
    rightIdx, rightPeel, rightOrientation = rightCf.peel_to_min()
    if leftIdx is None or rightIdx is None: 
        return (Result.PEEL_FAIL, (leftIdx, rightIdx))
        
    gain = dist - (leftPeel + rightPeel)
    if gain < gainThresh:
        #print(str(dist) + " vs " + str(gain))
        return (Result.EXCESSIVE_TRIM, ((leftPeel + rightPeel), dist))
    
    if unitigs is not None:
        inv = leftCf.create_interval(rightCf)
        #print(inv)
        hits = unitigs.all_hits(inv, overlap=1.0)
        #print(hits)
        
        if len(hits) < 1:
            return (Result.NO_UNITIG, None)
    
    leftFork = leftCf.maxFork
    rightFork = rightCf.minFork

    if leftOrientation == -1:
        leftFork = leftFork.flip_strands(lengthData, makeCopy=True)
    if rightOrientation == -1:
        rightFork = rightFork.flip_strands(lengthData, makeCopy=True)    
    
    if pathops.consistent_forks(leftFork, rightFork):
        return (Result.PASS, (leftFork, rightFork))
    
    #invalid fork pair, keep peeling!
    i1,i2 = keep_peeling(leftCf, rightCf, leftOrientation, rightOrientation, \
                 leftIdx, rightIdx, lengthData)
    
    if i1 is None or i2 is None:
        return (Result.PEEL_FAIL, None)
        
    leftFork = leftCf.path[i1]
    rightFork = rightCf.path[i2]
        
    leftIdx, leftPeel, leftOrientation = leftCf.peel_right(leftFork)
    rightIdx, rightPeel, rightOrientation = rightCf.peel_left(rightFork)

    if leftOrientation == -1:
        leftFork = leftFork.flip_strands(lengthData, makeCopy=True)
    if rightOrientation == -1:
        rightFork = rightFork.flip_strands(lengthData, makeCopy=True)    

    if leftIdx is None or rightIdx is None: 
        return (Result.PEEL_FAIL, (leftIdx, rightIdx))
        
    if pathops.consistent_forks(leftFork, rightFork):
        print("overlap=",overlap)
        return (Result.PASS, (leftFork, rightFork))

    return (Result.UNSURE, None)


def scaffold_tigs(paths, tigId, unitigs, lengthData, param):
    '''
    Takes a tigId and a list of paths. Attempts to scaffold every Path that
    contains a Fork with tigId on the end.
    Returns a tuple of:
    1) a list of scaffolded Paths
    2) a list of Paths that were discarded (ex. redundantly overlaps with scaffold)
    '''
    
    #build critical forks and sort left to right
    criticalForkList, excludePaths = [], []
    for path in paths: 
        criticalForkList.extend( critical_forks.build(path, tigId, lengthData) )

    criticalForkList.sort(key=lambda cf: cf.minPos)

    #plot critical segments
    #plotter.plot_segments(criticalForkList, tigId, lengthData)
    #return (paths, [])

    i=0
    results = []
    while i < len(criticalForkList) -1:
        print(criticalForkList[i].minTigId, criticalForkList[i+1].minTigId) 
        result = can_scaffold(criticalForkList, i, unitigs, lengthData, param)
        resultCode, resultHint = result
        
        #one segment fully overlaps another, remove the smaller one
        if resultCode == Result.FULL_OVERLAP:
            excludePaths.append(criticalForkList[resultHint].path)
            criticalForkList.pop(resultHint)
            if resultHint == i and i > 0:
                results.pop()
                i=i-1
            continue
            
        #segments are coming from the same path with nothing inbetween, merge
        if resultCode == Result.SAME_PATH:
            criticalForkList[i].merge(criticalForkList[i+1])
            criticalForkList.pop(i+1)
            continue
        
        if resultCode == Result.UNSURE:
            print("!!!!!")
            #todo: does this happen?
            input()
        
        results.append(result)
        i = i + 1
        
    scaffoldCandidates = [cf.path for cf in criticalForkList]

    print([Result.failDict[resultCode] for resultCode,_ in results])
    print([cf.minTigId for cf in criticalForkList])

    #need to create this so that a path is not duplicated if it has
    #multiple critical forks
    redundantDict = dict()
    for i,pth in enumerate(scaffoldCandidates):
        for j,otherpth in enumerate(scaffoldCandidates):
            if i < j and pth == otherpth: 
                redundantDict[j] = i 
                print("WARNING: WATCH OUT FOR REDUNDANCY!")
    
    for i in redundantDict.keys(): scaffoldCandidates[i] = None
    #scaffoldCopy = copy.deepcopy(scaffoldCandidates)

    i=0
    redoChance = True
    while i < len(scaffoldCandidates) -1:
        redoChance = True

        if scaffoldCandidates[i] is None:
            scaffoldCandidates[i] = scaffoldCandidates[redundantDict[i]]
            scaffoldCandidates[redundantDict[i]] = None
        if scaffoldCandidates[i+1] is None:
            scaffoldCandidates[i+1] = scaffoldCandidates[redundantDict[i+1]]
            scaffoldCandidates[redundantDict[i+1]] = None

        resultCode, resultHint = results[i]
        if resultCode == Result.PASS:
            fork1, fork2 = resultHint
            path1, path2 = scaffoldCandidates[i], scaffoldCandidates[i+1]
            if path1 is not None and path2 is not None:
                scaffold = scaffold_pair(path1, path2, fork1, fork2, lengthData, param)
    
                #try again
                if scaffold is None and redoChance:
                    cf1 = critical_forks.build(scaffoldCandidates[i], tigId, lengthData)
                    cf2 = critical_forks.build(scaffoldCandidates[i+1], tigId, lengthData)
                    result = can_scaffold([cf1[0],cf2[0]], 0, unitigs, lengthData, param)
                    results[i] = result
                    redoChance = False
                    continue
                if scaffold is not None:
                    scaffoldCandidates[i] = None
                    scaffoldCandidates[i+1] = scaffold
            
        i = i+1
    
    scaffolds = [s for s in scaffoldCandidates if s is not None]
    print(len(scaffolds), " scaffolds returned.")
        
    return (scaffolds, excludePaths)

