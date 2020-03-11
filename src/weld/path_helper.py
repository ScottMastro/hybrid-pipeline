from weld.paths import Path
from weld.paths import get_Nfork
import log as logger

REORIENTATION_FILE_NAME = "inversions.txt"

def interleave_paths(path1, path2):
    #assumes total ordering of query positions
    forks = path1.path + path2.path
    forks = sorted(forks, key=lambda fork: fork.qpos)
    
    #favours reference
    i=0
    while i < len(forks)-1:
        if forks[i].is_switch_reference() and forks[i+1].is_switch_reference():
            forks.pop(i+1)
            continue
        if forks[i].is_switch_query() and forks[i+1].is_switch_query():
            forks.pop(i)
            continue
        i += 1
    
    interleave = Path()
    interleave.extend(forks)
    return interleave

def filter_paths(paths, tigId, ends=True):
    '''
    Takes a list of Paths and returns the following two lists:
    1) Paths that contain tigId, sorted from smallest to largest
    2) Paths that do not start or end with tigId
    
    If ends is True, only start/end Fork will be checked, otherwise
    every Fork in the Path will be checked.    
    '''
    
    validPaths,invalidPaths = [],[]

    for path in paths:
    
        if ends:        
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


def clean_NNNs(path):
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
    
    #remove NNN that are not necessary
    toPop = []
    for i,fork in enumerate(path):
        if fork.is_Nfork() and \
            path[i-1].after_id() == path[i+1].before_id() and \
            path[i-1].after_strand() == path[i+1].before_strand() and \
            path[i-1].is_switch_reference() and path[i+1].is_switch_query() and \
            path[i-1].after_pos() <= path[i+1].before_pos():
                toPop.append(i)

    for i in reversed(toPop): path.pop(i)
    
    return path

def check_path(path):
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
            #input()
            return False
        
        if not startFork.after_strand() == endFork.before_strand():
            print("Strands do not match:")
            print(startFork)
            print(endFork)
            #input()
            return False

        if startFork.after_pos() > endFork.before_pos():
            print("Positional issue (going backwards):")
            print(startFork)
            print(endFork)
            #input()
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
                #input()
                return False

        startFork = endFork

    return True

def valid_fork_pair(beforeFork, afterFork):
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

def clean_overlapping_forks(path, param):
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
            startFork.before_pos() > endFork.after_pos()) and \
            startFork.is_switch_query():
           
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

def fix_path_orientation(path1, path2, lengthData, param, trimEnds=False):
    '''
    Takes in a pair of paths and compares strandedness of the last fork of path1
    and the first fork of path 2. Flips one or both of the paths to correct strandedness, 
    if necessary. Setting trimEnds=True will evaluate the second last fork of path1
    and second fork in path2 instead. Returns the corrected paths.
    '''
    index = 1 if trimEnds else 0
    
    fork1, fork2 = path1.tail(index), path2.head(index)
    fork1Flip = path1.head(index).flip_strands(lengthData, makeCopy=True)
    fork2Flip = path2.tail(index).flip_strands(lengthData, makeCopy=True)

    if not valid_fork_pair(fork1, fork2):

        print('attempting to flip path...')
        print("------")
        print(fork1)
        print(fork2)
        print("------")

        def path_qinfo(path):
            return str(path[0].qid) + ":" + str(path[0].qpos) + "-" + str(path[-1].qpos) + ":" + str(path[0].qstrand)
        def forks_rinfo(fork1, fork2):
            return str(fork1.rid) + ":" + str(fork1.qpos) + "-" + str(fork2.qpos) + ":" + str(fork1.qstrand)

        path1Info = path_qinfo(path1)
        path2Info = path_qinfo(path2)
        
        if valid_fork_pair(fork1, fork2Flip):
            print("flipping right")
            path2.flip_strands(lengthData)
            info = [path2Info, "LEFT", forks_rinfo(fork1, fork2Flip)]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)

        elif valid_fork_pair(fork1Flip, fork2):
            print("flipping left")
            path1.flip_strands(lengthData)
            info = [path1Info, "RIGHT", forks_rinfo(fork1Flip, fork2)]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)

        elif valid_fork_pair(fork1Flip, fork2Flip):
            print("flipping both")
            path1.flip_strands(lengthData)
            path2.flip_strands(lengthData)
            info = [path1Info, "RIGHT", forks_rinfo(fork1Flip, fork2Flip)]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)
            info = [path2Info, "LEFT", forks_rinfo(fork1Flip, fork2Flip)]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)

        else:
            print("no suitable flip found.")
            info = [path1Info, "FAIL", path2Info]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)


        logger.FileLogger().flush(REORIENTATION_FILE_NAME)

    return(path1, path2)

def add_Nforks(path, lengthData):
    '''
    Adds Nfork in path between invalid pairs of forks. Adds Nfork at the start
    or end of the path if the start/end does not correspond to the most extreme
    query position within the path. Returns path with Nforks inserted.
    '''
    #print(path)
    startNFlag, endNFlag = False, False
    if len(path) > 0:
        f1 = path[0]
        f2 = path[-1]
        
        firstqPos = f1.get_pos_norm(lengthData, q=True)
        lastqPos = f2.get_pos_norm(lengthData, q=True)
            
        for fork in path:
            if fork.is_Nfork(): continue
            qpos = fork.get_pos_norm(lengthData, q=True) 
            if qpos > lastqPos:
                print("ending with NNN")
                endNFlag = True
                break
            if firstqPos is not None and qpos < firstqPos:
                print("starting with NNN")
                startNFlag = True
                break
        
    if startNFlag: 
        path.add_fork_front(get_Nfork())
        print(f1)

    if endNFlag:
        path.add_fork(get_Nfork())
        print(f2)


    # add Nfork where necessary
    # todo: make sure there's no other possible fix here
    # does this happen?
    newPath = Path()
    for i in range(len(path)-1):
        newPath.add_fork(path[i])
        if not valid_fork_pair(path[i], path[i+1]):
            newPath.add_fork(get_Nfork())
            print("adding Nfork between:")
            print(path[i])
            print(path[i+1])
            #input()
            
    newPath.add_fork(path[-1])
    return newPath

def path_length(path, tigId=None):
    '''
    Calculates the total length of the path. If tigId is provided,
    only counts segments corresponding to tigId.
    '''

    if len(path) < 2: return 0
    pathSum = 0

    prevFork = path[0]
    for fork in path[1:]:
        if not fork.is_Nfork() and not prevFork.is_Nfork():
            
            if tigId is None or (prevFork.after_id() == tigId and fork.before_id() == tigId):
                pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
                
        prevFork = fork
        
    return pathSum


def path_overlap(path1, path2, lengthData, source=None, printInfo=False):
                                           #'r' for ref, 'q' for query, None = both
    
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
                pos = fork.get_pos_by_id_norm(tigId, lengthData)
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
            
            if printInfo:
                start = max(starts1[tigId], starts2[tigId])
                end = min(ends1[tigId], ends2[tigId])

                print(tigId + ": " + str(start) + " - " + str(end) + \
                      "(" + str(overlap) + ")")
            
        if tigId in starts1:
            total1 = total1 + abs(ends1[tigId] - starts1[tigId])
        if tigId in starts2:
            total2 = total2 + abs(ends2[tigId] - starts2[tigId])
        
    return (overlap/max(total1, 1), overlap/max(total2, 1))


def path_overlap2(path1, path2, lengthData, overlapOnly, source=None, printInfo=False):
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
    whichEnd1 = []
    whichEnd2 = [] 
    #total1 = 0.0
    #total2 = 0.0
    
    tigs = []
    posA = []
    posB = []
    sizeA = []
    sizeB = []
    s1 = list(starts1.keys())
    s2 = list(starts2.keys())
    for tigId in set(starts1.keys()).intersection(set(starts2.keys())):
        #print(tigId)
        if s1[0] == tigId:
            whichEnd1.append('Head')
        elif s1[-1] == tigId:
            whichEnd1.append('Tail')
        else:
            whichEnd1.append('Middle')
        if s2[0] == tigId:
            whichEnd2.append('Head')
        elif s2[-1] == tigId:
            whichEnd2.append('Tail')
        else:
            whichEnd2.append('Middle')
        tigs.append(tigId)
        posA.append(str(starts1[tigId]) + '-' + str(ends1[tigId]))
        posB.append(str(starts2[tigId]) + '-' + str(ends2[tigId]))
        sizeA.append(ends1[tigId] - starts1[tigId])
        sizeB.append(ends2[tigId] - starts2[tigId])
            
        
        if tigId in starts1 and tigId in starts2:
            start = max(starts1[tigId], starts2[tigId])
            end = min(ends1[tigId], ends2[tigId])
            overlap = overlap + max(0, end - start)
            
            if printInfo:
                #start = max(starts1[tigId], starts2[tigId])
                #end = min(ends1[tigId], ends2[tigId])
    
                print("path1: " + path1 + ", path2: " + path2 + ", " + tigId + ": " + str(start) + " - " + str(end) + \
                      "(" + str(overlap) + ")")
            
        #if tigId in starts1:
            #total1 = total1 + abs(ends1[tigId] - starts1[tigId])
        #if tigId in starts2:
            #total2 = total2 + abs(ends2[tigId] - starts2[tigId])
        
        
    if overlapOnly:
        return overlap
    return whichEnd1, whichEnd2, tigs, posA, posB, sizeA, sizeB


def get_tig_ids(path, source=None):
                       #'r' for ref, 'q' for query, None = both                
    tigIds = set()
    for fork in path:
        if source is None or source == 'r':
            tigIds.add(fork.rid)
        if source is None or source == 'q':
            tigIds.add(fork.qid)

    return tigIds