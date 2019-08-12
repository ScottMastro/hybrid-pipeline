import matplotlib.pyplot as plt
from paths import Path 
from paths import get_Nfork

def valid_fork_pair(beforeFork, afterFork):
    '''
    Checks if an adjacent pair of forks in a path are valid. Checks contig ID,
    strand and position. Returns True if valid pair and False otherwise.
    '''
    if beforeFork.is_Nfork() or afterFork.is_Nfork():
        return True
    
    if beforeFork.after_id() != afterFork.before_id(): return False
    if beforeFork.after_strand() != afterFork.before_strand(): return False
    if beforeFork.after_pos() > afterFork.before_pos(): return False
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

        if valid_fork_pair(fork1, fork2Flip):
            print("flipping right")
            path2.flip_strands(lengthData)

        elif valid_fork_pair(fork1Flip, fork2):
            print("flipping left")
            path1.flip_strands(lengthData)
                
        elif valid_fork_pair(fork1Flip, fork2Flip):
            print("flipping both")
            path1.flip_strands(lengthData)
            path2.flip_strands(lengthData)
        else:
            print("no suitable flip found.")

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
        
        if not f1.is_Nfork():
            firstqPos = f1.qpos if f1.qstrand == 1 else lengthData[f1.qid] - f1.qpos
        else:
            firstqPos = None
        if not f2.is_Nfork():
            lastqPos = f2.qpos if f2.qstrand == 1 else lengthData[f2.qid] - f2.qpos
        else:
            lastqPos = None
            
        for fork in path:
            if fork.is_Nfork(): continue
            qpos = fork.qpos if fork.qstrand == 1 else lengthData[fork.qid] - fork.qpos
            
            if lastqPos is not None and qpos > lastqPos:
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


def join_paths(paths, lengthData, param):
    
    if len(paths) < 1: return Path()

    for i in range(len(paths)-1):
        p1 = paths[i]
        p2 = paths[i+1]
        
        if p1[-1].is_Nfork() or p2[0].is_Nfork(): continue
        
        p2Flip = p2[-1].flip_strands(lengthData, makeCopy=True)
        p1Flip = p1[0].flip_strands(lengthData, makeCopy=True)

        
        if p1[-1].after_id() == p2[0].before_id() and \
            p1[-1].after_strand() != p2[0].before_strand():

            print('flipping block path...')
            print(p1[-1])
            print(p2[0])
            

            if p1[-1].after_id() == p2Flip.before_id() and \
                p1[-1].after_pos() < p2Flip.before_pos():
                    print("flipping right")
                    paths[i+1].flip_strands(lengthData)
                    
                    if i+2 >= len(paths): continue
                    p3 = paths[i+2]
                    if p3 is None or len(p3) < 1 or p3[0].is_Nfork(): 
                        continue
                    if p3[-1].after_id() != paths[i+1][-1].before_id():
                        paths[i+1].add_fork(get_Nfork())

                    
            elif p1Flip.after_id() == p2[0].before_id() and \
                p1Flip.after_pos() < p2[0].before_pos():
                    print("flipping left")
                    paths[i].flip_strands(lengthData)
                    
                    if i < 1: continue
                    p0 = paths[i-1]
                    if p0 is None or len(p0) < 1 or p0[-1].is_Nfork(): 
                        continue
                    if p0[-1].after_id() != paths[i][0].before_id():
                        paths[i].add_fork_front(get_Nfork())
            else:
                    print('skipped')
                    #input()
     
                    
        elif p1[-1].after_id() == p2[0].before_id() and \
            p1[-1].after_strand() == p2[0].before_strand():
        
            if p1Flip.after_id() == p2Flip.before_id() and \
                p1Flip.after_pos() < p2Flip.before_pos():
                    print("flipping both")
                    paths[i].flip_strands(lengthData)
                    paths[i+1].flip_strands(lengthData)
            
                    if i < 1: continue
                    p0 = paths[i-1]
                    if p0 is None or len(p0) < 1 or p0[-1].is_Nfork(): 
                        continue
                    if p0[-1].after_id() != paths[i][0].before_id():
                        paths[i].add_fork_front(path_helper.get_Nfork())
                        
                    if i > len(paths)-1: continue
                    p3 = paths[i+2]
                    if p3 is None or len(p3) < 1 or p3[0].is_Nfork(): 
                        continue
                    if p3[-1].after_id() != paths[i+1][-1].before_id():
                        paths[i+1].add_fork(path_helper.get_Nfork())

                        

    path = paths[0]
    lastIdx=0
    
    for blockPath in paths[1:]:
        
        if blockPath is None or len(blockPath) < 1: continue
     
        if path[-1].is_Nfork() or blockPath[0].is_Nfork():
            lastIdx=len(path)
            path.add_path(blockPath)
            continue
    
        #trying to switch between two different reference tigs
        #remove both forks
        if path[-1].after_id() != blockPath[0].before_id() and \
            path[-1].is_switch_reference() and blockPath[0].is_switch_query():
                #print("popping:")
                #print(path[-1])
                #print(blockPath[0])
                path.pop()
                blockPath.pop(0)
                
                if len(blockPath) < 1: continue
                if path[-1].is_switch_query() and blockPath[0].is_switch_reference() and \
                    path[-1].after_strand() != blockPath[0].before_strand():
                        path.add_fork(get_Nfork())
                        path.add_path(blockPath)
                        continue

    
        if path[-1].after_id() != blockPath[0].before_id():

            
            print('id issue')
            print(path[-1])
            print(blockPath[0])
            
           # if path[-1].is_switch_reference():
           #     print("removing last fork (resolved):")
           #     print(path[-1])
           #     path.pop()
            
            #input()
    
        if path[-1].after_strand() != blockPath[0].before_strand():
                        
            print('handle strand between blocks...')
            
            print(path[-1])
            print(blockPath[0])

            '''
            flippedFork = blockPath[-1].flip_strands(lengthData, makeCopy=True)
            if path[-1].after_id() == flippedFork.before_id() and \
                path[-1].after_pos() < flippedFork.before_pos():
                    print("flipping")
                    blockPath.flip_strands(lengthData)
            else:
                print('lastIdx='+str(lastIdx))
                print(path)
                flippedFork = path[lastIdx].flip_strands(lengthData, makeCopy=True)
                if flippedFork.after_id() == blockPath[0].before_id() and \
                    flippedFork.after_pos() < blockPath[0].before_pos():
                        print("flipping prev")
                        toFlip = []
                        while len(path) > lastIdx: toFlip.append(path.pop())
                        for flip in toFlip: 
                            flip.flip_strands(lengthData)
                            path.add_fork(flip)
                else:
                    print('unhandled')
                    input()
            '''
                
            if len(blockPath) > 1 and abs(blockPath[0].qpos - blockPath[-1].qpos) <= 6001:
                print(abs(blockPath[0].qpos - blockPath[-1].qpos))
                print("removing")
                continue
               
            #input()

        if path[-1].after_pos() > blockPath[0].before_pos():
            print('positional issue')
            print(path[-1])
            print(blockPath[0])
            #input()
 
        lastIdx=len(path)
        path.add_path(blockPath)
        
    return path



def join_path_list(paths, lengthData, param):
    
    while len(paths) > 0:
        if paths[0] is None or len(paths[0]) < 1: paths = paths[1:]
        else: break
    if len(paths) < 1: return Path()

    path = paths[0]
    for blockPath in paths[1:]:
        
        if blockPath is None or len(blockPath) < 1: continue
        
        if path[-1].after_id() != blockPath[0].before_id():
            print('id issue')
            print(path[-1])
            print(blockPath[0])
            
            if path[-1].is_switch_reference():
                print("removing last fork (resolved):")
                print(path[-1])
                path.pop()
            
            #input()
    
        if path[-1].after_strand() != blockPath[0].before_strand():
                        
            print('handle strand between blocks...')
            
            print(path[-1])
            print(blockPath[0])

            print("flipping")
            blockPath.flip_strands(lengthData)
            
         #   if abs(blockPath[0].qpos - blockPath[-1].qpos) <= 6001:
         #       print(abs(blockPath[0].qpos - blockPath[-1].qpos))
          #      print("removing")
          #      continue
               
            #input()

            
           # if path[-1].is_switch_reference() and len(blockPath) <= 2:
           #     print("removing")
               # continue
            
            #input()
            
        if path[-1].after_pos() > blockPath[0].before_pos():
            print('positional issue')
            print(path[-1])
            print(blockPath[0])
            #input()

        path.add_path(blockPath)
        
    return path












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
            if not fork.is_Nfork() and not prevFork.is_Nfork():
                pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
            prevFork = fork
        
    return max(1, pathSum)

def path_length_id(path, tigId):

    if len(path) >= 2: 
        first = None
        last = None
        for i in range(0, len(path)):
            if str(path[i].before_id()) == str(tigId):
                if first is None:
                    first = path[i].before_pos()
            if str(path[i].after_id()) == str(tigId):
                last = path[i].after_pos()
        pathSum = abs(last - first)
        '''
        prevFork = path[0]
        
        for fork in path[1:]:
            if not fork.is_Nfork() and not prevFork.is_Nfork():
                if prevFork.after_id() == tigId and fork.before_id() == tigId:
                    pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
            prevFork = fork
        '''
    return max(1, pathSum)

def plot_length_all(paths, lengthData, outputPath):
    
    outputPath="/media/scott/HDD/sickkids/NA24385/correlation/"
    for path in paths:
        plot_length(path, lengthData, printout=False, outputPath=outputPath)

def plot_length(path, lengthData, printout=True, outputPath=None):

    pathSum = 0
    x=[]
    y=[]
    if len(path) >= 2: 
        prevFork = path[0]
        
        for fork in path[1:]:
            if not fork.is_Nfork() and not prevFork.is_Nfork():
                pos = fork.qpos
                if fork.qstrand == -1:
                    pos = lengthData[str(fork.qid)] - pos
                pathSum = pathSum + abs(prevFork.after_pos() - fork.before_pos())
                
                x.append(pos)
                y.append(pathSum)
                
            if printout:
                print(fork)
                print(str(pathSum) + " / " + str(pos) + "  " + str(round(pathSum/pos,2)))
            
            prevFork = fork
            
    plt.scatter(x,y)
    if outputPath is not None:
        plt.savefig(outputPath + str(path[0].qid) + ".png")
        plt.close('all')
    else:
        plt.show()
        plt.close('all')


def path_overlap(path1, path2, lengthData, source=None, printInfo=False):
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