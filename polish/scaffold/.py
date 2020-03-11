import sys
sys.path.append("..")
import log.log as logger

import weld.path_helper as path_helper

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

