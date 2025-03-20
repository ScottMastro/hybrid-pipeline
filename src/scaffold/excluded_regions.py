import copy
import sys
sys.path.append("..")

import utils.log as logger
import structures.region as rgn
from structures.region import region_difference
from structures.path import Fork

def get_unused_regions(paths, tigIds, lengthData, minSize=10000):

    fullRegions = [ rgn.SimpleRegion(tigId, 0, lengthData[tigId]-1) for tigId in tigIds ]
    usedRegions = { tigId : [] for tigId in tigIds }

    def get_region(prevFork, fork, tigId):
        start = prevFork.get_pos_by_id_norm(tigId, lengthData)
        end = fork.get_pos_by_id_norm(tigId, lengthData)
        if start is None or end is None: return None
        return rgn.SimpleRegion(tigId, min(start, end), max(start, end))

    for path in paths:
        if len(path) < 2: continue
    
        for i,fork in enumerate(path[1:]):
            prevFork = path[i]

            if not fork.is_Nfork() and not prevFork.is_Nfork():
                r1 = get_region(prevFork, fork, prevFork.after_id())
                r2 = get_region(prevFork, fork, prevFork.before_id())

                if r1 is not None and r1.chrom in usedRegions:
                    usedRegions[r1.chrom].append(r1)
                    
                if r2 is not None and r2.chrom in usedRegions:
                    if (len(r2) - len(r1) > minSize):
                        continue
                    else:
                        usedRegions[r2.chrom].append(r2)

    unusedRegions = { tigId : [] for tigId in tigIds }
    for region in fullRegions:
        unusedRegions[region.chrom] = region_difference(region, usedRegions[region.chrom], minSize)
    
    return unusedRegions


def get_extension(fork, unusedDict, lengthData, maxDist, start=True):
    
    def get_extension_region(before=True):
        if before:
            strand = fork.before_strand()
            pos = fork.before_pos_norm(lengthData)
            tigId = fork.before_id()
        else:
            strand = fork.after_strand()
            pos = fork.after_pos_norm(lengthData)
            tigId = fork.after_id()
        
        if tigId is None: return None
        d = strand if start else -strand
        for region in unusedDict[tigId]:
            
            condition1 = region.start < pos if d == 1 else \
                         region.end > pos
            condition2 = abs(region.end - pos) < maxDist if d == 1 else \
                         abs(region.start - pos) < maxDist

            if condition1 and condition2:
                return region
        return None
    
    region1 = get_extension_region(before=True)
    region2 = get_extension_region(before=False)
    
    region = region1
    if region1 is None or (region2 is not None and len(region1) < len(region2)):
        region = region2
    return region
        

def get_fork_cap(region, fork, lengthData, start=True):
    q = (fork.qid == region.chrom)
    strand = fork.get_strand(q)
    d = strand if start else -strand
    
    pos = region.start if d==1 else region.end
    pos = pos if strand == 1 else lengthData[fork.get_id(q)]-pos
    
    if q:  cap = Fork(region.chrom, pos, strand, None, None, None)
    else:  cap = Fork(None, None, None, region.chrom, pos, strand)
        
    if q == start:  cap.switch_query()
    else:           cap.switch_reference()
    
    #a buffer is needed if the cap does not
    bufferFork = None
    if start and cap.after_id() != fork.before_id():
        if cap.is_switch_query():
            bufferFork = Fork(cap.qid, fork.qpos-fork.qstrand, cap.qstrand, fork.rid, fork.rpos, fork.rstrand)
            bufferFork.switch_reference()
        if cap.is_switch_reference():
            bufferFork = Fork(fork.qid, fork.qpos, fork.qstrand, cap.rid, fork.rpos-fork.rstrand, cap.rstrand,)
            bufferFork.switch_query()
    if not start and fork.after_id() != cap.before_id():
        if fork.is_switch_query():
            bufferFork = Fork(fork.qid, fork.qpos, fork.qstrand, cap.rid, fork.rpos-fork.rstrand, cap.rstrand,)
            bufferFork.switch_reference()
        if cap.is_switch_reference():
            bufferFork = Fork(cap.qid, fork.qpos-fork.qstrand, cap.qstrand, fork.rid, fork.rpos, fork.rstrand)
            bufferFork.switch_query()

    return (cap, bufferFork)

def extend_path_ends(paths, unusedDict, lengthData, maxDist=10):
    
    ud = copy.deepcopy(unusedDict)
         
    for path in paths:
        logger.log("Attempting to extend path: (" + path.pid + ").", logger.LOG_DEBUG)

        startFork = path[0]
        logger.log("Starting fork: (" + str(startFork) + ").", logger.LOG_DEBUG, indent=1)

        region = get_extension(startFork, ud, lengthData, maxDist, start=True)
        if region is not None:
            cap, bufferFork = get_fork_cap(region, startFork, lengthData, start=True)
            ud[region.chrom].remove(region)
            
            if bufferFork is not None:
                path.add_fork_front(bufferFork)
            path.add_fork_front(cap)
            
        logger.log("New starting fork: (" + str(path[0]) + ").", logger.LOG_DEBUG, indent=1)

        endFork = path[-1]
        logger.log("Ending fork: (" + str(endFork) + ").", logger.LOG_DEBUG, indent=1)

        region = get_extension(endFork, ud, lengthData, maxDist, start=False)
        if region is not None:
            cap, bufferFork = get_fork_cap(region, endFork, lengthData, start=False)
            ud[region.chrom].remove(region)
            
            if bufferFork is not None:
                path.add_fork(bufferFork)
            path.add_fork(cap)
            
        logger.log("New ending fork: (" + str(path[-1]) + ").", logger.LOG_DEBUG, indent=1)

            
    return paths, ud
