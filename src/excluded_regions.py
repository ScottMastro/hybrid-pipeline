from structures.region import SimpleRegion
from structures.region import region_difference

from structures.fork import Fork
import copy

def parse_block_map(fileName, param):
    qreader = open(param.OUTPUT_DIR + "/q_" + fileName, "r")
    rreader = open(param.OUTPUT_DIR + "/r_" + fileName, "r")

    def parse_block(qline, rline):
        qparts = qline.split("\t")
        rparts = rline.split("\t")
        q = SimpleRegion(qparts[0], int(qparts[1]), int(qparts[2]))
        r = SimpleRegion(rparts[0], int(rparts[1]), int(rparts[2]))
        return (q,r)
    
    regionBlocks = [parse_block(qline, rline) for qline, rline in zip(qreader, rreader)]
    qreader.close() ; rreader.close()

    blockDict = dict()
    for x,y in regionBlocks:
        if x.chrom not in blockDict:
            blockDict[x.chrom] = []
        if y.chrom not in blockDict:
            blockDict[y.chrom] = []

        blockDict[x.chrom].append((x,y))
        blockDict[y.chrom].append((y,x))

    return blockDict

def get_unused_contig(contig, lengthData, minSize=None, megablockLevel=True):

    tigId = contig.id
    fullRegion = SimpleRegion(tigId, 0, lengthData[tigId]-1) 
    
    usedRegions = []
    for mblock in contig.mblocks:
        if not megablockLevel:
                for block in mblock:
                    usedRegion = SimpleRegion(tigId, block.left(q=True), block.right(q=True)) 
                    usedRegions.append(usedRegion)            
        else:
            usedRegion = SimpleRegion(tigId, mblock.left(q=True), mblock.right(q=True)) 
            usedRegions.append(usedRegion)
            
    unusedRegions = region_difference(fullRegion, usedRegions, minSize)



def get_unused_regions(paths, tigIds, lengthData, minSize=10000):
    
    fullRegions = [ SimpleRegion(tigId, 0, lengthData[tigId]-1) for tigId in tigIds ]
    usedRegions = { tigId : [] for tigId in tigIds }
    
    def normalized_pos(fork, tigId):
        pos = fork.get_pos_by_id(tigId)
        if pos is None: return None
        if fork.get_strand_by_id(tigId) == -1:
            pos = lengthData[str(tigId)] - pos
        return pos
    
    def get_region(prevFork, fork, tigId):
        start = normalized_pos(prevFork, tigId)
        end = normalized_pos(fork, tigId)
        if start is None or end is None: return None
        return SimpleRegion(tigId, min(start, end), max(start, end))

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

def valid_fork_pair(beforeFork, afterFork):
    if beforeFork.is_Nfork() or afterFork.is_Nfork():
        return True
    
    if beforeFork.after_id() != afterFork.before_id(): return False
    if beforeFork.after_strand() != afterFork.before_strand(): return False
    if beforeFork.after_pos() > afterFork.before_pos(): return False
    return True

def _get_extension(fork, unusedDict, lengthData, maxDist, start=True):
    
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
        

def _get_fork_cap(region, fork, lengthData, start=True):
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
        #print("~~~~~~")
        startFork = path[0]
        region = _get_extension(startFork, ud, lengthData, maxDist, start=True)
        if region is not None:
            cap, bufferFork = _get_fork_cap(region, startFork, lengthData, start=True)
            ud[region.chrom].remove(region)
            
            if bufferFork is not None:
                path.add_fork_front(bufferFork)
            path.add_fork_front(cap)
            '''
            flag=True
            if bufferFork is not None:
                flag = flag and (valid_fork_pair(cap, bufferFork))
                flag= flag and (valid_fork_pair(bufferFork, startFork))
            else:
                flag= flag and (valid_fork_pair(cap, startFork))
        
            if not flag:
                print(cap)
                if bufferFork is not None:
                    print(bufferFork.before_id(), bufferFork.before_strand(), bufferFork.before_pos())
                    print(bufferFork.after_id(), bufferFork.after_strand(), bufferFork.after_pos())
                else:
                    print("-")
                print(startFork.before_id(), startFork.before_strand(), startFork.before_pos())
                print(startFork.after_id(), startFork.after_strand(), startFork.after_pos())
                '''

    for path in paths:
        #print("~~~~~~")
        endFork = path[-1]
        region = _get_extension(endFork, ud, lengthData, maxDist, start=False)
        if region is not None:
            cap, bufferFork = _get_fork_cap(region, endFork, lengthData, start=False)
            ud[region.chrom].remove(region)
            
            if bufferFork is not None:
                path.add_fork(bufferFork)
            path.add_fork(cap)

            '''
            flag=True
            if bufferFork is not None:
                flag = flag and (valid_fork_pair(endFork, bufferFork))
                flag= flag and (valid_fork_pair(bufferFork, cap))
            else:
                flag= flag and (valid_fork_pair(endFork, cap))
        
            if not flag:
                print(endFork.before_id(), endFork.before_strand(), endFork.before_pos())
                print(endFork.after_id(), endFork.after_strand(), endFork.after_pos())
                
                if bufferFork is not None:
                    print(bufferFork.before_id(), bufferFork.before_strand(), bufferFork.before_pos())
                    print(bufferFork.after_id(), bufferFork.after_strand(), bufferFork.after_pos())
                else:
                    print("-")
    
                print(cap)
            '''
            
    return paths, ud
