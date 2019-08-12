import welder_helper as helper
import copy

import sys
sys.path.append('../')
import log
from paths import Path
import path_helper

def weld_megablock(megablock, seqData, param):
    '''
    Welds together a megablock into a Path. Blocks in megablock have their ends
    anchored between query and reference. NNNs in query sequence are replaced 
    using reference sequence. A list of Paths is returned, one for each block in
    megablock. Each Path in the list should be internally consistent.
    '''
    megaPath = []
    
    rid, qid = megablock.rid, megablock.qid
    rSeq, qSeq = seqData[str(rid)], seqData[str(qid)]
    
    log.out("Merging contigs " + str(rid) + " and " + str(qid) + ".", 1, param)
    log.out("----------------------------------------------------", 1, param)
            
    for block in megablock:
        
        # anchor ends
        startFork, endFork = helper.anchor_block_ends(qSeq, rSeq, block, param)

        if startFork is not None and endFork is not None and \
            path_helper.valid_fork_pair(startFork, endFork):                        
                log.out("Successfully anchored block.", 2, param)
        else:        
            log.out("Could not anchor block, skipping.", 1, param)
            continue

        # find and fill NNNs
        Nstart, Nend = block.left(q=True), block.right(q=True)
        blockPath = helper.fill_NNNs(qSeq, rSeq, block, Nstart, Nend, param)     
        # combine overlapping NNN filled regions
        blockPath = path_helper.clean_overlapping_forks(blockPath, param)
        
        # merge block ends with NNN path
        if len(blockPath) < 1:
            blockPath.add_fork(startFork)
            blockPath.add_fork(endFork)
        else:
            if (startFork.after_pos() > blockPath[0].before_pos() or \
                startFork.before_pos() > blockPath[0].after_pos()):
                #remove startFork and first path fork            
                blockPath.pop(0)
            else:
                blockPath.add_fork_front(startFork)
                
            if (blockPath[-1].after_pos() > endFork.before_pos() or \
                blockPath[-1].before_pos() > endFork.after_pos()):
                #remove endFork and last path fork            
                blockPath.pop()
            else:
                blockPath.add_fork(endFork)
        
        megaPath.append(blockPath)
        
    return megaPath

def join_blockpaths(paths, lengthData, param):
    
    if len(paths) < 1:
        return Path()
    
    # debug: verify input is ok
    for path in paths:
        ok = path_helper.check_path(path)
        if not ok:
            print("Not ok")
            #input()  
    
    # fix strand orientation
    for i in range(len(paths)-1):
        paths[i], paths[i+1] = path_helper.fix_path_orientation( \
             paths[i], paths[i+1], lengthData, param, False)
    
    # join paths
    path = Path()
    for p in paths: path.add_path(p)

    # add Nforks if necessary
    path = path_helper.add_Nforks(path, lengthData)   

    return path
        
def join_megablockpaths(paths, lengthData, param):
    
    # exclude paths of length <= 2 if trimEnds = True
    filterEmpty = []
    for path in paths:
        if len(path) > 2:
            filterEmpty.append(path)
    paths = filterEmpty

    if len(paths) < 1:
        return Path()
    
    # debug: verify input is ok
    for path in paths:
        ok = path_helper.check_path(path)
        if not ok:
            print("Not ok")
            #input()  
    
    # fix strand orientation
    for i in range(len(paths)-1):
        paths[i], paths[i+1] = path_helper.fix_path_orientation( \
             paths[i], paths[i+1], lengthData, param, True)
    
    # join paths
    path = Path()
    firstFork = paths[0][0]
    lastFork = paths[-1][-1]

    if not firstFork.is_Nfork(): path.add_fork(firstFork)

    for p in paths: 
        i,j = None, None
        if not p[-1].is_Nfork(): j=-1
        if not p[0].is_Nfork(): i=1
            
        path.extend(p[i:j])
    if not lastFork.is_Nfork(): path.add_fork(lastFork)

    # add Nforks if necessary
    path = path_helper.add_Nforks(path, lengthData)   

    return path


def weld(contig, seqData, lengthData, param):
    '''
    ...
    '''
    contigPaths = []
    
    plot = False
    plotPaths = []
        
    for megablock in reversed(contig.mblocks):
        megaPath = weld_megablock(megablock, seqData, param)
        
        if plot:
            plotPaths.append(copy.deepcopy(megaPath))
        
        contigPath = join_blockpaths(megaPath, lengthData, param)
        
        if contigPath is not None and len(contigPath) > 0:
            contigPaths.append(contigPath)
        
     
    path = join_megablockpaths(contigPaths, lengthData, param)
        
    ok = path_helper.check_path(path)
    
    if not ok:
        print("Not ok")
        input()
        
    return path










        
        




        

    
    
