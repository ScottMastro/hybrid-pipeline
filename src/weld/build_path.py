import sys
sys.path.append("..")

import utils.log as logger
import structures.path as pathtools

import weld.aligner as aligner
import weld.gap_filler as filler

REORIENTATION_FILE_NAME = "inversions.txt"

#==================================================
# Aligning functions
#==================================================

def anchor_block_ends(qSeq, rSeq, block, param):
    '''
    Identifies and creates Forks at beginning and end of block, using the block
    to guide the alignments. Returns a tuple of (startFork, endFork).
    If alignment fails, None will be returned in place of a valid Fork.
    '''
    qid, rid = block.qid, block.rid
    qstart, qend = block.start(q=True), block.end(q=True)
    rstart, rend = block.start(q=False), block.end(q=False)
    rdir = block.get_dir(q=False)
    alignBuffer = param.ALIGN_BUFFER

    startFork = aligner.align_right(qid, rid, qSeq, rSeq, qstart, rstart, param, alignBuffer, rdir)
    if startFork is not None:
        startFork.switch_query()
    #todo: retry alignment?

    endFork = aligner.align_left(qid, rid, qSeq, rSeq, qend, rend, param, alignBuffer, rdir)
    if endFork is not None: 
        endFork.switch_reference()
    #todo: retry alignment?

    return (startFork, endFork)

def weld_megablock(megablock, seqData, param):
    '''
    Welds together a Megablock into a Paths. Ends of blocks are anchored between
    query and reference. NNNs in q sequence are replaced using r sequence.
    A list of Paths is returned, each corresponding to a Block.
    '''
    
    megaPath = []
    replaceChunks = param.TRUST_REF_CHUNKS
    
    rid, qid = megablock.rid, megablock.qid
    rSeq, qSeq = seqData[str(rid)], seqData[str(qid)]
    
    logger.log("Welding " + str(qid) + " and " + str(rid) + " (" + str(len(megablock)) + " blocks).", logger.LOG_DETAILS)
            
    for block in megablock:
        
        logger.log("Anchoring start and end of " + str(block.iid), logger.LOG_DETAILS, indent=1)

        # anchor ends
        startFork, endFork = anchor_block_ends(qSeq, rSeq, block, param)

        if startFork is not None and endFork is not None and \
            pathtools.consistent_forks(startFork, endFork):                        
                logger.log("Successfully anchored block.", logger.LOG_DETAILS, indent=1)
        else:        
            logger.log("Could not anchor start and end of " + str(block.iid) +", skipping.", logger.LOG_PROGRESS, indent=1)
            continue

        logger.log("Filling gaps in " + str(block.iid), logger.LOG_DETAILS, indent=1)

        # find and fill NNNs
        blockPath = filler.fill_NNNs(qSeq, rSeq, block, param)
        if replaceChunks:
            diffPath = filler.fill_diffs(qSeq, rSeq, block, param)
            blockPath = pathtools.interleave_paths(blockPath, diffPath, q=True, preferRef=True)
            
        # combine overlapping NNN filled regions
        blockPath = pathtools.clean_overlapping_forks(blockPath, param, q=True)
        
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

#==================================================
# Path functions
#==================================================

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

    if not pathtools.consistent_forks(fork1, fork2):

        message = ["Problematic fork pair detected:", "------", \
                   "Left  = " + str(fork1), "Right = " + str(fork2), "------"]
        logger.log("\n".join(message), logger.LOG_DETAILS, indent=1)

        def path_qinfo(path):
            return str(path[0].qid) + ":" + str(path[0].qpos) + "-" + str(path[-1].qpos) + ":" + str(path[0].qstrand)
        def forks_rinfo(fork1, fork2):
            return str(fork1.rid) + ":" + str(fork1.qpos) + "-" + str(fork2.qpos) + ":" + str(fork1.qstrand)

        path1Info = path_qinfo(path1)
        path2Info = path_qinfo(path2)
        
        if pathtools.consistent_forks(fork1, fork2Flip):
            logger.log("Correction found, flipping right path.", logger.LOG_DETAILS, indent=1)
            path2.flip_strands(lengthData)
            info = [path2Info, "LEFT", forks_rinfo(fork1, fork2Flip)]

        elif pathtools.consistent_forks(fork1Flip, fork2):
            logger.log("Correction found, flipping left path.", logger.LOG_DETAILS, indent=1)
            path1.flip_strands(lengthData)
            info = [path1Info, "RIGHT", forks_rinfo(fork1Flip, fork2)]

        elif pathtools.consistent_forks(fork1Flip, fork2Flip):
            logger.log("Correction found, flipping both paths.", logger.LOG_DETAILS, indent=1)
            path1.flip_strands(lengthData)
            path2.flip_strands(lengthData)
            info = [path1Info, "RIGHT", forks_rinfo(fork1Flip, fork2Flip)]
            logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)
            info = [path2Info, "LEFT", forks_rinfo(fork1Flip, fork2Flip)]
        else:
            logger.log("No suitable correction found!", 2)
            info = [path1Info, "FAIL", path2Info]

        logger.FileLogger().write_cols(REORIENTATION_FILE_NAME, info)
        
    logger.FileLogger().flush(REORIENTATION_FILE_NAME)
    return(path1, path2)

def join_blockpaths(paths, lengthData, param):
    '''
    Joins together a list of Paths constructed from Blocks.
    Returns a single Path.
    '''
    if len(paths) < 1: return pathtools.Path()

    # fix strand orientation
    for i in range(len(paths)-1):
        paths[i], paths[i+1] = fix_path_orientation( \
             paths[i], paths[i+1], lengthData, param, False)
    
    # join paths
    path = pathtools.Path()
    for p in paths: path.add_path(p)
    # add Nforks if necessary
    path = pathtools.make_ends_consistent(path, lengthData, q=True)
    path = pathtools.make_consistent(path)   

    return path
        
def join_mblockpaths(paths, lengthData, param):
    '''
    Joins together a list of Paths constructed from Megablocks.
    Returns a single Path.
    '''

    # exclude paths of length <= 2 if trimEnds = True
    filterEmpty = []
    for path in paths:
        if len(path) > 2:
            filterEmpty.append(path)
    paths = filterEmpty

    if len(paths) < 1:
        return pathtools.Path()
        
    # fix strand orientation
    for i in range(len(paths)-1):
        paths[i], paths[i+1] = fix_path_orientation( \
             paths[i], paths[i+1], lengthData, param, True)
    
    # join paths
    path = pathtools.Path()
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
    path = pathtools.make_ends_consistent(path, lengthData, q=True)
    path = pathtools.make_consistent(path)   

    path = pathtools.clean_Nforks(path)
    path = pathtools.clean_Nforks_if_consistent(path, q=False)

    return path
