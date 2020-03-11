import sys
sys.path.append("..")
import log.log as logger

from weld.paths import Path
import weld.path_helper as path_helper
import weld.aligner as aligner
import weld.N_filler as nfill
import weld.diff_filler as dfill

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
    else: 
        logger.out("Failure in aligning start of block. Block will be skipped. " + \
                "Reference=" + rid + " Query=" + qid, 1, param)
        #todo: retry alignment?

    endFork = aligner.align_left(qid, rid, qSeq, rSeq, qend, rend, param, alignBuffer, rdir)
    if endFork is not None: endFork.switch_reference()
    else:
        logger.out("Failure in aligning end of block. Block will be skipped. " + \
                "Reference=" + rid + " Query=" + qid, 1, param)        
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
    
    logger.out("Merging contigs " + str(rid) + " and " + str(qid) + ".", 1, param)
    logger.out("----------------------------------------------------", 1, param)
            
    for block in megablock:
        
        # anchor ends
        startFork, endFork = anchor_block_ends(qSeq, rSeq, block, param)

        if startFork is not None and endFork is not None and \
            path_helper.valid_fork_pair(startFork, endFork):                        
                logger.out("Successfully anchored block.", 2, param)
        else:        
            logger.out("Could not anchor block, skipping.", 1, param)
            continue

        # find and fill NNNs
        blockPath = nfill.fill_NNNs(qSeq, rSeq, block, param)
        if replaceChunks:
            diffPath = dfill.fill_diffs(qSeq, rSeq, block, param)
            blockPath = path_helper.interleave_paths(blockPath, diffPath)
            
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
    '''
    Joins together a list of Paths constructed from Blocks.
    Returns a single Path.
    '''

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