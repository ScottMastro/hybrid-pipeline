import sys
sys.path.append('..')
import log as logger
import re
from structures.path import Path
import weld.aligner as aligner

OUTPUT_FILE_NAME = "gap_fill.bed"

def find_NNNs(sequence, minIdx=-1, maxIdx=-1, buffer=1500):
    """
    Generator function for iterating over a sequence and finding consecutive NNNs
    between minIdx and maxIdx. Groups of NNNs separated by less than buffer 
    are grouped together.
    """

    prevStart, prevEnd = -1, -1
    
    p = re.compile("N+")
    for match in p.finditer( \
            sequence[None if minIdx < 0 else minIdx: None if maxIdx < 0 else maxIdx]):
        start = match.start() + (0 if minIdx < 0 else minIdx)
        end = match.end() + (0 if minIdx < 0 else minIdx)
        
        if prevEnd > 0:
            
            #check if distance between NNNs is less than buffer
            if buffer > 0 and abs(start - prevEnd) <= buffer:
                prevEnd = end
                continue
            
            yield(prevStart, prevEnd)
       
        prevStart, prevEnd = start, end

    #returns last set of NNNs
    if prevEnd > 0:
        yield(prevStart, prevEnd)
        
def fill_NNNs_side(qSeq, rSeq, block, qpos, side, extension, param):
    """
    Aligns one side of the sequence flanking NNNs. Returns fork if alignment successful.
    Returns None if alignment fails.
    """
    qid = block.qid
    rid = block.rid

    sameStrand = (block.get_dir(q=True) == block.get_dir(q=False))

    alignBuffer = param.ALIGN_BUFFER
    qpos = qpos + extension*(-1 if side == 'l' else 1)
    rdir = block.get_dir(q=False)

    rpos = block.closest_corresponding_position(qpos, q=True, side=side)
    qp = block.closest_corresponding_position(rpos, q=False)

    if abs(qp - qpos) > 20000:     
        print(extension)
        logger.out("No corresponding sequence in reference was found on " + side + " side. Skipping N filling.", 1, param)
        return None

    if side == 'l':
        fork = aligner.align_left(qid, rid, qSeq, rSeq, qp, rpos, param, alignBuffer, rdir)
    else:
        fork = aligner.align_right(qid, rid, qSeq, rSeq, qp, rpos, param, alignBuffer, rdir)

    if fork is None: return None
    
    sameStrandStart = fork.before_strand() == fork.after_strand()
    if sameStrand != sameStrandStart: return None

    return fork
 
def fill_NNNs(qSeq, rSeq, block, param):
    '''
    Identifies and creates Forks around NNNs in the query sequence between Nstart and
    Nend using the block to guide the alignments. Returns a Path that fills in NNNs.
    '''
    
    Nstart, Nend = block.left(q=True), block.right(q=True)

    minPos = min(Nstart, Nend)
    maxPos = max(Nstart, Nend)
    pathN = Path()

    logger.out("Finding Ns from position " + str(minPos) + " to " + str(maxPos) + ".", 2, param)
        
    for qstart, qend in find_NNNs(qSeq, minPos, maxPos):

        logger.out("Ns found between position " + str(qstart) + " to " + str(qend) + ".", 2, param)        
        logger.out("Attempting to fill Ns...", 2, param)

        startFork, endFork = None, None
        retries = 3
        extensionIncrease = 1500
        extension = 0
        
        while True:
            startFork = fill_NNNs_side(qSeq, rSeq, block, qstart, 'l', extension, param)
            endFork = fill_NNNs_side(qSeq, rSeq, block, qend, 'r', extension, param)

            if startFork is None or endFork is None:
                if retries > 0:
                    retries = retries - 1
                    extension = extension + extensionIncrease
                    continue
                else:
                    logger.out("N filling unsuccessful. Giving up.", 2, param)
        
                    rgb =  logger.rbg_generator(None)
                    info = [block.qid, qstart, qend, ".", 0, "+", qstart, qend, rgb]
                    logger.FileLogger().write_cols(OUTPUT_FILE_NAME, info)

                    break
                
            startFork.switch_reference()
            endFork.switch_query()

            #check for weird alignment due to repeated segments
            if startFork.after_pos() > endFork.before_pos() or \
                startFork.before_pos() > endFork.after_pos():
                extension = extension + max(200, abs(startFork.after_pos() - endFork.before_pos()))
                continue

            pathN.add_fork(startFork)
            pathN.add_fork(endFork)
            logger.out("N filling successful.", 2, param)
    
            rgb =  logger.rbg_generator(block.rid)
            rInfo = str(block.rid) + ":" + str(startFork.rpos) + "-" + str(endFork.rpos)
            qInfo = [block.qid, startFork.qpos, endFork.qpos, \
                     rInfo, 0, startFork.rstrand, \
                    startFork.qpos, endFork.qpos, rgb]
            logger.FileLogger().write_cols(OUTPUT_FILE_NAME, qInfo)


            break
        
    logger.out("N filling complete.", 2, param)
    return pathN