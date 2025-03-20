import re
import sys
sys.path.append('..')

import utils.log as logger
from structures.path import Path
from structures.path import Fork
import weld.aligner as aligner

OUTPUT_FILE_NAME = "gap_fill.bed"

#==================================================
# Fill gaps denoted by NNNs
#==================================================

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

    #no corresponding sequence found
    if abs(qp - qpos) > 20000:     
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
        
    for qstart, qend in find_NNNs(qSeq, minPos, maxPos):

        logger.log("Attempting to fill gap, " + block.qid + ":" + str(qstart) + "-" + str(qend), \
                   logger.LOG_DETAILS, indent=1)        

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
                    logger.log("Gap filling unsuccessful. Skipping.", logger.LOG_DETAILS, indent=1)
        
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
            logger.log("Gap filling successful.", logger.LOG_DETAILS, indent=1)
    
            rgb =  logger.rbg_generator(block.rid)
            rInfo = str(block.rid) + ":" + str(startFork.rpos) + "-" + str(endFork.rpos)
            qInfo = [block.qid, startFork.qpos, endFork.qpos, \
                     rInfo, 0, startFork.rstrand, \
                    startFork.qpos, endFork.qpos, rgb]
            logger.FileLogger().write_cols(OUTPUT_FILE_NAME, qInfo)

            break
        
    return pathN

#==================================================
# Replace region with low similarity
#==================================================

def fill_diffs(qSeq, rSeq, block, param):
    '''
    Identifies smaller differences between query and reference and creates Forks 
    around the query sequence between using the chunks to guide the alignments.
    Returns a Path that fills in diffs.
    '''
    #TODO: output log info
    
    pcidThresh = 94
    anchorSize = 100
    gapDetectSize = 5
    buffer = 60
    printAlignments=False
    qid = block.qid
    rid = block.rid
    pathDiff = Path()

    for i,chunk in enumerate(block[1:]):
        i = i + 1
        
        if chunk.id -1 != block[i-1].id:
            chunk1, chunk2 = block[i-1], chunk
            
            qRange = (chunk1.right(q=True) - anchorSize,  chunk2.left(q=True) + anchorSize)
            
            r1 = min( max(chunk1.left(q=False), chunk1.right(q=False)), \
                      max(chunk2.left(q=False), chunk2.right(q=False)) )
            r2 = max( min(chunk1.left(q=False), chunk1.right(q=False)), \
                      min(chunk2.left(q=False), chunk2.right(q=False)) )
            if r1 >= r2: 
                r1 = min( min(chunk1.left(q=False), chunk1.right(q=False)), \
                          min(chunk2.left(q=False), chunk2.right(q=False)) )
                r2 = max( max(chunk1.left(q=False), chunk1.right(q=False)), \
                          max(chunk2.left(q=False), chunk2.right(q=False)) )
            rRange = (r1 - anchorSize, r2 + anchorSize)
            rReverse = (chunk1.get_dir(q=False) == -1)
            
            alignment = aligner.ssw_align(qSeq, rSeq, qRange, rRange, False, rReverse, printAlignments)
            
            x1 = alignment.alignment[1][:anchorSize].count("|")
            x2 = alignment.alignment[1][-anchorSize:].count("|")
            
            if 1.0*x1 /anchorSize > pcidThresh or 1.0*x2 /anchorSize > pcidThresh:
                continue   
            
            start, end = None, None
            count = 0
            for l,cigar in alignment.iter_cigar:
                if (cigar == 'I' or cigar == 'D') and l >= gapDetectSize:
                    if start is None:
                        start = count            
                    end = count + l
                count += l
            if start is None or end is None: continue
            leftCutResult, rightCutResult = None, None  
            while start >= buffer:
                start = start - buffer
                leftCutResult = aligner.cut_alignment(qSeq, rSeq, alignment, 'l', param, oneTry=start)
                if leftCutResult is not None: break
            
            if leftCutResult is None: continue

            length = len(alignment.alignment[1])
            while end <= length-buffer:
                end = end + buffer
                rightCutResult = aligner.cut_alignment(qSeq, rSeq, alignment, 'r', param, oneTry=end)
                if rightCutResult is not None: break
                            
            if rightCutResult is None: continue
  
            lq, lr = leftCutResult
            rq, rr = rightCutResult
    
            leftFork =  Fork(qid, lq, 1, rid, lr, -1 if rReverse else 1)
            rightFork = Fork(qid, rq, 1, rid, rr, -1 if rReverse else 1)
            leftFork.switch_reference()
            rightFork.switch_query()

            pathDiff.add_fork(leftFork)
            pathDiff.add_fork(rightFork)
            
    return pathDiff
