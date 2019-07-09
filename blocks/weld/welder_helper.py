import re
import sys
from paths import Path
import aligner
sys.path.append('../')
import log

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
        log.out("Failure in aligning start of block. Block will be skipped. " + \
                "Reference=" + rid + " Query=" + qid, 1, param)
        #todo: retry alignment?

    endFork = aligner.align_left(qid, rid, qSeq, rSeq, qend, rend, param, alignBuffer, rdir)
    if endFork is not None: endFork.switch_reference()
    else:
        log.out("Failure in aligning end of block. Block will be skipped. " + \
                "Reference=" + rid + " Query=" + qid, 1, param)        
        #todo: retry alignment?

    return (startFork, endFork)

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
        log.out("No corresponding sequence in reference was found on " + side + " side. Skipping N filling.", 1, param)
        return None

    if side == 'l':
        fork = aligner.align_left(qid, rid, qSeq, rSeq, qp, rpos, param, alignBuffer, rdir)
    else:
        fork = aligner.align_right(qid, rid, qSeq, rSeq, qp, rpos, param, alignBuffer, rdir)

    if fork is None: return None
    
    sameStrandStart = fork.before_strand() == fork.after_strand()
    if sameStrand != sameStrandStart: return None

    return fork
 
def fill_NNNs(qSeq, rSeq, block, Nstart, Nend, param):
    '''
    Identifies and creates Forks around NNNs in the query sequence between Nstart and
    Nend using the block to guide the alignments. Returns a Path that fills in NNNs.
    '''
    minPos = min(Nstart, Nend)
    maxPos = max(Nstart, Nend)
    pathN = Path()

    log.out("Finding Ns from position " + str(minPos) + " to " + str(maxPos) + ".", 2, param)
        
    for qstart, qend in find_NNNs(qSeq, minPos, maxPos):

        log.out("Ns found between position " + str(qstart) + " to " + str(qend) + ".", 2, param)        
        log.out("Attempting to fill Ns...", 2, param)

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
                    log.out("N filling unsuccessful. Giving up.", 2, param)
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
            log.out("N filling successful.", 2, param)
            break
        
    log.out("N filling complete.", 2, param)
    return pathN