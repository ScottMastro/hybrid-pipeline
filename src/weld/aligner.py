import sys
sys.path.append('..')
import utils.log as logger

import weld.alignment_cutter as cutter
from structures.path import Fork

#---- ssw aligner required for fast alignments ----
import ssw
alnr = ssw.Aligner()
#---- ssw aligner required for fast alignments ----

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

#==================================================
# Main Smith-Waterman alignment function
#==================================================

def ssw_align(qSeq, rSeq, qRange, rRange, qReverse=False, rReverse=False, printAlignment=False):
    '''
    Runs a single striped Smith-Waterman alignment between qSeq and rSeq and the
    specfied ranges (startPos, endPos). Reverse complement is used if specified.
    
    Returns an alignment object which includes the start and end position
    of alignment on qSeq and rSeq (qstart, qend, rstart, rend)
    
    https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
    '''
    qSubSeq = qSeq[qRange[0]:qRange[1]]
    rSubSeq = rSeq[rRange[0]:rRange[1]]
    if qReverse: qSubSeq = reverse_complement(qSubSeq)
    if rReverse: rSubSeq = reverse_complement(rSubSeq)

    alignment = alnr.align(qSubSeq, rSubSeq, revcomp=False)  

    if qReverse:
        alignment.qstart = len(qSeq) - (len(qSubSeq) - alignment.query_begin + qRange[0])
        alignment.qend = len(qSeq) - (len(qSubSeq) - alignment.query_end + qRange[0])
        alignment.qstrand = -1

    else:
        alignment.qstart = alignment.query_begin + qRange[0]
        alignment.qend = alignment.query_end + qRange[0]
        alignment.qstrand = 1
        
    if rReverse:        
        alignment.rstart = len(rSeq) - (len(rSubSeq) - alignment.reference_begin + rRange[0])
        alignment.rend = len(rSeq) - (len(rSubSeq) - alignment.reference_end + rRange[0])
        alignment.rstrand = -1
    else:
        alignment.rstart = alignment.reference_begin + rRange[0]
        alignment.rend = alignment.reference_end + rRange[0]
        alignment.rstrand = 1
    
    alignment.qend = alignment.qend + 1
    alignment.rend = alignment.rend + 1

    if printAlignment:
        print(alignment.alignment_report())

    return alignment

#==================================================
# Flank aligning functions
#==================================================

def get_positions(seq, pos, alignBuffer, extension, side):
    '''
    Gets a left/right pair of positions to attempt alignment on.
    '''
    leftBuffer = 0 if side == 'r' else alignBuffer + extension
    rightBuffer = 0 if side == 'l' else alignBuffer + extension
    
    left= max(pos - leftBuffer, 0) 
    right = min(pos + rightBuffer, len(seq))
    if left == 0:
        right = min(rightBuffer + leftBuffer, len(seq))
    if right == len(seq):
        left = max(len(seq) - rightBuffer - leftBuffer, 0) 

    return (left, right)

def align(qSeq, rSeq, qpos, rpos, side, param, alignBuffer=200, hint=None):
    '''
    Aligns qSeq and rSeq. Optionally pass reference strand hint as a parameter. 
    Attempts to cut alignment at high-confidence position. Returns tuple of
    (qCutPos, qStrand, rCutPos, rStrand) or None if alignment fails.
    '''
    verboseAlignment = False
    
    extension = 0
    retryCount = 0
    sideStr="right" if side == 'r' else "left"
    
    while extension < 20000:
        
        logger.log("Aligning " + sideStr + " side, attempt=" + str(retryCount+1), logger.LOG_DEBUG, indent=2)

        qRange = get_positions(qSeq, qpos, alignBuffer, extension, side)
        rRange = get_positions(rSeq, rpos, alignBuffer, extension, side)
        rRangeRev = get_positions(rSeq, rpos, alignBuffer, extension, \
                           "l" if side == 'r' else 'r')        

        if hint == -1:
            alignment = ssw_align(qSeq, rSeq, qRange, rRangeRev, \
                            qReverse=False, rReverse=True, printAlignment=verboseAlignment)
        else:
            alignment = ssw_align(qSeq, rSeq, qRange, rRange, \
                            qReverse=False, rReverse=False, printAlignment=verboseAlignment)
            
        cut = cutter.cut_alignment(qSeq, rSeq, alignment, side, param)
            
        if cut is None:
            if hint == -1:
                alignment = ssw_align(qSeq, rSeq, qRange, rRange, \
                                qReverse=False, rReverse=False, printAlignment=verboseAlignment)
            else:
                alignment = ssw_align(qSeq, rSeq, qRange, rRangeRev, \
                                qReverse=False, rReverse=True, printAlignment=verboseAlignment)

            cut = cutter.cut_alignment(qSeq, rSeq, alignment, side, param)
                        
        if cut is not None:
            logger.log("Alignment successful. Cut point = q:" + str(cut[0]) + " / r:" + str(cut[1]), logger.LOG_DEBUG, indent=2)
            return (cut[0], alignment.qstrand, cut[1], alignment.rstrand)

        retryCount = retryCount + 1
        extension=max(alignBuffer, extension*2)
        logger.log("Alignment failed. Cut site not identified. Trying again with extension = " + str(extension) + "bp", logger.LOG_DEBUG, indent=2)
    
    logger.log("Giving up on alignment, " + sideStr + " side failed to align).", logger.LOG_DEBUG, indent=2)
    return None
    
def align_left(qid, rid, qSeq, rSeq, qpos, rpos, param, alignBuffer=200, hint=None):
    '''
    Aligns qSeq and rSeq on the left side of qpos/rpos. Optionally pass reference
    strand hint as a parameter. Returns a Fork is alignment succesful, otherwise None.
    '''
    result = align(qSeq, rSeq, qpos, rpos, 'l', param, alignBuffer, hint)  

    if result is None: return None
    (qp, qdir, rp, rdir) = result
    return Fork(qid, qp, qdir, rid, rp, rdir)

def align_right(qid, rid, qSeq, rSeq, qpos, rpos, param, alignBuffer=200, hint=None):
    '''
    Aligns qSeq and rSeq on the right side of qpos/rpos. Optionally pass reference
    strand hint as a parameter. Returns a Fork is alignment succesful, otherwise None.
    '''
    result = align(qSeq, rSeq, qpos, rpos, 'r', param, alignBuffer, hint)  
    
    if result is None: return None
    (qp, qdir, rp, rdir) = result
    return Fork(qid, qp, qdir, rid, rp, rdir)
