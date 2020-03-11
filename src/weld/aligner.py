import sys
sys.path.append('../')
import log as logger

import weld.aligner_helper as helper
from structures.fork import Fork

#---- ssw aligner required for fast alignments ----
import ssw
alnr = ssw.Aligner()
#---- ssw aligner required for fast alignments ----

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

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

def cut_alignment(qSeq, rSeq, alignment, side, param, buffer=50, minPercent=90, smallBuffer=5, oneTry=None):
    '''
    Attempts to cut alignment alignment at high-confidence position. Checks
    percent identity of alignment and entropy around buffer, check that
    alignment is exact within small buffer. Returns tuple of (qCutPos, rCutPos)
    or None if no cut position could be found.
    '''
    if len(alignment.alignment[1]) < buffer:
        return None

    terminate = False
    attempt = 0
    
    while True:
        
        if terminate: return None
            
        #grab position to attempt to cut the alignment at
        if oneTry is not None:
            cutPos, terminate = oneTry, True
        else:
            cutPos, terminate = helper.get_cut_position(alignment, attempt, side, buffer)
        attempt += 1
        logger.out("Attempting to find cut site at position: " + str(cutPos), 5, param)


        #check the quality of the alignment in this region
        idPass = helper.check_identity(alignment, cutPos, side, buffer, minPercent)
        if not idPass: 
            logger.out("Percent identity around cut region too low.", 5, param)
            continue

        cutPos = helper.ajdust_cut_pos(alignment, cutPos, side, buffer, smallBuffer)
        if cutPos is None: 
            logger.out("Could not find high-entropy exact match for cut site.", 5, param) 
            continue
        
    
        qCut = alignment.qstart + cutPos - alignment.alignment[2][:cutPos].count('-')
        rCut = alignment.rstart + cutPos - alignment.alignment[0][:cutPos].count('-') 
        
        if param.VERBOSE >= 5:
            s = 15
            qs = reverse_complement(qSeq[len(qSeq)-qCut:len(qSeq)-(qCut-s)]) + "|" + \
                      reverse_complement(qSeq[len(qSeq)-qCut:len(qSeq+s)-qCut])\
                if alignment.qstrand == -1 else \
                qSeq[qCut-s:qCut] + "|" + qSeq[qCut:qCut+s]
        
            rs = reverse_complement(rSeq[len(rSeq)-rCut:len(rSeq)-(rCut-s)]) + "|" + \
                      reverse_complement(rSeq[len(rSeq)-(rCut+s):len(rSeq)-rCut]) \
                if alignment.rstrand == -1 else \
                rSeq[rCut-s:rCut] + "|" + rSeq[rCut:rCut+s]
        
            
            print("Cut result:")
            print(qs + " query")
            print(rs + " reference")
        
        return (qCut, rCut)

    return None

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

    verbose = param.VERBOSE >= 4
    extension = 0
    retryCount = 0
    sideStr="right" if side == 'r' else "left"
            
    while extension < 20000:
        qRange = get_positions(qSeq, qpos, alignBuffer, extension, side)
        rRange = get_positions(rSeq, rpos, alignBuffer, extension, side)
        rRangeRev = get_positions(rSeq, rpos, alignBuffer, extension, \
                           "l" if side == 'r' else 'r')        

        '''
        print(qRange)
        print(qSeq[qRange[0]:qRange[1]])
        print(rRange)
        print(qSeq[rRange[0]:rRange[1]])
        print(rRangeRev)
        print(qSeq[rRangeRev[0]:rRangeRev[1]])
        '''

        if hint == -1:
            alignment = ssw_align(qSeq, rSeq, qRange, rRangeRev, \
                            qReverse=False, rReverse=True, printAlignment=verbose)
        else:
            alignment = ssw_align(qSeq, rSeq, qRange, rRange, \
                            qReverse=False, rReverse=False, printAlignment=verbose)
            
        cut = cut_alignment(qSeq, rSeq, alignment, side, param)
            
        if cut is None:
            if hint == -1:
                alignment = ssw_align(qSeq, rSeq, qRange, rRange, \
                                qReverse=False, rReverse=False, printAlignment=verbose)
            else:
                alignment = ssw_align(qSeq, rSeq, qRange, rRangeRev, \
                                qReverse=False, rReverse=True, printAlignment=verbose)

            cut = cut_alignment(qSeq, rSeq, alignment, side, param)
                        
        if cut is not None:
            logger.out("Alignment of " + sideStr + " side successful.\n" + \
                       "Cut point found, q:" + str(cut[0]) + " / r:" + str(cut[1]), 4, param, wait=True)
            return (cut[0], alignment.qstrand, cut[1], alignment.rstrand)

        retryCount = retryCount + 1
        extension=max(alignBuffer, extension*2)
        logger.out("Cut site not identified. Trying " + sideStr + " again with extension = " + str(extension), 4, param, wait=True)
        
    logger.out("Could not align (" + sideStr + " failure).", 1, param)

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


def align_chunk(chunk1, chunk2, qSeq, rSeq, param, alignBuffer=500, hint=None):
    qRange = (chunk1.left(q=True),  chunk2.right(q=True))
    rRange = (min(chunk1.left(q=False), chunk1.right(q=False),
                 chunk2.left(q=False), chunk2.right(q=False) ) ,
              max(chunk1.left(q=False), chunk1.right(q=False),
              chunk2.left(q=False), chunk2.right(q=False) ) )
        
    print(rRange, qRange)
    rReverse = (chunk1.get_dir(q=False) == -1)
    print(rReverse)
    return ssw_align(qSeq, rSeq, qRange, rRange, False, rReverse, True)

