import sys
import ssw_aligner
from reversecomp import reverse_complement

sys.path.append('../')
import log
from path_helper import Fork

def align(seqA, seqB, posA, posB, alignBuffer, extension, side, printAlignment=False):
    #todo: speedup? pass direction as hint?    
    '''
    seqA, seqB = qSeq, rSeq,
    posA, posB = qstart, rstart
    posA, posB = qend, rend

    alignBuffer = 20000
    printAlignment=True
    '''

    leftBuffer = 0 if side == 'r' else alignBuffer + extension
    rightBuffer = 0 if side == 'l' else alignBuffer + extension
    
    leftA = max(posA - leftBuffer, 0) 
    rightA = min(posA + rightBuffer, len(seqA))
    if leftA == 0:
        rightA = min(rightBuffer + leftBuffer, len(seqA))
    if rightA == len(seqA):
        leftA = max(len(seqA) - rightBuffer - leftBuffer, 0) 
    subSeqA = seqA[leftA:rightA]
    
    leftB = max(posB - leftBuffer, 0)
    rightB = min(posB + rightBuffer, len(seqB))
    if leftB == 0:
        rightB = min(rightBuffer + leftBuffer, len(seqB))
    if rightB == len(seqB):
        leftB = max(len(seqB) - rightBuffer - leftBuffer, 0) 
    subSeqB = seqB[leftB:rightB]

    results = ssw_aligner.align(subSeqB, subSeqA, "B", "A", forward=True, reverse=False)    
    
    leftBReverse = max(posB - rightBuffer, 0)
    rightBReverse = min(posB + leftBuffer, len(seqB))    
    if leftBReverse == 0:
        rightBReverse = min(rightBuffer + leftBuffer, len(seqB))
    if rightBReverse == len(seqB):
        leftBReverse = max(len(seqB) - rightBuffer - leftBuffer, 0) 
    subSeqBReverse = seqB[leftBReverse:rightBReverse]

    reverseResults= ssw_aligner.align(subSeqBReverse, subSeqA, "B", "A", forward=False, reverse=True)    

    if results[-1] > reverseResults[-1]:
        if printAlignment:
            ssw_aligner.align(subSeqB, subSeqA, "B", "A", forward=True, reverse=False, printResults=True)    

        alnB, alnStartB, alnEndB, aln, \
        alnA, alnStartA, alnEndA, score = results
        direction = 1
    else:
        if printAlignment:
            ssw_aligner.align(subSeqBReverse, subSeqA, "B", "A", forward=False, reverse=True, printResults=True)    
        
        alnB, alnStartB, alnEndB, aln, \
        alnA, alnStartA, alnEndA, score = reverseResults
        direction = -1
        (leftB, rightB) = (leftBReverse, rightBReverse)

    if(side == 'l'):
        A = leftA + alnStartA
        
        if direction == 1: B = leftB + alnStartB
        if direction == -1: B = len(seqB) - (rightB - alnStartB)
            
        if(aln[:100].count('|') < 92):
            return None
        
    if(side == 'r'):
        A = rightA - (len(subSeqA) - alnEndA)
    
        if direction == 1: B = rightB - (len(subSeqB) - alnEndB)
        if direction == -1: B = len(seqB) - (leftB + len(subSeqB) - alnEndB)
        
        if(aln[-100:].count('|') < 92 and side == 'r'):
            return None

    splitSize=5
    splitShift=[0, 3, 5, 8, 10, 13, 15]

    for shift in splitShift:
        
        if side == 'r': 
            s = -shift+splitSize
            shiftA = -shift + (0 if s >= 0 else alnA[s:].count('-'))
            shiftB = (-shift + (0 if s >= 0 else alnB[s:].count('-')))*direction

        if side == 'l': 
            s = shift-splitSize
            shiftA = shift - (0 if s <= 0 else alnA[:s].count('-'))
            shiftB = (shift - (0 if s <= 0 else alnB[:s].count('-')))*direction

        splitA = seqA[max(0, A-splitSize+shiftA):A+shiftA] + "-" + seqA[A+shiftA:A+splitSize+shiftA]
        B_ = B
        if direction == -1:
            B_ = len(seqB) - B_
    
        splitB = seqB[max(0, B_-splitSize+shiftB):B_+shiftB] + "-" + seqB[B_+shiftB:B_+splitSize+shiftB]
        if direction == -1:
            splitB = reverse_complement(splitB)
    
        if printAlignment:
            print(splitA)
            print(splitB)
            print(splitA == splitB)
        
        if splitA == splitB:
            break
        
    if splitA != splitB:
        return None
    
    return (A+shiftA, 1, B+(shiftB*direction), direction)

def align_side(rSeq, qSeq, rpos, qpos, side, param, alignBuffer=500):
    extension = 0
    retryCount = 0
    sideStr="right" if side == 'r' else "left"
    
    while extension < 20000:
        cut = align(qSeq, rSeq, qpos, rpos, alignBuffer, extension, side, param.VERBOSE >= 3)
        if cut is not None:
            return cut
        else:
            retryCount = retryCount + 1
            extension=max(alignBuffer, extension*2)
            log.out("Trying " + sideStr + " again with extension = " + str(extension), 3, param, wait=True)

    log.out("Could not create bubble (" + sideStr + " failure).", 1, param)
    return None

def align_left(rid, qid, rSeq, qSeq, rpos, qpos, param, alignBuffer=500):
    result = align_side(rSeq, qSeq, rpos, qpos, 'l', param, alignBuffer)  

    if result is None: return None
    (qp, qdir, rp, rdir) = result
    return Fork(qid, qp, qdir, rid, rp, rdir)

def align_right(rid, qid, rSeq, qSeq, rpos, qpos, param, alignBuffer=500):
    result = align_side(rSeq, qSeq, rpos, qpos, 'r', param, alignBuffer)  
    
    if result is None: return None
    (qp, qdir, rp, rdir) = result
    return Fork(qid, qp, qdir, rid, rp, rdir)
