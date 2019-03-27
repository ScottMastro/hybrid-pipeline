import sys
sys.path.append('./ssw_aligner')
import ssw_aligner as ssw
from reversecomp import reverse_complement

from path import Fork

'''
def find_cut_right(aln, leftPos, rightPos, minThresh=0.95):
    if(rightPos - leftPos == 1):
        if aln[leftPos] == "|": 
            return leftPos
        else: 
            return rightPos
    
    mid = round((leftPos+rightPos)/2)
    leftAln = aln[leftPos:mid]
    leftPercent = leftAln.count("|")/max(1,len(leftAln))
    
    if(leftPercent > minThresh):
        return find_cut_right(aln, mid, rightPos, minThresh)
    else:
        return find_cut_right(aln, leftPos, mid, minThresh)

def cut(masterAln, slaveAln, aln, side, baseBuffer):
    buffer = baseBuffer
    while buffer > 0:
        pos = (buffer-1) if side == 'l' else (len(aln)-buffer)
        if aln[pos] == "|":
            break
        else:
            buffer = buffer - 1
        
    if side == 'l':
        ngapsMaster = masterAln[:buffer].count("-")
        ngapsSlave = slaveAln[:buffer].count("-")
    if side == 'r':
            ngapsMaster = masterAln[-buffer:].count("-")
            ngapsSlave = slaveAln[-buffer:].count("-")
            
    return (buffer - ngapsMaster, buffer - ngapsSlave)
        
def cut_left(masterAln, slaveAln, aln, baseBuffer):
    return cut(masterAln, slaveAln, aln, "l", baseBuffer)
def cut_right(masterAln, slaveAln, aln, baseBuffer):
    return cut(masterAln, slaveAln, aln, "r", baseBuffer)
    
    
    
def cut_test():
    
    seqA = "GGGGGAAAAAGCGCGCGCGCATTATTCCCCC"
    seqB = "GGGGGATTATTGCGCGCGCGCAAAAACCCCC"

    A:  GGGGGA--AAAA--GCGCGCGCGCATT---ATTCCCCC
    B:  GGGGGATT---ATTGCGCGCGCGCA--AAAA--CCCCC
    aln: ||||||-----|--|||||||||||-----|--|||||
              **                        *-*
    alnA = "GGGGGA--AAAA--GCGCGCGCGCATT---ATTCCCCC"
    alnB = "GGGGGATT---ATTGCGCGCGCGCA--AAAA--CCCCC"
    aln = "||||||-----|--|||||||||||-----|--|||||"
    
    
    alnB, alnStartB, AlnEndB, aln, alnA, alnStartA, alnEndA,= \
        ssw.align(seqB, seqA, "B", "A", \
                      reverse=False, printResults=True)

    



    
    masterAln = masterAln[::-1]
    slaveAln = slaveAln[::-1]
    aln = aln[::-1]

    s1^-1:  CCCCCTTA---TTACGCGCGCGCG--AAAA--AGGGGG
    s2^-1:  CCCCC--AAAA--ACGCGCGCGCGTTA---TTAGGGGG
    aln^-1: |||||--|-----|||||||||||--|-----||||||
                *-*                        **

    l=cut_left(masterAln, slaveAln, aln, baseBuffer)
    print("(l^-1) expected: (5, 5)" + " seen: " + str(l))

    r=cut_right(masterAln, slaveAln, aln, baseBuffer)
    print("(r^-1) expected: (6, 6)" + " seen: " + str(r))
'''

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
    rightA = min(posA + rightBuffer, len(seqA)-1)
    subSeqA = seqA[leftA:rightA]
    
    #print("query sequence:")
    #print(subSeqA)
    #input()
    
    leftB = max(posB - leftBuffer, 0)
    rightB = min(posB + rightBuffer, len(seqB)-1)    
    subSeqB = seqB[leftB:rightB]
    
    results = ssw.align(subSeqB, subSeqA, "B", "A", forward=True, reverse=False)    
    
    leftBReverse = max(posB - rightBuffer, 0)
    rightBReverse = min(posB + leftBuffer, len(seqB)-1)    
    subSeqBReverse = seqB[leftBReverse:rightBReverse]
    
    reverseResults= ssw.align(subSeqBReverse, subSeqA, "B", "A", forward=False, reverse=True)    
            
    if results[-1] > reverseResults[-1]:
        if printAlignment:
            ssw.align(subSeqB, subSeqA, "B", "A", forward=True, reverse=False, printResults=True)    

        alnB, alnStartB, alnEndB, aln, \
        alnA, alnStartA, alnEndA, score = results
        direction = 1
    else:
        if printAlignment:
            ssw.align(subSeqBReverse, subSeqA, "B", "A", forward=False, reverse=True, printResults=True)    
        
        alnB, alnStartB, alnEndB, aln, \
        alnA, alnStartA, alnEndA, score = reverseResults
        direction = -1
        (leftB, rightB) = (leftBReverse, rightBReverse)
        
        
        
    if(side == 'l'):
        A = leftA + alnStartA
        
        if direction == 1: B = leftB + alnStartB
        if direction == -1: B = len(seqB) - (rightB - alnStartB)
            
        if(aln[:100].count('|') < 90):
            return None
        
    if(side == 'r'):
        A = rightA - (len(subSeqA) - alnEndA)
    
        if direction == 1: B = rightB - (len(subSeqB) - alnEndB)
        if direction == -1: B = len(seqB) - (leftB + len(subSeqB) - alnEndB)
        
        if(aln[-100:].count('|') < 90 and side == 'r'):
            return None
        
    splitSize=6
    splitA = seqA[max(0, A-splitSize):A] + "-" + seqA[A:A+splitSize]
    B_ = B
    if direction == -1:
        B_ = len(seqB) - B_

    splitB = seqB[max(0, B_-splitSize):B_] + "-" + seqB[B_:B_+splitSize]
    if direction == -1:
        splitB = reverse_complement(splitB)
    
    if printAlignment:
        print(splitA)
        print(splitB)
        print(splitA == splitB)
        
    if splitA != splitB:
        return None
    
    return (A, 1, B, direction)
    
def align_left(seqA, seqB, posA, posB, alignBuffer=2000, verbose=False):
    extension=0
    retryCount = 0
    while True:
        leftCut = align(seqA, seqB, posA, posB, alignBuffer, extension, 'l', verbose)
        if leftCut is not None:
            if retryCount > 0:
                print(retryCount)
            return leftCut
        else:
            retryCount = retryCount + 1
            extension=max(alignBuffer, extension*2)
            if verbose:
                print("trying left again with extension=" + str(extension))
                input()
            
def align_right(seqA, seqB, posA, posB, alignBuffer=2000, verbose=False):
    extension=0
    retryCount = 0
    while True:
        rightCut = align(seqA, seqB, posA, posB, alignBuffer, extension, 'r', verbose)
        if rightCut is not None:
            if retryCount > 0:
                print(retryCount)
            return rightCut
        else:
            retryCount = retryCount + 1
            extension=max(alignBuffer, extension*2)
            if verbose:
                print("trying right again with extension=" + str(extension))
                input()


def create_bubble(rid, qid, rSeq, qSeq, rstart, rend, qstart, qend, alignBuffer=2000, verbose=False):
    
    if verbose:
        print("Starting alignment, left")
        input()
    
    #block start
    (qst, qdirSt, rst, rdirSt) = align_left(qSeq, rSeq, qstart, rstart, alignBuffer, verbose)
    
    if verbose:
        print("Starting alignment, right")
        input()
    
    
    #block end
    (qed, qdirEd, red, rdirEd) = align_right(qSeq, rSeq, qend, rend, alignBuffer, verbose) 
    
    if qst > qed:
        print("bad alignment??")

    startFork = Fork(qid, qst, qdirSt, rid, rst, rdirSt)
    startFork.switch_reference()
    endFork = Fork(qid, qed, qdirEd, rid, red, rdirEd)
    endFork.switch_query()


    #import time
    #time.sleep(5)
    return (startFork, endFork)



