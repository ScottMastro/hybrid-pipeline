import sys
sys.path.append('./ssw_aligner')
import ssw_aligner as ssw
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

def align(seqA, seqB, posA, posB, alignBuffer, printAlignment=False):
    #todo: speedup? pass direction as hint?
    
    '''
    seqA, seqB = qSeq, rSeq,
    posA, posB = qstart, rstart
    posA, posB = qend, rend

    alignBuffer = 20000
    printAlignment=True
    
    '''

    leftA = max(posA - alignBuffer, 0)
    rightA = min(posA + alignBuffer, len(seqA)-1)
    leftB = max(posB - alignBuffer, 0)
    rightB = min(posB + alignBuffer, len(seqB)-1)
    
    subSeqA = seqA[leftA:rightA]
    subSeqB = seqB[leftB:rightB]

    alnB, alnStartB, alnEndB, aln, \
    alnA, alnStartA, alnEndA, direction = \
        ssw.align(subSeqB, subSeqA, "B", "A", \
                      reverse=True, printResults=printAlignment)
    
    startA = leftA + alnStartA
    endA = rightA - (len(subSeqA) - alnEndA)

    if direction == 1:
        startB = leftB + alnStartB
        endB = rightB - (len(subSeqB) - alnEndB)
    if direction == -1:
        startB = len(seqB) - (rightB - alnStartB)
        endB = len(seqB) - (leftB + len(subSeqB) - alnEndB)
    
    leftCut = (startA, 1, startB, direction)
    rightCut = (endA, 1, endB, direction)
    return (leftCut, rightCut)

def align_left(seqA, seqB, posA, posB, alignBuffer=2000, printAlignment=False):
    (leftCut, rightCut) = align(seqA, seqB, posA, posB, alignBuffer, printAlignment)
    return leftCut
def align_right(seqA, seqB, posA, posB, alignBuffer=2000, printAlignment=False):
    (leftCut, rightCut) = align(seqA, seqB, posA, posB, alignBuffer, printAlignment)
    return rightCut


def create_bubble(rid, qid, rSeq, qSeq, rstart, rend, qstart, qend, alignBuffer=2000, printAlignment=False):
        
    #block start
    (qst, qdirSt, rst, rdirSt) = align_left(qSeq, rSeq, qstart, rstart, alignBuffer, printAlignment)
    #block end
    (qed, qdirEd, red, rdirEd) = align_right(qSeq, rSeq, qend, rend, alignBuffer, printAlignment) 
    
    if qst > qed:
        print("bad alignment??")

    startFork = Fork(qid, qst, qdirSt, rid, rst, rdirSt)
    startFork.switch_reference()
    endFork = Fork(qid, qed, qdirEd, rid, red, rdirEd)
    endFork.switch_query()

    return (startFork, endFork)



