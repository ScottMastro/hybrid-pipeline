#!/usr/bin/env python

import ctypes as ct
import ssw_lib
import os

#requires libssw.so to be in the same directory as this script
ssw = ssw_lib.CSsw(os.path.dirname(os.path.realpath(__file__)))
#---------------------------


lEle = ['A', 'C', 'G', 'T', 'N']
dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
nMatch = 4
nMismatch = -6
nOpen = 4
nExt = 1

dEle2Int = {}
dInt2Ele = {}
for i,ele in enumerate(lEle):
    dEle2Int[ele] = i
    dEle2Int[ele.lower()] = i
    dInt2Ele[i] = ele
    
nEleNum = len(lEle)
lScore = [0 for i in range(nEleNum**2)]
for i in range(nEleNum-1):
    for j in range(nEleNum-1):
        if lEle[i] == lEle[j]:
            lScore[i*nEleNum+j] = nMatch
        else:
            lScore[i*nEleNum+j] = nMismatch
            
mat = (len(lScore) * ct.c_int8) ()
mat[:] = lScore


def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num


def align_one(ssw, qProfile, rNum, nRLen, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, 2, 0, 0, ct.c_int32(int(nMaskLen)))

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)


def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sA = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sA += ''.join(['|' if q[nQOff+j] == r[nROff+j] else '*' for j in range(n)])
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += r[nROff : nROff+n]
            nROff += n

    return sCigar, sQ, sA, sR


def reverse_complement(seq):
    return ''.join([dRc[x] for x in seq[::-1]])

def align(qSeq, rSeq, qId, rId, reverse=False, printResults=False):

# build query profile
    qNum = to_int(qSeq, lEle, dEle2Int)
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(qSeq)), mat, nEleNum, 2)
# build rc query profile
    if reverse:
        rqSeq = reverse_complement(qSeq)
        rqNum = to_int(rqSeq, lEle, dEle2Int)
        rqProfile = ssw.ssw_init(rqNum, ct.c_int32(len(rqSeq)), mat, nEleNum, 2)
        
# set mask len
    if len(qSeq) > 30:
        nMaskLen = len(qSeq) / 2
    else:
        nMaskLen = 15

# iter target sequence
    rNum = to_int(rSeq, lEle, dEle2Int)

# format ofres: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
    res = align_one(ssw, qProfile, rNum, len(rSeq), nMaskLen)
# align rc query
    resRc = None
    if reverse:
        resRc = align_one(ssw, rqProfile, rNum, len(rSeq), nMaskLen)

# build cigar and trace back path
    strand = 0
    if resRc == None or res[0] > resRc[0]:
        resPrint = res
        strand = 0
        sCigar, sQ, sA, sR = buildPath(qSeq, rSeq, res[4], res[2], res[8])
    else:
        resPrint = resRc
        strand = 1
        sCigar, sQ, sA, sR = buildPath(rqSeq, rSeq, resRc[4], resRc[2], resRc[8])

    rbegin=resPrint[2]
    rend=resPrint[3]
    qbegin=resPrint[4]
    qend=resPrint[5]

# print results
    if printResults:
        print(sCigar)
        print('target_name: {}\nquery_name: {}\noptimal_alignment_score: {}\t'.format(rId, qId, resPrint[0]))
        if resPrint[1] > 0:
            print('suboptimal_alignment_score: {}\t'.format(resPrint[1]))
        if strand == 0:
            print('strand: +\t')
        else: 
            print('strand: -\t')
        if resPrint[2] + 1:
            print('target_begin: {}\t'.format(resPrint[2] + 1))
        print('target_end: {}\t'.format(resPrint[3] + 1))
        if resPrint[4] + 1:
            print('query_begin: {}\t'.format(resPrint[4] + 1))
        print('query_end: {}\n'.format(resPrint[5] + 1))
        if resPrint[-2] > 0:
            n1 = 1 + resPrint[2]
            n2 = min(60,len(sR)) + resPrint[2] - sR.count('-',0,60)
            n3 = 1 + resPrint[4]
            n4 = min(60,len(sQ)) + resPrint[4] - sQ.count('-',0,60)
            for i in range(0, len(sQ), 60):
                print( 'Target:{:>8}\t{}\t{}'.format(n1, sR[i:i+60], n2))
                n1 = n2 + 1
                n2 = n2 + min(60,len(sR)-i-60) - sR.count('-',i+60,i+120)
    
                print ('{: ^15}\t{}'.format('', sA[i:i+60]))
    
                print('Query:{:>9}\t{}\t{}\n'.format(n3, sQ[i:i+60], n4))
                n3 = n4 + 1
                n4 = n4 + min(60,len(sQ)-i-60) - sQ.count('-',i+60,i+120)

    ssw.init_destroy(qProfile)
    if reverse:
        ssw.init_destroy(rqProfile)

    direction = 1 if strand == 0 else -1

    if reverse:
        return [sQ, qbegin, qend, sA, sR, rbegin, rend, direction]
        
    return [sQ, qbegin, qend, sA, sR, rbegin, rend]