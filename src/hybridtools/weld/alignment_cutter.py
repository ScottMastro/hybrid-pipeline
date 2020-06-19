import sys
sys.path.append('../../')
import utils.log as logger

import re
import math

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

#==================================================
# Functions to ensure productive alignment cut
#==================================================

def check_identity(alignment, cutPos, side, buffer, minPercent):
    '''
    Checks if percent identity around buffer is greater than minPercent.
    Also checks if the rest of the alignment continuing on side is greater than minPercent.
    Returns True if within tolerance, False otherwise.
    '''
    length = len(alignment.alignment[1])
    bufferL = max(0, cutPos - buffer)
    bufferR = min(length, cutPos + buffer)
    
    def check(start, end):
        seq = alignment.alignment[1][start:end]
        percentId = 100.0 * seq.count('|') / len(seq)

        return percentId > minPercent
    
    #check if region around cut site is a good match
    if not check(bufferL, bufferR): return False

    #check if region before/after cut side is a good match
    start, end = None, None
    if side == 'l': end = cutPos
    if side == 'r': start = -(length-cutPos)

    if check(start, end): return True
    else: return False

def count_kmers(seq, k):
    '''
    Counts the number of times a kmer appears in seq. 
    Returns dictionary of kmer -> count 
    '''
    kmers = dict()
    for i in range(len(seq) - k + 1):
        kmr = seq[i:i + k]

        if kmr in kmers: kmers[kmr] += 1
        else: kmers[kmr] = 1
    return kmers

binDict = {"A": "00",  "T": "01", "C": "10","G": "11", \
           "a": "00",  "t": "01", "c": "10","g": "11", }
def clean_bases(seq): return re.sub('[^ATCGatcg]+', '', seq)

def kmer_entropy(seq, k=4):
    '''
    Returns the kmer entropy of a given sequence.
    '''
    letters = clean_bases(seq)
    kmers = count_kmers(letters, k)
    count = {int(''.join([binDict[nt] for nt in kmer]), 2): n \
             for kmer,n in kmers.items()}
    total = sum(kmers.values())
    entropy = 0
    for x in range(k**4):
        if x in count:
            p_x = float(count[x])/total
            if p_x > 0:
                entropy += - p_x*math.log(p_x, 2)
    return entropy

def check_entropy(alignment, cutPos, buffer):
    '''
    Compares entropy of kmers with and without kmers surrounding cutPos.
    If including kmers around cutPos increases entropy, returns True. False otherwise.
    '''
    kmer = 4
    kmerBuffer = 5
    length = len(alignment.alignment[1])
    bufferL = min(length, cutPos - buffer)
    bufferR = max(0, cutPos + buffer)


    seqFull = alignment.alignment[2][bufferL:bufferR]
    seqCut = alignment.alignment[2][bufferL:cutPos-kmerBuffer] + \
             alignment.alignment[2][cutPos+kmerBuffer:bufferR]
    
    fullEntropy = kmer_entropy(seqFull, k=kmer)
    cutEntropy = kmer_entropy(seqCut, k=kmer)

    return fullEntropy > cutEntropy

def ajdust_cut_pos(alignment, cutPos, side, buffer, smallBuffer):
    '''
    Checks alignment around cutPos of size smallBuffer. If alignment is not exact
    match, attempts to shift cutPos. If kmer entropy is low around cutPos (repetitive
    region), attempts to shift cutPos. Returns shifted cutPos or None if no shift was found.
    '''
    length = len(alignment.alignment[1])
    splitShift=[0, 3, 5, 8, 10, 13, 15]
    smallBufferL = min(length, cutPos-smallBuffer )
    smallBufferR = max(0, cutPos+smallBuffer)
    bufferL = min(length, cutPos-buffer)
    bufferR = max(0, cutPos+buffer)
    
    for shift in splitShift:
        shift = shift * (1 if side == 'r' else -1)
    
        s1 = clean_bases(alignment.alignment[0][smallBufferL + shift:smallBufferR + shift])
        s2 = clean_bases(alignment.alignment[2][smallBufferL + shift:smallBufferR + shift])
        
        if len(s1) < smallBuffer + 1 or s1 != s2: continue
        
        s1 = clean_bases(alignment.alignment[0][bufferL + shift:bufferR + shift])
        s2 = clean_bases(alignment.alignment[2][bufferL + shift:bufferR + shift])

        if not check_entropy(alignment, cutPos + shift, buffer): continue

        return shift + cutPos
    
    return None

#==================================================
# Main cut functions
#==================================================

def get_cut_position(alignment, attempt, side, buffer):
    '''
    Gets position where a cut will be attempted. Shifts position further 
    as number of attempts increases. Signals attempt termination if buffer 
    overlaps with the start/end of the alignment. Returns a tuple of
    (cutPos, terminationSignal).
    '''
    length = len(alignment.alignment[1])
    pushFactor = 6
    push = math.floor(length/pushFactor)

    #cut moves closer to end of alignment after every iteration
    cutPos = length - math.ceil(length/2) + attempt*push \
        if side == 'r' else math.floor(length/2) - attempt*push
        
    #check if we have reached the end of alignment
    if side == 'l' and cutPos <= buffer:
        cutPos = buffer
        return (cutPos, True)  
    if side == 'r' and cutPos >= length - buffer:
        cutPos = length - buffer
        return (cutPos, True)  
    
    return (cutPos, False)

def cut_alignment(qSeq, rSeq, alignment, side, param,
                  #todo: make these parameters changable by command line?
                  buffer=50, minPercent=90, smallBuffer=5, oneTry=None):
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
            cutPos, terminate = get_cut_position(alignment, attempt, side, buffer)
        attempt += 1
        logger.Logger().out("Attempting to find cut site at position: " + str(cutPos), 5)

        #check the quality of the alignment in this region
        idPass = check_identity(alignment, cutPos, side, buffer, minPercent)
        if not idPass: 
            logger.Logger().out("Percent identity around cut region too low.", 5)
            continue

        cutPos = ajdust_cut_pos(alignment, cutPos, side, buffer, smallBuffer)
        if cutPos is None: 
            logger.Logger().out("Could not find high-entropy exact match for cut site.", 5) 
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
