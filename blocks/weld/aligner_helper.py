import re
import math

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

def check_identity(alignment, cutPos, side, buffer, minPercent):
    '''
    Checks if percent identity around buffer is greater than minPercent.
    Also checks if the rest of the alignment continuing on side is greater than minPercent.
    Returns True if within tolerance, False otherwise.
    '''
    length = len(alignment.alignment[1])
    bufferL = min(length, cutPos - buffer)
    bufferR = max(0, cutPos + buffer)
    
    def check(start, end):
        seq = alignment.alignment[1][start:end]
        percentId = 100.0 * seq.count('|') / len(seq)

        #print(str(percentId) + "% match (min " + str(minPercent) + "%)." )
        #print("------------------------\n")
        #print(alignment.alignment[0][start:end])
        return percentId > minPercent
    
    #check if region around cut site is a good match
    if not check(bufferL, bufferR): return False

    #check if region before/after cut side is a good match
    start, end = None, None
    if side == 'l': end = cutPos
    if side == 'r': start = -cutPos

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
    
    '''
    print("Checking entropy. \n--------------------------\n")
    print(alignment.alignment[2][bufferL:cutPos-kmerBuffer] + "|" + \
          alignment.alignment[2][cutPos-kmerBuffer:cutPos+kmerBuffer] + "|" + \
          alignment.alignment[2][cutPos+kmerBuffer:bufferR])
    print("Partial entropy = " + str(cutEntropy))
    print("Full entropy = " + str(fullEntropy))
    if fullEntropy <= cutEntropy:
        print("Low entropy")
        #input()
    '''
    
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
        
        '''
        print("Checking if cut point is perfect.")
        print(s1)
        print(s2)
        '''
        if len(s1) < smallBuffer + 1 or s1 != s2: continue
        
        s1 = clean_bases(alignment.alignment[0][bufferL + shift:bufferR + shift])
        s2 = clean_bases(alignment.alignment[2][bufferL + shift:bufferR + shift])

        if not check_entropy(alignment, cutPos + shift, buffer): continue

        return shift + cutPos
    
    return None



'''
#alignment check

    print("----------------")
    print("query")
    #print(qSubSeq + " q0")
    check1 = qSubSeq[alignment.query_begin:alignment.query_end+1] == \
          alignment.alignment[2].replace("-", "")
    print(check1)
    if not check1:
        print(alignment.alignment[2].replace("-", ""))
        print(qSubSeq[alignment.query_begin:alignment.query_end+1])
        
    check2 = (qSeq[alignment.qstart:alignment.qend] if not qReverse else \
            helper.reverse_complement(qSeq)[alignment.qstart:alignment.qend]) == \
            alignment.alignment[2].replace("-", "")
    
    print(check2)
    if not check2:
        print(alignment.alignment[2].replace("-", ""))
        print(qSeq[alignment.qstart:alignment.qend])
    
        
    print("ref")
    check3 = rSubSeq[alignment.reference_begin:alignment.reference_end+1] == \
      alignment.alignment[0].replace("-", "")
    print(check3)
    if not check3:
        print(alignment.alignment[0].replace("-", ""))
        print(rSubSeq[alignment.reference_begin:alignment.reference_end+1])
        
    check4 = (rSeq[alignment.rstart:alignment.rend]  if not rReverse else \
            helper.reverse_complement(rSeq)[alignment.rstart:alignment.rend]) == \
          alignment.alignment[0].replace("-", "")
    
    print(check4)
    if not check4:
        print(alignment.alignment[0].replace("-", ""))
        print(rSeq[alignment.rstart:alignment.rend])
    
    if not check1 or not check2 or not check3 or not check4:
        input()

    print(qSeq[alignment.qstart:alignment.qstart+20] if not qReverse else \
            helper.reverse_complement(qSeq)[alignment.qstart:alignment.qstart+20])
    print(rSeq[alignment.rstart:alignment.rstart+20]  if not rReverse else \
            helper.reverse_complement(rSeq)[alignment.rstart:alignment.rstart+20])
    input()
'''