import aligner
from contig_matrix import construct_matrix

#functions
#----------------

def find_Ns(sequence):
    """
    Find the (start,end) of Ns in sequence
    """
    idx = -1
    startIdx = -1
    lastIdx = -1

    while True:
        idx = sequence.find('N', idx+1)

        if(idx != lastIdx + 1):
            if startIdx == -1:
                startIdx = idx
            else:           
                yield (startIdx, lastIdx)
                startIdx = idx
            
        if idx == -1:
            break
        
        lastIdx = idx

#finds canu contigs that match region of supernova contig
def match_canu_block(blockdf, matrix, novaId, novaLeft, novaRight, chunk_size):
    
    # lookup matches in matrix
    try:
        matrix_match = set(matrix[novaId].dropna().keys())
    except KeyError:
        return []
    
    blocks = blockdf.loc[(blockdf['contig'] == novaId)]
    if len(blocks) < 1:
        return []
    
    total_chunk=blocks["total_chunk"].iloc[0]
    
    chunk_start = max(1, int(novaLeft/chunk_size) + 1)
    chunk_end = min(total_chunk, int(novaRight/chunk_size) + 1)

    blocks_start = blocks[blocks['chunk_start'] <= chunk_start]
    blocks_start = blocks_start[blocks_start['chunk_end'] >= chunk_start]
    
    blocks_end = blocks[blocks['chunk_start'] <= chunk_end]
    blocks_end = blocks_end[blocks_end['chunk_end'] >= chunk_end]

    return list(set(blocks_start['chrom']) & set(blocks_end['chrom']) & set(matrix_match))

#finds start and end of canu sequence to extract
def get_canu_range(aligndf, novaId, canuId, chunk_size, start, end, buffer=0):
    
    chunk_start = max(1, int(start/chunk_size + 1))
    chunk_end = int(end/chunk_size + 1)
    
    alignments = aligndf[(aligndf[2] == canuId) & (aligndf[0] == int(novaId))]
    chunks = alignments.iloc[0,4].split(",")

    while str(chunk_start) not in chunks:
        chunk_start = chunk_start-1
        if chunk_start < 1:
            return ""

    while str(chunk_end) not in chunks:
        chunk_end = chunk_end+1
        if chunk_end > int(chunks[-1]):
            return ""

    #todo: there are potentially multiple alignments for each block!!!
    #possible solution = take all starts and ends and resolve?? use block information???
        
    if chunks.count(str(chunk_start)) > 1:
        print("WARNING: multiple aligned chunks when getting canu sequence (start)")
    if chunks.count(str(chunk_end)) > 1:
        print("WARNING: multiple aligned chunks when getting canu sequence (end)")

    chunk_start_index = chunks.index(str(chunk_start))
    chunk_end_index = chunks.index(str(chunk_end))
    
    canu_starts = alignments.iloc[0,5].split(",")
    canu_ends = alignments.iloc[0,6].split(",")
   
    canuStart = min( int(canu_starts[chunk_start_index]), int(canu_starts[chunk_end_index]), int(canu_ends[chunk_start_index]), int(canu_ends[chunk_end_index]) )
    canuEnd = max( int(canu_starts[chunk_start_index]), int(canu_starts[chunk_end_index]), int(canu_ends[chunk_start_index]), int(canu_ends[chunk_end_index]) )

    return (canuStart-buffer, canuEnd+buffer)

def trim(novaAln, canuAln, base_buffer):

    while(novaAln[0] == "-" or canuAln[0] == "-"):
        canuAln=canuAln[1:]
        novaAln=novaAln[1:]

    while(novaAln[-1] == "-" or canuAln[-1] == "-"):
        canuAln=canuAln[:-1]
        novaAln=novaAln[:-1]
        
    nidx = novaAln.index("N")
    counter = base_buffer
    left = novaAln[:nidx]
    left_buffer = ""
    index = nidx-1
    while(index >= 0 and counter > 0):
        if not canuAln[index] == "-":
            left_buffer = canuAln[index] + left_buffer
            counter=counter-1
            
        left=left[:-1]
        index=index-1
    
    left = left + left_buffer
    
    counter = base_buffer
    right = novaAln[nidx:]
    right_buffer = ""
    index = nidx-1
    while(index < len(canuAln) and counter > 0):
        if not canuAln[index] == "-":
            right_buffer = right_buffer +  canuAln[index]
            counter=counter-1
            
        right=right[1:]
        index=index+1
        
    right = right_buffer + right
    trimmed = (left + right).replace("N", "").replace("-", "")
    return trimmed
        

def stitch(novaAln, canuAln, pre, post):
    base_buffer=7
    
    pre = pre[:-base_buffer]
    post = post[base_buffer:]

    _pre = pre
    _post = post

    while(len(pre) > 0):
        if pre[0] == novaAln[0]:
            pre = pre[1:]
    
        canuAln = canuAln[1:]
        novaAln = novaAln[1:]
        
    while(len(post) > 0):
        if post[-1] == novaAln[-1]:
            post = post[:-1]
    
        canuAln = canuAln[:-1]
        novaAln = novaAln[:-1]

    canuAln = canuAln.replace("-", "")
    
    #check if canu is smaller than supernova
    if(len(canuAln) <= base_buffer*2 ):
        return trim(novaAln, canuAln, base_buffer)
    return _pre + canuAln + _post


def align_fill(novaId, novaSeq, novaLeft, novaRight, canuId, canuSeq, canuLeft, canuRight):
    
    novaSeq_= novaSeq[novaLeft:novaRight]
    canuSeq_ = canuSeq[canuLeft:canuRight]
    
    mid=int(len(novaSeq_)/2)
    
    # extract prefix and suffix of sequence
    if novaSeq_[mid] != 'N':
       print("WARNING: expected N to be in the middle of the sequence!")
       return None
   
    pre = novaSeq_[0:mid]
    post = novaSeq_[mid:len(novaSeq_)]
    
    canuAlnPre, canuStartPre, canuEndPre, alnPre, novaAlnPre, novaStartPre, novaEndPre, dir= \
        aligner.align(canuSeq_, pre, "canu", "prenova", reverse=True, printResults=False)    
    canuAlnPost, canuStartPost, canuEndPost, alnPost, novaAlnPost, novaStartPost, novaEndPost, postDir =  \
        aligner.align(canuSeq_, post, "canu", "postnova", reverse=True, printResults=False)    
    
    if dir != postDir:
       print("WARNING: pre and post align in opposite directions")
       return None

    if canuStartPre > canuStartPost:
        print("WARNING: pre and post align unexpectedly")
        return None
    
    ngapsNova = novaAlnPre.count("-")
    ngapsCanu = canuAlnPre.count("-")
    
    if ngapsNova + ngapsCanu > (len(novaAlnPre) + len(canuAlnPre)) * 0.08:
        print("WARNING: poor alignment at left end of NNNNs")
        return None

    ngapsNova = novaAlnPost.count("-")
    ngapsCanu = canuAlnPost.count("-")
    
    if ngapsNova + ngapsCanu > (len(novaAlnPost) + len(canuAlnPost)) * 0.08:
        print("WARNING: poor alignment at right end of NNNNs")
        return None
    
    base_buffer=10
 
    #verify we cut at a non-gap position
    #buffer bases will be included as canu
    #--------------------------------------
    bufferPre = base_buffer
    while bufferPre < len(alnPre):
        if alnPre[-bufferPre] == "|":
            break
        else:
            bufferPre = bufferPre + 1
            
    bufferPost = base_buffer
    while bufferPost < len(alnPost):
        if alnPost[bufferPost-1] == "|":
            break
        else:
            bufferPost = bufferPost + 1
    #--------------------------------------
    
    ngapsNovaPost = novaAlnPost[:bufferPost].count("-")
    ngapsNovaPre = novaAlnPre[-bufferPre:].count("-")
    ngapsCanuPre = canuAlnPre[:bufferPre].count("-")
    ngapsCanuPost = canuAlnPost[-bufferPost:].count("-")
    
    nlPre = novaLeft + novaStartPre + base_buffer
    nrPre = novaLeft + novaEndPre - bufferPre + ngapsNovaPre
    
    nlPost = novaLeft + mid + novaStartPost + bufferPost - ngapsNovaPost
    nrPost = novaLeft + mid + novaEndPost - base_buffer

    cl = canuLeft + canuEndPre - bufferPre + ngapsCanuPre
    cr = canuLeft + canuStartPost + bufferPost - ngapsCanuPost
    if dir == -1:
        cl = cl - canuLeft + len(canuSeq)-canuRight
        cr = cr - canuLeft + len(canuSeq)-canuRight

    '''
    size=20
    size2=10
    print(novaSeq[nrPre-size:nrPre] + "-"*size2 + " " + novaId)
    
    if dir == 1:
        print(canuSeq[cl-size:cl + size2] + " " + canuId)
        print("...1...")
        print(canuSeq[cr-size2:cr+size] + " " + canuId)
    else:
        print(aligner.reverse_complement(canuSeq)[cl-size:cl + size2] + " " + canuId)
        print("...-1...")
        print(aligner.reverse_complement(canuSeq)[cr-size2:cr+size] + " " + canuId)
        
    print("-"*size2 + novaSeq[nlPost:nlPost+size] + " " + novaId)
    '''
      
    return [(novaId, nlPre, nrPre, 1), (canuId, cl, cr, dir), (novaId, nlPost, nrPost, 1)]


def fix_Ns(novaTigs, canuTigs, novaData, canuData, blockdf, aligndf,
           novaBuffer, canuBuffer, blockThreshold, chunkSize):
    
    matrix = construct_matrix(novaTigs, canuTigs, blockdf, threshold=blockThreshold, asDict=False)

    info = {x : [] for x in novaTigs}
    nfixes = {x : [] for x in novaTigs}
    
    for novaId in novaTigs:
    
        novaSeq = str(novaData[novaId])
        print("Filling Ns on contig " + novaId)
    
        # iterate over supernova contigs --> find all Ns in each contig
        for start,end in find_Ns(novaSeq):
            
            #extracts subsequence at specified range + buffer
            novaLeft = max(0, start - novaBuffer)
            novaRight = min(end + novaBuffer, len(novaSeq))
            
            #find matching canu contigs using blocks and matrix
            canuContigs = match_canu_block(blockdf, matrix, novaId, novaLeft, novaRight, chunkSize)        
            
            #verify only one contig matches (for now...)
            if(len(canuContigs) > 1):
                info[novaId].append((start,end,"multimatch"))
                continue
            
            if(len(canuContigs) < 1):
                info[novaId].append((start,end,"unmatched"))
                continue
            
            try:
                #get matching canu sequence
                canuId = canuContigs[0]
                canuLeft, canuRight = get_canu_range(aligndf, novaId, canuId, chunkSize, start, end, canuBuffer)
                canuSeq = str(canuData[canuId])
                
                if(canuRight-canuLeft > 20000):
                    info[novaId].append((start,end,"longmatch"))
                    continue
                if(canuRight-canuLeft < 10):
                    info[novaId].append((start,end,"shortmatch"))
                    continue
                
                #fix nova with canu
                fill = align_fill(novaId, novaSeq, novaLeft, novaRight, canuId, canuSeq, canuLeft, canuRight)
                
                if fill is None:
                    info[novaId].append((start,end,"fixerror"))
                    continue
    
                nfixes[novaId].append(fill)
                    
            except ValueError as e:
                print(str(e))
                info[novaId].append((start,end,"error"))
                
    return nfixes