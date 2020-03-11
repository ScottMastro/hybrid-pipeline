import sys
sys.path.append('..')
import log as logger

from structures.path import Path
from structures.fork import Fork
import weld.aligner as aligner


def fill_diffs(qSeq, rSeq, block, param):
    '''
    Identifies smaller differences between query and reference and creates Forks 
    around the query sequence between using the chunks to guide the alignments.
    Returns a Path that fills in diffs.
    '''
    #TODO: output log info
    
    pcidThresh = 94
    anchorSize = 100
    gapDetectSize = 5
    buffer = 60
    printAlignments=False
    qid = block.qid
    rid = block.rid
    pathDiff = Path()

    for i,chunk in enumerate(block[1:]):
        i = i + 1
        
        if chunk.id -1 != block[i-1].id:
            chunk1, chunk2 = block[i-1], chunk
            
            qRange = (chunk1.right(q=True) - anchorSize,  chunk2.left(q=True) + anchorSize)
            
            r1 = min( max(chunk1.left(q=False), chunk1.right(q=False)), \
                      max(chunk2.left(q=False), chunk2.right(q=False)) )
            r2 = max( min(chunk1.left(q=False), chunk1.right(q=False)), \
                      min(chunk2.left(q=False), chunk2.right(q=False)) )
            if r1 >= r2: 
                r1 = min( min(chunk1.left(q=False), chunk1.right(q=False)), \
                          min(chunk2.left(q=False), chunk2.right(q=False)) )
                r2 = max( max(chunk1.left(q=False), chunk1.right(q=False)), \
                          max(chunk2.left(q=False), chunk2.right(q=False)) )
            rRange = (r1 - anchorSize, r2 + anchorSize)
            rReverse = (chunk1.get_dir(q=False) == -1)
            
            alignment = aligner.ssw_align(qSeq, rSeq, qRange, rRange, False, rReverse, printAlignments)
            
            x1 = alignment.alignment[1][:anchorSize].count("|")
            x2 = alignment.alignment[1][-anchorSize:].count("|")
            
            if 1.0*x1 /anchorSize > pcidThresh or 1.0*x2 /anchorSize > pcidThresh:
                continue   
            
            start, end = None, None
            count = 0
            for l,cigar in alignment.iter_cigar:
                if (cigar == 'I' or cigar == 'D') and l >= gapDetectSize:
                    if start is None:
                        start = count            
                    end = count + l
                count += l
            if start is None or end is None: continue
            leftCutResult, rightCutResult = None, None  
            while start >= buffer:
                start = start - buffer
                leftCutResult = aligner.cut_alignment(qSeq, rSeq, alignment, 'l', param, oneTry=start)
                if leftCutResult is not None: break
            
            if leftCutResult is None: continue

            length = len(alignment.alignment[1])
            while end <= length-buffer:
                end = end + buffer
                rightCutResult = aligner.cut_alignment(qSeq, rSeq, alignment, 'r', param, oneTry=end)
                if rightCutResult is not None: break
                            
            if rightCutResult is None: continue
  
            lq, lr = leftCutResult
            rq, rr = rightCutResult
    
            leftFork =  Fork(qid, lq, 1, rid, lr, -1 if rReverse else 1)
            rightFork = Fork(qid, rq, 1, rid, rr, -1 if rReverse else 1)
            leftFork.switch_reference()
            rightFork.switch_query()
            
            #print(qSeq[lq:rq])
            #print(rSeq[lr:rr])

            pathDiff.add_fork(leftFork)
            pathDiff.add_fork(rightFork)
            
    return pathDiff