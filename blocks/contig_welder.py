from path import Path
from path import clean_path

import aligner

def weld_ends(path, contig, seqData, param, verbose=False):
            
    startBlock = contig.mblocks[0]
    endBlock = contig.mblocks[-1]

    rid = startBlock.rid
    qid = startBlock.qid
    qSeq = seqData[str(qid)]
    rSeq = seqData[str(rid)]
    
    rstart = startBlock.start(q=False)
    qstart = startBlock.start(q=True)

    startFork = aligner.align_single_right(rid, qid, rSeq, qSeq, rstart, qstart, \
                                  alignBuffer=100, verbose=verbose)

    rid = endBlock.rid
    qid = endBlock.qid
    qSeq = seqData[str(qid)]
    rSeq = seqData[str(rid)]
    
    rend = endBlock.end(q=False)
    qend = startBlock.end(q=True)

    endFork = aligner.align_single_left(rid, qid, rSeq, qSeq, rend, qend, \
                                  alignBuffer=100, verbose=verbose)
    
    if len(path) == 0:
        if startFork is not None:
            startFork.switch_query()
            path.add_fork(startFork)
        if endFork is not None:
            endFork.switch_reference()
            path.add_fork(endFork)
        return path

    if startFork is not None:
        if startFork.after_id() == path.first_fork().before_id():
            path.add_fork_front(startFork)
        else:
            startFork.flip_switch()
            if startFork.after_id() == path.first_fork().before_id():
                path.add_fork_front(startFork)
    if endFork is not None:
        if path.last_fork().after_id() == endFork.before_id():
            path.add_fork(endFork)
        else:
            endFork.flip_switch()
            if path.last_fork().after_id() == endFork.before_id():
                path.add_fork(endFork)
        
    return path


def find_Ns(sequence, minIdx=-1, maxIdx=-1):
    """
    Find the (start,end) of Ns in sequence
    """
    idx = minIdx
    startIdx = -1
    lastIdx = -1

    while True:
        idx = sequence.find('N', idx+1)
        
        if maxIdx > 0 and idx > maxIdx:
            if startIdx < lastIdx:
                yield (startIdx, lastIdx)
            break

        if(idx != lastIdx + 1):
            if startIdx == -1:
                startIdx = idx
            else:           
                yield (startIdx, lastIdx)
                startIdx = idx
            
        if idx == -1:
            break
        
        lastIdx = idx

def weld_contig(contig, seqData, param, verbose=False):
    path = Path()
    alignBuffer=500
    
    def fill_Ns(block, Nstart, Nend):

        minPos = min(Nstart, Nend)
        maxPos = max(Nstart, Nend)
        pathN = Path()
        retryIncrease=1500
        sameDir = (block.get_dir(q=False) == block.get_dir(q=True))
        
        for qstart, qend in find_Ns(qSeq, minPos, maxPos):
            #print(str(qstart) +"-" + str(qend) + " " + str(Nend- Nstart))

            retryAttempt=0           
            while retryAttempt < 3:
                startFork, endFork = None, None
                rstart = block.closest_corresponding_position(qstart, q=True, side='l')
                rend = block.closest_corresponding_position(qend, q=True, side='r')
    
                qs = block.closest_corresponding_position(rstart, q=False)
                qe = block.closest_corresponding_position(rend, q=False)
                
                if abs(qs - qstart) > 20000 or abs(qe - qend) > 20000:
                    print("skipping N filling")
                    break
    
                (startFork, endFork) = aligner.create_bubble(rid, qid, rSeq, qSeq, \
                    rstart, rend, qstart, qend, alignBuffer, verbose)
                    
                redo = False
                if sameDir:
                    if startFork.before_strand() != startFork.after_strand():
                        qstart = qstart - retryIncrease
                        redo = True
                    if endFork.before_strand() != endFork.after_strand():
                        qend = qend + retryIncrease
                        redo = True
                if not sameDir:
                    if startFork.before_strand() == startFork.after_strand():
                        qstart = qstart - retryIncrease
                        redo = True
                    if endFork.before_strand() == endFork.after_strand():
                        qend = qend + retryIncrease
                        redo = True
                if startFork is None:
                        qstart = qstart - retryIncrease
                        redo = True
                if endFork is None:
                        qend = qend + retryIncrease
                        redo = True

                if redo:
                    retryAttempt = retryAttempt + 1
                    print(startFork)
                    print(endFork)
                    print("retrying N filling")
                    continue

                break
            
            if startFork is not None and endFork is not None:
                pathN.add_fork(startFork)
                pathN.add_fork(endFork)

        return pathN

    '''
    megablock = contig.mblocks[0]
    prevBlock = megablock[1]
    block = megablock[2]
    verbose=False
    '''
    
    for megablock in contig.mblocks:
        prevBlock = None
        rid = megablock.rid
        qid = megablock.qid
        qSeq = seqData[str(qid)]
        rSeq = seqData[str(rid)]
    
        for block in megablock:
            if prevBlock is not None:
                rstart = prevBlock.end(q=False)
                rend = block.start(q=False)
                qstart = prevBlock.end(q=True)
                qend = block.start(q=True)
                
                '''
                print(">novaseq: " + str(qstart) + "-" + str(qend))
                print (qSeq[qstart-3000:qstart])            
                print (qSeq[min(qstart,qend):max(qstart,qend)])
                print (qSeq[qend:qend+3000])
                print ("")
                print ("")
    
                print(">canuseq: " + str(rstart) + "-" + str(rend))
                print (rSeq[rstart-100:rstart])
                print (rSeq[min(rstart,rend):max(rstart,rend)])
                print (rSeq[rend:rend+100])
                
                input()
                '''
                
                #aligns the boundary points between blocks
                #creates a ref/query bubble for the space between blocks
                (startFork, endFork) = aligner.create_bubble(rid, qid, \
                    rSeq, qSeq, rstart, rend, qstart, qend, alignBuffer,
                    verbose)
                    
                if startFork is not None and endFork is not None:

                    # find and fill Ns
                    Nstart = 0 if path.last_fork() is None else path.last_fork().get_pos(q=True)
                    Nend = startFork.get_pos(q=True)
                    pathN = fill_Ns(prevBlock, Nstart, Nend)
                    path.add_path(pathN)
                    
                    path.add_fork(startFork)
                    path.add_fork(endFork)
    
            prevBlock = block

    return path

    
def weld(contig, seqData, lengthData, param, verbose=False):
    path = weld_contig(contig, seqData, param, verbose)
    path = weld_ends(path, contig, seqData, param, verbose)                        


    cleanPath = clean_path(path, lengthData, verbose)
    return cleanPath

        
    
    
    
