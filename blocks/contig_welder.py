from path import Path
from path import clean_path
import log

import aligner

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

def fill_Ns(rid, qid, rSeq, qSeq, block, Nstart, Nend, param):

    alignBuffer = 500
    retryIncrease = 1500
    retryAllowed = 3
    
    sameStrand = (block.get_dir(q=False) == block.get_dir(q=True))

    minPos = min(Nstart, Nend)
    maxPos = max(Nstart, Nend)
    pathN = Path()
    
    log.out("Finding Ns from position " + str(minPos) + " to " + str(maxPos) + ".", 1, param)
    
    for qstart, qend in find_Ns(qSeq, minPos, maxPos):
        #print(str(qstart) +"-" + str(qend) + " " + str(Nend- Nstart))
        
        log.out("Ns found between position " + str(qstart) + " to " + str(qend) + ".", 2, param)        
        log.out("Attempting to fill Ns...", 2, param)
        
        startFork, endFork = None, None

        retryAttempt=0  
        while retryAttempt < retryAllowed:
            
            rstart = block.closest_corresponding_position(qstart, q=True, side='l')
            qs = block.closest_corresponding_position(rstart, q=False)

            if abs(qs - qstart) > 20000:           
                log.out("No corresponding sequence in reference was found on left side. Skipping N filling.", 1, param)
                break

            startFork = aligner.align_left(rid, qid, rSeq, qSeq, \
                rstart, qstart, param, alignBuffer)

            #check if aligments make sense
            redo = False
            if startFork is None: redo = False
            else:   
                sameStrandStart = startFork.before_strand() == startFork.after_strand()
                redo = sameStrand != sameStrandStart
            if redo:
                qstart = qstart - retryIncrease
                startFork = None
                log.out("Shifting start position and retrying alignment.", 1, param)
                retryAttempt = retryAttempt + 1
                continue
            
            #everything is good
            log.out("Left alignment successful.", 3, param)
            break
       
        if startFork is None:
            log.out("N filling failed on left side. Skipping.", 1, param)
            continue

        retryAttempt=0  
        while retryAttempt < retryAllowed:

            rend = block.closest_corresponding_position(qend, q=True, side='r')
            qe = block.closest_corresponding_position(rend, q=False)
                        
            if abs(qe - qend) > 20000:           
                log.out("No corresponding sequence in reference was found on right side. Skipping N filling.", 1, param)
                break

            endFork = aligner.align_right(rid, qid, rSeq, qSeq, \
                rend, qend, param, alignBuffer)

            #check if aligments make sense            
            if endFork is None: redo = True
            else:
                sameStrandEnd = endFork.before_strand() == endFork.after_strand()
                redo = sameStrand != sameStrandEnd
            if redo:
                qend = qend + retryIncrease
                endFork = None
                log.out("Shifting end position and retrying alignment.", 1, param)
                retryAttempt = retryAttempt + 1
                continue

            #everything is good
            log.out("Right alignment successful.", 3, param)
            break
              
        if endFork is None:
            log.out("N filling failed on right side. Skipping.", 1, param)
            continue

        if startFork is not None and endFork is not None:
            startFork.switch_reference()
            pathN.add_fork(startFork)
            endFork.switch_query()
            pathN.add_fork(endFork)
            log.out("N filling successful.", 2, param)

    log.out("N filling complete.", 2, param)

    return pathN


def weld_ends(path, contig, seqData, param):
    
    alignBuffer=100
    retryIncrease = 500
    retryAllowed = 5

    if len(contig.mblocks) < 1:
        return path
    
    startBlock = contig.mblocks[0]
    endBlock = contig.mblocks[-1]

    rid = startBlock.rid
    qid = startBlock.qid
    qSeq = seqData[str(qid)]
    rSeq = seqData[str(rid)]

    rstart = startBlock.start(q=False)
    qstart = startBlock.start(q=True)
   
    log.out("Aligning start of contig " + str(qid) + " to " + str(rid) + ".", 1, param)

    retryAttempt=0  
    while retryAttempt < retryAllowed:
        
        if retryAttempt > 0:
            qstart = qstart + retryIncrease
            rstart = startBlock.closest_corresponding_position(qstart, q=True, side='l')
            log.out("Shifting start position and retrying alignment.", 1, param)


        startFork = aligner.align_right(rid, qid, rSeq, qSeq, rstart, qstart, \
                                  param, alignBuffer)
        
        if startFork is not None:
            break
        else:
            retryAttempt=retryAttempt+1

    if startFork is not None:
        startFork.switch_query()
        path.add_fork_front(startFork)
        log.out("Start position alignment successful.", 2, param)
    else:
        log.out("Start position alignment failed. Skipping.", 1, param)

    rid = endBlock.rid
    qid = endBlock.qid
    qSeq = seqData[str(qid)]
    rSeq = seqData[str(rid)]
    
    rend = endBlock.end(q=False)
    qend = endBlock.end(q=True)
    
    log.out("Aligning end of contig " + str(qid) + " to " + str(rid) + ".", 1, param)

    retryAttempt=0  
    while retryAttempt < retryAllowed:
        
        if retryAttempt > 0:
            qend = qend - retryIncrease
            rend = startBlock.closest_corresponding_position(rend, q=True, side='r')
            log.out("Shifting end position and retrying alignment.", 1, param)

        endFork = aligner.align_left(rid, qid, rSeq, qSeq, rend, qend, \
                                     param, alignBuffer)
        
        if endFork is not None:
            break
        else:
            retryAttempt=retryAttempt+1


    if endFork is not None:
        endFork.switch_reference()
        path.add_fork(endFork)
        log.out("End position alignment successful.", 2, param)
    else:
        log.out("End position alignment failed. Skipping.", 1, param)

    return path

def weld_contig(contig, seqData, param):
    path = Path()
    alignBuffer=500

    '''
    megablock = contig.mblocks[0]
    prevBlock = megablock[1]
    block = megablock[2]
    '''
    
    for megablock in contig.mblocks:
        prevBlock = None
        rid = megablock.rid
        qid = megablock.qid
        qSeq = seqData[str(qid)]
        rSeq = seqData[str(rid)]
    
        log.out("Merging contigs " + str(rid) + " and " + str(qid) + ".", 1, param)
        log.out("----------------------------------------------------", 1, param)

        for block in megablock:
            if prevBlock is not None:
                
                log.out("Filling space between chunk " + str(prevBlock.right_id()) + " and chunk " + \
                        str(block.left_id()) + ".", 1, param)
                
                rstart = prevBlock.end(q=False)
                rend = block.start(q=False)
                qstart = prevBlock.end(q=True)
                qend = block.start(q=True)
                
                #aligns the boundary points between blocks
                #creates a ref/query bubble for the space between blocks
                (startFork, endFork) = aligner.create_bubble(rid, qid, \
                    rSeq, qSeq, rstart, rend, qstart, qend, param, alignBuffer)
                    
                #todo: retry alignment?
                #alignment should work out though...
                if startFork is None:
                    print("failure in start fork!!!!!!!!!")
                    
                if endFork is None:
                    print("failure in end fork!!!!!!!!!")

                if startFork is not None and endFork is not None:

                    log.out("Successfully replaced gap.", 2, param)

                    # find and fill Ns
                    Nstart = 0 if path.last_fork() is None else path.last_fork().get_pos(q=True)
                    Nend = startFork.get_pos(q=True)
                    pathN = fill_Ns(rid, qid, rSeq, qSeq, prevBlock, Nstart, Nend, param)
                    path.add_path(pathN)
                    
                    path.add_fork(startFork)
                    path.add_fork(endFork)
                    
            prevBlock = block
            
        # find and fill Ns to the end of the megablock
        Nstart = 0 if path.last_fork() is None else path.last_fork().get_pos(q=True)
        Nend = megablock.right(q=True)
        pathN = fill_Ns(rid, qid, rSeq, qSeq, prevBlock, Nstart, Nend, param)
        path.add_path(pathN)
            
        log.out("----------------------------------------------------", 1, param)

    return path

    
def weld(contig, seqData, lengthData, param):
    path = weld_contig(contig, seqData, param)
    path = weld_ends(path, contig, seqData, param)                        

    cleanPath = clean_path(path, lengthData, param)
    return cleanPath

        
    
    
    
