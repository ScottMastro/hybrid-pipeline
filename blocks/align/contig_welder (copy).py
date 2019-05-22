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

def fill_Ns_fork(side, rid, qid, rSeq, qSeq, block, qpos, param):
    #side is left ('l') or right ('r')
    sideName = "left" if side == "l" else "right"
    sameStrand = (block.get_dir(q=False) == block.get_dir(q=True))

    retryAllowed = 3
    alignBuffer = 500
    retryIncrease = 1500
    retryAttempt = 0  

    report = log.Report(log.ALIGNMENT_ATTEMPT)
    report.add_detail(log.QID, str(qid))
    report.add_detail(log.RID, str(qid))
    report.add_detail(log.SIDE, log.LEFT if side == "l" else log.RIGHT)

    while retryAttempt < retryAllowed:
        
        rpos = block.closest_corresponding_position(qpos, q=True, side=side)
        qp = block.closest_corresponding_position(rpos, q=False)

        report.add_detail_list(log.QPOS, qpos)
        report.add_detail_list(log.RPOS, rpos)

        if abs(qp - qpos) > 20000:        
            log.out("No corresponding sequence in reference was found on " + sideName + " side. Skipping N filling.", 1, param)
            report.set_fail()
            report.add_detail(log.REASON, log.MISSING)
            return (None, report)
            break

        if side == 'l':
            fork = aligner.align_left(rid, qid, rSeq, qSeq, rpos, qp, param, alignBuffer)
        else:
            fork = aligner.align_right(rid, qid, rSeq, qSeq, rpos, qp, param, alignBuffer)

        #check if aligments make sense
        redo = False
        if fork is None: redo = True
        else:   
            sameStrandStart = fork.before_strand() == fork.after_strand()
            redo = sameStrand != sameStrandStart
            
        if redo:
            qpos = qpos + (retryIncrease*(-1 if side == 'l' else 1))
            log.out("Shifting position and retrying alignment.", 1, param)
            retryAttempt = retryAttempt + 1
            continue
        
        #everything is good
        log.out("Alignment successful (" + sideName + ").", 3, param)
        break
   
    if fork is None:
        log.out("N filling failed on " + sideName + " side. Skipping.", 1, param)
        report.set_fail()
        report.add_detail(log.REASON, log.ALIGNMENT)
    
    return (fork, report)


def fill_Ns(rid, qid, rSeq, qSeq, block, Nstart, Nend, param):

    minPos = min(Nstart, Nend)
    maxPos = max(Nstart, Nend)
    pathN = Path()
    
    log.out("Finding Ns from position " + str(minPos) + " to " + str(maxPos) + ".", 1, param)
    
    for qstart, qend in find_Ns(qSeq, minPos, maxPos):
        #print(str(qstart) +"-" + str(qend) + " " + str(Nend- Nstart))
        
        log.out("Ns found between position " + str(qstart) + " to " + str(qend) + ".", 2, param)        
        log.out("Attempting to fill Ns...", 2, param)

        startFork, endFork = None, None

        startFork, leftReport = fill_Ns_fork('l', rid, qid, rSeq, qSeq, block, qstart, param)
        endFork, rightReport = fill_Ns_fork('r', rid, qid, rSeq, qSeq, block, qend, param)

        report = log.ReportSet(log.NFIX_ATTEMPT, [leftReport, rightReport])
        param.add_report(report)

        if startFork is None or endFork is None:
            continue

        startFork.switch_reference()
        pathN.add_fork(startFork)
        endFork.switch_query()
        pathN.add_fork(endFork)
        log.out("N filling successful.", 2, param)

    log.out("N filling complete.", 2, param)
    return pathN


def weld_end(side, path, block, seqData, param):
    #side is left ('l') or right ('r')
    sideName = "start" if side == "l" else "end"

    alignBuffer=100
    retryIncrease = 500
    retryAllowed = 5
    
    rid = block.rid
    qid = block.qid
    qSeq = seqData[str(qid)]
    rSeq = seqData[str(rid)]

    if side == 'l':
        rpos = block.start(q=False)
        qpos = block.start(q=True)
    else:
        rpos = block.end(q=False)
        qpos = block.end(q=True)
        
    log.out("Aligning " + sideName + " of contig " + str(qid) + " to " + str(rid) + ".", 1, param)

    report = log.Report(log.ALIGNMENT_ATTEMPT)
    report.add_detail(log.QID, str(qid))
    report.add_detail(log.RID, str(qid))
    report.add_detail(log.SIDE, log.LEFT if side == "l" else log.RIGHT)

    retryAttempt=0  
    while retryAttempt < retryAllowed:
        
        if retryAttempt > 0:
            qpos = qpos + (retryIncrease*(1 if side == 'l' else -1))
            rpos = block.closest_corresponding_position(qpos, q=True, side=side)
            log.out("Shifting " + sideName + " position and retrying alignment.", 1, param)

        report.add_detail_list(log.QPOS, qpos)
        report.add_detail_list(log.RPOS, rpos)

        if side == 'l':
            fork = aligner.align_right(rid, qid, rSeq, qSeq, rpos, qpos, param, alignBuffer)
        else:
            fork = aligner.align_left(rid, qid, rSeq, qSeq, rpos, qpos, param, alignBuffer)

        if fork is None:
            retryAttempt=retryAttempt+1
        else:
            break

    if fork is not None:
        log.out("Alignment successful (" + sideName + ").", 2, param)
    else:
        log.out("Alignment failed (" + sideName + "). Skipping.", 1, param)
        report.set_fail()
        report.add_detail(log.REASON, log.ALIGNMENT)

    return (fork, report)


def weld_ends(path, contig, seqData, param):
    
    if len(contig.mblocks) < 1:
        return path
   
    startFork, startReport = weld_end('l', path, contig.mblocks[0], seqData, param)
    endFork, endReport = weld_end('r', path, contig.mblocks[-1], seqData, param)

    if startFork is not None:
        startFork.switch_query()
        path.add_fork_front(startFork)

    if endFork is not None:
        endFork.switch_reference()
        path.add_fork(endFork)
        
    report = log.ReportSet(log.NFIX_ATTEMPT, [startReport, endReport])
    param.add_report(report)

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

        
    
    
    
