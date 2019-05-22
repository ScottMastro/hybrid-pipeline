import aligner
import sys
sys.path.append('../')
import log
from path_helper import Path

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

def fill_Ns_side(side, rid, qid, rSeq, qSeq, block, qpos, param):
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

        startFork, leftReport = fill_Ns_side('l', rid, qid, rSeq, qSeq, block, qstart, param)
        endFork, rightReport = fill_Ns_side('r', rid, qid, rSeq, qSeq, block, qend, param)

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
    
    
    
