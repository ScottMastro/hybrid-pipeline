import aligner
import sys
sys.path.append('../')
import log

def anchor_ends_side(side, path, block, seqData, param):
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


def anchor_ends(path, contig, seqData, param):
    
    if len(contig.mblocks) < 1:
        return path
   
    startFork, startReport = anchor_ends_side('l', path, contig.mblocks[0], seqData, param)
    endFork, endReport = anchor_ends_side('r', path, contig.mblocks[-1], seqData, param)

    if startFork is not None:
        startFork.switch_query()
        path.add_fork_front(startFork)

    if endFork is not None:
        endFork.switch_reference()
        path.add_fork(endFork)
        
    report = log.ReportSet(log.END_ANCHOR_ATTEMPT, [startReport, endReport])
    param.add_report(report)

    return path


        
    
    
    
