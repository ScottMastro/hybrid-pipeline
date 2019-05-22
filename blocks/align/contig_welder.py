import aligner
from fill_n import fill_Ns
from anchor_ends import anchor_ends

import sys
sys.path.append('../')
import log
from path_helper import Path
from path_helper import clean_path

def create_bubble(rid, qid, rSeq, qSeq, rstart, rend, qstart, qend, param, alignBuffer=500):

    startReport = log.Report(log.ALIGNMENT_ATTEMPT)
    startReport.add_detail(log.QID, str(qid))
    startReport.add_detail(log.RID, str(qid))
    startReport.add_detail(log.SIDE, log.LEFT)
    startReport.add_detail_list(log.QPOS, qstart)
    startReport.add_detail_list(log.RPOS, rstart)

    #block start    
    log.out("Starting alignment, left.", 3, param, wait=True)
    startFork = aligner.align_left(rid, qid, rSeq, qSeq, rstart, qstart, param, alignBuffer)

    endReport = log.Report(log.ALIGNMENT_ATTEMPT)
    endReport.add_detail(log.QID, str(qid))
    endReport.add_detail(log.RID, str(qid))
    endReport.add_detail(log.SIDE, log.LEFT)
    endReport.add_detail_list(log.QPOS, qend)
    endReport.add_detail_list(log.RPOS, rend)

    #block end
    log.out("Starting alignment, right.", 3, param, wait=True)
    endFork = aligner.align_right(rid, qid, rSeq, qSeq, rend, qend, param, alignBuffer)
    
    if startFork is None or endFork is None:
        return (None, startReport, None, endReport)
    
    startFork.switch_reference()
    endFork.switch_query()

    return (startFork, startReport, endFork, endReport)

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
                (startFork, startReport, endFork, endReport) = create_bubble(rid, qid, \
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
                    
                report = log.ReportSet(log.BLOCK_FILL_ATTEMPT, [startReport, endReport])
                param.add_report(report)
                    
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
    path = anchor_ends(path, contig, seqData, param)                        

    cleanPath = clean_path(path, lengthData, param)
    return cleanPath

        
    
    
    
