import aligner
import sys
sys.path.append('../')
import log
from path_helper import Path

def find_Ns(sequence, minIdx=-1, maxIdx=-1, buffer=1500):
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
                
                if buffer > 0 and 'N' in sequence[lastIdx+1:lastIdx+buffer]:
                    lastIdx = idx
                    continue
                
                yield (startIdx, lastIdx)
                startIdx = idx
            
        if idx == -1:
            break
        
        lastIdx = idx

def fill_Ns_side(side, rid, qid, rSeq, qSeq, block, qpos, extend, param):
    #side is left ('l') or right ('r')
    sideName = "left" if side == "l" else "right"
    sameStrand = (block.get_dir(q=False) == block.get_dir(q=True))

    retryAllowed = 3
    alignBuffer = 200
    retryIncrease = 1500
    retryAttempt = 0

    report = log.Report(log.ALIGNMENT_ATTEMPT)
    report.add_detail(log.QID, str(qid))
    report.add_detail(log.RID, str(qid))
    report.add_detail(log.SIDE, log.LEFT if side == "l" else log.RIGHT)

    qpos = qpos + (extend*(-1 if side == 'l' else 1))


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



def check_N_path(path):
    
    if len(path) < 1: return True
        
    if path[0].is_switch_query():
        print('First fork starts from reference:')
        print(path[0])
        input()
        return False
    
    if path[-1].is_switch_reference():
        print('Last fork ends on reference:')
        print(path[-1])
        input()
        return False
    
    startFork = path[0]
    
    for endFork in path[1:]:
                        
        if startFork.after_id() != endFork.before_id():  
            print("Contig IDs do not match:")
            print(startFork)
            print(endFork)
            input()
            return False
        
        if not startFork.after_strand() == endFork.before_strand():
            print("Strands do not match:")
            print(startFork)
            print(endFork)
            input()
            return False

        if startFork.after_pos() > endFork.before_pos():
            print("Positional issue (going backwards):")
            print(startFork)
            print(endFork)
            input()
            return False
        
        if startFork.before_pos() > endFork.after_pos():
            print("Positional issue (ending backwards):")
            print(startFork)
            print(endFork)
            input()
            return False

        startFork = endFork

    return True

    
    
def clean_N_path(path, param):

    if len(path) < 1: return path
    
    cleanPath = Path()
    startFork = path[0]
    skip = False
    
    for i in range(1, len(path)):
        if skip:
            skip=False
            continue
        
        endFork = path[i]

        if (startFork.after_pos() > endFork.before_pos() or \
            startFork.before_pos() > endFork.after_pos()) and \
            startFork.is_switch_query():

            if len(cleanPath) > 0:
                startFork = cleanPath.pop()
            elif len(path) > i+1:
                startFork = path[i+1]
                skip=True
                
            continue
            

        cleanPath.add_fork(startFork)
        startFork = endFork
        
    cleanPath.add_fork(endFork)
    return cleanPath

def fill_Ns(rid, qid, rSeq, qSeq, block, Nstart, Nend, param):

    minPos = min(Nstart, Nend)
    maxPos = max(Nstart, Nend)
    pathN = Path()
    
    log.out("Finding Ns from position " + str(minPos) + " to " + str(maxPos) + ".", 2, param)
    
    for qstart, qend in find_Ns(qSeq, minPos, maxPos):
        #print(str(qstart) +"-" + str(qend) + " " + str(Nend- Nstart))
        
        log.out("Ns found between position " + str(qstart) + " to " + str(qend) + ".", 2, param)        
        log.out("Attempting to fill Ns...", 2, param)

        startFork, endFork = None, None

        extend = 0
        while True:
    
            startFork, leftReport = fill_Ns_side('l', rid, qid, rSeq, qSeq, block, qstart, extend, param)
            endFork, rightReport = fill_Ns_side('r', rid, qid, rSeq, qSeq, block, qend, extend, param)

            if startFork is not None and endFork is not None:
                startFork.switch_reference()
                endFork.switch_query()

                #check for weird alignment due to repeated segments
                if startFork.after_pos() > endFork.before_pos() or \
                    startFork.before_pos() > endFork.after_pos():
                    extend = extend + max(200, abs(startFork.after_pos() - endFork.before_pos()))
                    continue
            
            break

        report = log.ReportSet(log.NFIX_ATTEMPT, [leftReport, rightReport])
        param.add_report(report)

        if startFork is None or endFork is None:
            continue

        pathN.add_fork(startFork)
        pathN.add_fork(endFork)
        log.out("N filling successful.", 2, param)

    if len(pathN) > 0:
        pathN = clean_N_path(pathN, param)
        #check_N_path(pathN)

        
    log.out("N filling complete.", 2, param)
    return pathN
    
    
    
