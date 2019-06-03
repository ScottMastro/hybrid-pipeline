import aligner
from fill_n import fill_Ns
from anchor_ends import anchor_ends

import sys
sys.path.append('../')
import log
from path_helper import Path
from path_helper import clean_path
import path_helper

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
    endReport.add_detail(log.SIDE, log.RIGHT)
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

def join_block_paths(startFork, path, endFork, param):
    
    if len(path) < 1:
        if startFork is not None: path.add_fork(startFork)
        if endFork is not None: path.add_fork(endFork)
        return path

    if startFork is not None:
        if (startFork.after_strand() != path[0].before_strand()):
            print('strand issue within block')
            print(startFork)
            print(path[0])
            input()
            
        if (startFork.after_pos() > path[0].before_pos() or \
            startFork.before_pos() > path[0].after_pos()):
            path.pop(0)
        else:
            path.add_fork_front(startFork)

    if endFork is not None:
        if path[-1].after_strand() != endFork.before_strand():
            print('strand issue within block')
            print(endFork)
            print(path[-1])

            input()
            
        if (path[-1].after_pos() > endFork.before_pos() or \
            path[-1].before_pos() > endFork.after_pos()):
            path.pop()
        else:
            path.add_fork(endFork)

    return path

def join_path_list(paths, lengthData, param):
    
    while len(paths) > 0:
        if paths[0] is None or len(paths[0]) < 1: paths = paths[1:]
        else: break
    if len(paths) < 1: return Path()

    path = paths[0]
    for blockPath in paths[1:]:
        
        if blockPath is None or len(blockPath) < 1: continue
        
        if path[-1].after_id() != blockPath[0].before_id():
            print('id issue')
            print(path[-1])
            print(blockPath[0])
            
            if path[-1].is_switch_reference():
                print("removing last fork (resolved):")
                print(path[-1])
                path.pop()
            
            input()
    
        if path[-1].after_strand() != blockPath[0].before_strand():
                        
            print('handle strand between blocks...')
            
            print(path[-1])
            print(blockPath[0])

            print("flipping")
            blockPath.flip_strands(lengthData)
            
         #   if abs(blockPath[0].qpos - blockPath[-1].qpos) <= 6001:
         #       print(abs(blockPath[0].qpos - blockPath[-1].qpos))
          #      print("removing")
          #      continue
               
            #input()

            
           # if path[-1].is_switch_reference() and len(blockPath) <= 2:
           #     print("removing")
               # continue
            
            #input()
            
        if path[-1].after_pos() > blockPath[0].before_pos():
            print('positional issue')
            print(path[-1])
            print(blockPath[0])
            input()

        path.add_path(blockPath)
        
    return path
        
        
def check_path(path):
    
    if len(path) < 1: return True
    
    startFork = path[0]
    
    for endFork in path[1:]:
        if startFork.is_Nfork() or endFork.is_Nfork():
            startFork = endFork
            continue
                                
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
        
        if startFork.before_id() == endFork.after_id() and \
        startFork.before_strand() == endFork.after_strand() and \
        startFork.before_pos() > endFork.after_pos():
            print("Positional issue (ending backwards):")
            print(startFork)
            print(endFork)
            
            if startFork.is_switch_reference():
                print("Might be intentional...")
            else:
                input()
                return False

        startFork = endFork

    return True

def weld_contig(contig, seqData, param):
    paths = []
    alignBuffer=200

    '''
    megablock = contig.mblocks[0]
    prevBlock = megablock[0]
    block = megablock[1]
    '''
    
    for megablock in contig.mblocks:
        prevBlock = None
        prevEndFork = None
        rid = megablock.rid
        qid = megablock.qid
        qSeq = seqData[str(qid)]
        rSeq = seqData[str(rid)]
    
        log.out("Merging contigs " + str(rid) + " and " + str(qid) + ".", 1, param)
        log.out("----------------------------------------------------", 1, param)
        
        for block in megablock:
            startFork, endFork = None, None
            if prevBlock is not None:
                
                log.out("Filling space between chunk " + str(prevBlock.right_id()) + " and chunk " + \
                        str(block.left_id()) + ".", 2, param)
                   
                rstart = prevBlock.end(q=False)
                rend = block.start(q=False)
                qstart = prevBlock.end(q=True)
                qend = block.start(q=True)
             
                #flipBuffer=True
                while True:
                
                    #aligns the boundary points between blocks
                    #creates a ref/query bubble for the space between blocks
                    (startFork, startReport, endFork, endReport) = create_bubble(rid, qid, \
                        rSeq, qSeq, rstart, rend, qstart, qend, param, alignBuffer)
                        

                    if startFork is None or endFork is None:
                        #todo: retry alignment?
                        print("Failure in bridging blocks!!!!!!!!!")
                        startFork = None
                        endFork = None                        

                    if startFork is not None and endFork is not None:
    
                        #if forks result in backwards movement
                       # if startFork.rstrand == endFork.rstrand and startFork.rpos > endFork.rpos:
                       #     dist = abs(startFork.rpos - endFork.rpos)
                      #      qstart = max(0, int(startFork.qpos - dist))
                      #      qend = min(len(qSeq), int(endFork.qpos + dist))
                      #      rstart = megablock.closest_corresponding_position(qstart, q=True, side='l')
                      #      rend = megablock.closest_corresponding_position(qend, q=True, side='r')
                       #     continue
                  #      if startFork.rstrand == endFork.rstrand and startFork.qpos > endFork.qpos:
                 #          dist = abs(startFork.qpos - endFork.qpos)
                  #          qstart = max(0, int(startFork.qpos - dist))
                 #           qend = min(len(qSeq), int(endFork.qpos + dist))
                 #           rstart = megablock.closest_corresponding_position(qstart, q=True, side='l')
                  #          rend = megablock.closest_corresponding_position(qend, q=True, side='r')
                #            continue
                        #if strand gets flipped, give some buffer space
                    #    elif startFork.rstrand != endFork.rstrand and flipBuffer:
                   #         dist = 500
                   #         qstart = max(0, int(startFork.qpos - dist))
                   #         qend = min(len(qSeq), int(endFork.qpos + dist))
                    #        rstart = prevBlock.closest_corresponding_position(qstart, q=True, side='l')
                    #        rend = block.closest_corresponding_position(qend, q=True, side='r')
                    #        flipBuffer = False
                  #          continue
                        
                        log.out("Successfully replaced gap.", 2, param)
                                        
                    # find and fill Ns
                    Nstart = prevBlock.left(q=True)
                    Nend = prevBlock.right(q=True)

                    pathN = fill_Ns(rid, qid, rSeq, qSeq, prevBlock, Nstart, Nend, param)
                    blockPath = join_block_paths(prevEndFork, pathN, startFork, param)
                    
                    paths.append(blockPath)

                    report = log.ReportSet(log.BLOCK_FILL_ATTEMPT, [startReport, endReport])
                    param.add_report(report)
                    break
            
            prevEndFork = endFork
            prevBlock = block
            
        # find and fill Ns to the end of the megablock
        #Nstart = 0 if len(paths) < 1 or paths[-1].last_fork() is None \
        #    else paths[-1].last_fork().get_pos(q=True)
        #Nstart = max(block.left(q=True), Nstart)
        #Nend = megablock.right(q=True)
        Nstart = block.left(q=True)
        Nend = block.right(q=True)

        pathN = fill_Ns(rid, qid, rSeq, qSeq, prevBlock, Nstart, Nend, param)
        blockPath = join_block_paths(prevEndFork, pathN, None, param)

        #print("~~~~~")
        #print(blockPath)
        paths.append(blockPath)
                
        log.out("----------------------------------------------------", 1, param)

    return paths


def weld_contig2(contig, seqData, param):
    paths = []
    alignBuffer=200

    '''
    megablock = contig.mblocks[0]
    block = megablock[0]
    '''
    
    for megablock in contig.mblocks:
        rid = megablock.rid
        qid = megablock.qid
        qSeq = seqData[str(qid)]
        qLen = len(qSeq)
        rSeq = seqData[str(rid)]
    
        log.out("Merging contigs " + str(rid) + " and " + str(qid) + ".", 1, param)
        log.out("----------------------------------------------------", 1, param)
        
        mblockPaths = []
        
        for block in megablock:
            
            rstart = block.start(q=False)
            rend = block.end(q=False)
            qstart = block.start(q=True)
            qend = block.end(q=True)
        
            startFork = aligner.align_right(rid, qid, rSeq, qSeq, rstart, qstart, param, alignBuffer)
            if startFork is not None: startFork.switch_query()
            else: print("Failure in aligning startFork!!!!!!!!!")
            #todo: retry alignment?

            endFork = aligner.align_left(rid, qid, rSeq, qSeq, rend, qend, param, alignBuffer)
            if endFork is not None: endFork.switch_reference()
            else: print("Failure in aligning endFork!!!!!!!!!")
            #todo: retry alignment?

            if startFork is not None and endFork is not None:                        
                log.out("Successfully replaced gap.", 2, param)
                                
            # find and fill Ns
            Nstart = block.left(q=True)
            Nend = block.right(q=True)

            pathN = fill_Ns(rid, qid, rSeq, qSeq, block, Nstart, Nend, param)            
            blockPath = join_block_paths(startFork, pathN, endFork, param)
            
            if len(blockPath) > 0:
                mblockPaths.append(blockPath)
            
        
        startNFlag, endNFlag = False, False
        if len(mblockPaths) > 0:
            firstqPos = mblockPaths[0][0].qpos if mblockPaths[0][0].qstrand == 1 \
                else qLen - mblockPaths[0][0].qpos
            lastqPos = mblockPaths[-1][-1].qpos if mblockPaths[-1][-1].qstrand == 1 \
                else qLen - mblockPaths[-1][-1].qpos
                
            for mblockForks in mblockPaths:
                for fork in mblockForks:
                    qpos = fork.qpos if fork.qstrand == 1 else qLen - fork.qpos
                    if qpos > lastqPos:
                        print("ending with NNN")
                        endNFlag = True
                        break
                    if qpos < firstqPos:
                        print("starting with NNN")
                        startNFlag = True
                        break
                    
        if startNFlag: paths.append(path_helper.get_Npath())
        paths.extend(mblockPaths)
        if endNFlag: paths.append(path_helper.get_Npath())

        log.out("----------------------------------------------------", 1, param)

    return paths

def join_paths2(paths, lengthData, param):
    
    while len(paths) > 0:
        if paths[0] is None or len(paths[0]) < 1: paths = paths[1:]
        else: break
    if len(paths) < 1: return Path()

    if len(paths) == 1 and len(paths[0]) == 2 and \
        paths[0][0].before_id() == paths[0][-1].after_id() and \
        paths[0][0].before_strand() == paths[0][-1].after_strand() and \
        paths[0][0].before_pos() > paths[0][-1].after_pos():
            path = Path()
            return path


    for i in range(len(paths)-1):
        p1 = paths[i]
        p2 = paths[i+1]
        
        if p1 is None or p2 is None or len(p1) < 1 or len(p2) < 1 or \
            p1[-1].is_Nfork() or p2[0].is_Nfork(): 
            continue
        
        p2Flip = p2[-1].flip_strands(lengthData, makeCopy=True)
        p1Flip = p1[0].flip_strands(lengthData, makeCopy=True)

        
        if p1[-1].after_id() == p2[0].before_id() and \
            p1[-1].after_strand() != p2[0].before_strand():

            print('flipping block path...')
            print(p1[-1])
            print(p2[0])
            

            if p1[-1].after_id() == p2Flip.before_id() and \
                p1[-1].after_pos() < p2Flip.before_pos():
                    print("flipping right")
                    paths[i+1].flip_strands(lengthData)
                    
                    if i+2 >= len(paths): continue
                    print(paths)
                    print(i)
                    p3 = paths[i+2]
                    if p3 is None or len(p3) < 1 or p3[0].is_Nfork(): 
                        continue
                    if p3[-1].after_id() != paths[i+1][-1].before_id():
                        paths[i+1].add_fork(path_helper.get_Nfork())

                    
            elif p1Flip.after_id() == p2[0].before_id() and \
                p1Flip.after_pos() < p2[0].before_pos():
                    print("flipping left")
                    paths[i].flip_strands(lengthData)
                    
                    if i < 1: continue
                    p0 = paths[i-1]
                    if p0 is None or len(p0) < 1 or p0[-1].is_Nfork(): 
                        continue
                    if p0[-1].after_id() != paths[i][0].before_id():
                        paths[i].add_fork_front(path_helper.get_Nfork())
            else:
                    print('skipped')
                    input()
     
                    
        elif p1[-1].after_id() == p2[0].before_id() and \
            p1[-1].after_strand() == p2[0].before_strand():
        
            if p1Flip.after_id() == p2Flip.before_id() and \
                p1Flip.after_pos() < p2Flip.before_pos():
                    print("flipping both")
                    paths[i].flip_strands(lengthData)
                    paths[i+1].flip_strands(lengthData)
            
                    if i < 1: continue
                    p0 = paths[i-1]
                    if p0 is None or len(p0) < 1 or p0[-1].is_Nfork(): 
                        continue
                    if p0[-1].after_id() != paths[i][0].before_id():
                        paths[i].add_fork_front(path_helper.get_Nfork())
                        
                    if i > len(paths)-1: continue
                    p3 = paths[i+2]
                    if p3 is None or len(p3) < 1 or p3[0].is_Nfork(): 
                        continue
                    if p3[-1].after_id() != paths[i+1][-1].before_id():
                        paths[i+1].add_fork(path_helper.get_Nfork())

                        

    path = paths[0]
    lastIdx=0
    
    for blockPath in paths[1:]:
        
        if blockPath is None or len(blockPath) < 1: continue
     
        if path[-1].is_Nfork() or blockPath[0].is_Nfork():
            lastIdx=len(path)
            path.add_path(blockPath)
            continue
    
        #trying to switch between two different reference tigs
        #remove both forks
        if path[-1].after_id() != blockPath[0].before_id() and \
            path[-1].is_switch_reference() and blockPath[0].is_switch_query():
                print("popping:")
                print(path[-1])
                print(blockPath[0])
                path.pop()
                blockPath.pop(0)
                
                if len(blockPath) < 1: continue
                if path[-1].is_switch_query() and blockPath[0].is_switch_reference() and \
                    path[-1].after_strand() != blockPath[0].before_strand():
                        path.add_fork(path_helper.get_Nfork())
                        path.add_path(blockPath)
                        continue

    
        if path[-1].after_id() != blockPath[0].before_id():

            
            print('id issue')
            print(path[-1])
            print(blockPath[0])
            
           # if path[-1].is_switch_reference():
           #     print("removing last fork (resolved):")
           #     print(path[-1])
           #     path.pop()
            
            input()
    
        if path[-1].after_strand() != blockPath[0].before_strand():
                        
            print('handle strand between blocks...')
            
            print(path[-1])
            print(blockPath[0])

            '''
            flippedFork = blockPath[-1].flip_strands(lengthData, makeCopy=True)
            if path[-1].after_id() == flippedFork.before_id() and \
                path[-1].after_pos() < flippedFork.before_pos():
                    print("flipping")
                    blockPath.flip_strands(lengthData)
            else:
                print('lastIdx='+str(lastIdx))
                print(path)
                flippedFork = path[lastIdx].flip_strands(lengthData, makeCopy=True)
                if flippedFork.after_id() == blockPath[0].before_id() and \
                    flippedFork.after_pos() < blockPath[0].before_pos():
                        print("flipping prev")
                        toFlip = []
                        while len(path) > lastIdx: toFlip.append(path.pop())
                        for flip in toFlip: 
                            flip.flip_strands(lengthData)
                            path.add_fork(flip)
                else:
                    print('unhandled')
                    input()
            '''
                
            if len(blockPath) > 1 and abs(blockPath[0].qpos - blockPath[-1].qpos) <= 6001:
                print(abs(blockPath[0].qpos - blockPath[-1].qpos))
                print("removing")
                continue
               
            #input()

        if path[-1].after_pos() > blockPath[0].before_pos():
            print('positional issue')
            print(path[-1])
            print(blockPath[0])
            input()
 
        lastIdx=len(path)
        path.add_path(blockPath)
        
    return path
        



def weld(contig, seqData, lengthData, param):
    paths = weld_contig2(contig, seqData, param)
    path = join_paths2(paths, lengthData, param)
    #path = anchor_ends(path, contig, seqData, param)
    ok = check_path(path)
    print(ok)
    
    if not ok:
        input()
        
    return path
    
    #cleanPath = clean_path(path, lengthData, param)
    #return cleanPath

        
    
    
    
