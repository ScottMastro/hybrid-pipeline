import re
import log
from path_helper import valid_fork_pair
from path_helper import path_overlap
from path_helper import path_length
import pandas as pd

def linkedReadsDict(canutigs, arks_output):
    lrdict_dict = dict()
    file = open(arks_output, 'r')
    for line in file:
        sline = line.split('\t')
        first = sline[1]
        second = sline[2]
        tig1_r = first[1:] #tig1_r = tig1_renamed --> contig in the arks renamed fasta file
        tig2_r = second[1:]
        tig1 = canutigs[(int(tig1_r))-1] #finding the original contig from the renamed one
        tig2 = canutigs[(int(tig2_r))-1]
        bars = int(sline[3])
        
        if tig1 not in lrdict_dict.keys():
            lrdict = dict()
            lrdict[tig2] = bars
            lrdict_dict[tig1] = lrdict
        else:
            lrdict_dict[tig1][tig2] = bars
    return lrdict_dict

def linkedReadsDF(canutigs, lengthData, arks_output):
    #Make df with linked reads, currently unused
    lrdf = pd.DataFrame(columns=['tig1', 'tig1len', 'dir2', 'tig2', 'tig2len', 'dir2', 'n_barcodes'])
    file = open(arks_output, 'r')
    i = 0
    freq_dict = dict()
    for tig in canutigs:
        freq_dict[tig] = 0 
    
    for line in file:
        sline = line.split('\t')
        first = sline[1]
        second = sline[2]
        tig1_r = first[1:] #tig1_r = tig1_renamed --> contig in the arks renamed fasta file
        tig2_r = second[1:]
        
        tig1 = canutigs[(int(tig1_r))-1] #finding the original contig from the renamed one
        tig1len = lengthData[tig1]
        dir1 = first[0]
        tig2 = canutigs[(int(tig2_r))-1]
        tig2len = lengthData[tig2]
        dir2 = second[0]
        bars = int(sline[3])
        
        freq_dict[tig1] += bars
        freq_dict[tig2] += bars
        
        lrdf.loc[i] = tig1, tig1len, dir1, tig2, tig2len, dir2, bars
        i += 1
    #tig1_freq and tig2_freq = total number of barcodes for this contig that connect to contigs other than the one in the current row
    list1 = []
    list2 = []
    for index, row in lrdf.iterrows():
        list1.append(freq_dict[row[0]]-row[6]) 
        list2.append(freq_dict[row[3]]-row[6])
    lrdf.insert(2, 'tig1_freq', list1)
    lrdf.insert(6, 'tig2_freq', list2)
    
    file.close()
    lrdf.to_csv('linked_reads.csv')
    
    return lrdf

def checkLinkedReads(path, nextPath, linkedReads, threshold):   
    tig1a = path[0].rid
    tig1b = path[-1].rid
    tig2a = nextPath[0].rid
    tig2b = nextPath[-1].rid
    
    if tig1a in linkedReads.keys():
        if tig2a in linkedReads[tig1a].keys():
            if linkedReads[tig1a][tig2a] >= threshold:
                return True
        if tig2b in linkedReads[tig1a].keys():
            if linkedReads[tig1a][tig2b] >= threshold:
                return True
    if tig1b in linkedReads.keys():
        if tig2a in linkedReads[tig1b].keys():
            if linkedReads[tig1b][tig2a] >= threshold:
                return True
        if tig2b in linkedReads[tig1b].keys():
            if linkedReads[tig1b][tig2b] >= threshold:
                return True
    return False

def output(path):
    (sequence, source) = path_to_sequence(path, seqData)
    (sequenceInvert, sourceInvert) = path_to_sequence(path, seqData, invert=True)

    f = open("canu_corrected2.fasta", "w+")
    g = open("supernova_uncorrected2.fasta", "w+")

    for i in range(len(sequence)):
        if i % 2 == 1 and i != len(sequence):

            print(str(i) + " canu = " + str(len(sequence[i])))
            print(str(i) + " nova = " + str(len(sequenceInvert[i])))
            f.write(">" + str(i) + "_canu\n")
            f.write(re.sub("(.{164})", "\\1\n", sequence[i], 0, re.DOTALL) + "\n")
            g.write(">" + str(i) + "_nova\n")
            g.write(re.sub("(.{164})", "\\1\n", sequenceInvert[i], 0, re.DOTALL) + "\n")

    f.close()
    g.close()

    f = open("hybrid.fasta", "w+")
     
    f.write(">" + str(contig.id) + "\n")
    f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
    
    f.close()

    '''
    g.write(">" + str(megablock.qid) + "_" + str(megablock.rid) + "\n")
    g.write(re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")
    '''
    
    #validate_forks(path, seqData, 25)
    #qforks = path.get_fork_sequence(seqData, 3000, q=True)
    #rforks = path.get_fork_sequence(seqData, 3000, q=False)        
        
    
def scaffold_pair(path, nextPath, side, nextSide, lengthData, param):

    reports =  []
        
    for shift in [0,1,2,3]:
        report = log.Report(log.SCAFFOLD_ATTEMPT)
        reports.append(report)
        idx = 0 if side == 1 else -1
        nextIdx = 0+shift if nextSide == 1 else -1-shift
        flip = side == nextSide
        report.add_detail(log.FLIPPED, flip)
        okay = False
        
        try:
            fork = path[idx]
            nextFork = nextPath[nextIdx]
        except IndexError:
            report.set_fail(log.INVALID_INDEX)
            break
        
        if flip:
            nextFork = nextPath[nextIdx].flip_strands(lengthData, makeCopy=True)
        if side == 2: 

            if fork.after_id() != fork.rid or \
            nextFork.before_id() != nextFork.rid or \
            fork.after_id() != nextFork.before_id():
                report.set_fail(log.INVALID_ID)
                continue

            if valid_fork_pair(fork, nextFork):
                if flip:
                    nextPath.flip_strands(lengthData)
                path.add_path(nextPath)
                print("success" + ( ", shift=" + str(shift) if shift > 0 else ""))
                okay = True
                break
            else:
                report.set_fail(log.JOIN_FAIL)
                continue

                        
        elif side == 1: 
            
            if fork.before_id() != fork.rid or \
            nextFork.after_id() != nextFork.rid or \
            fork.before_id() != nextFork.after_id():
                report.set_fail(log.INVALID_ID)
                continue

            if valid_fork_pair(nextFork, fork):
                if flip:
                    nextPath.flip_strands(lengthData)
                nextPath.add_path(path)
                path = nextPath
                print("success" + ( ", shift=" + str(shift) if shift > 0 else ""))
                okay = True
                break
            else:
                report.set_fail(log.JOIN_FAIL)
                continue
        else:
            report.set_fail(log.INVALID)
            break

    reportSet = log.ReportSet(log.SCAFFOLD_ATTEMPT, reports)
    param.add_report(reportSet)

    import copy
    reportSet.add_detail(log.PATH, copy.deepcopy(path))
    reportSet.add_detail(log.NEXT_PATH, copy.deepcopy(nextPath))
    reportSet.add_detail(log.SIDE, side)
    reportSet.add_detail(log.NEXT_SIDE, copy.deepcopy(nextSide))

    #if okay: continue
    if okay: 
        return path

    pcOverlap = path_overlap(path, nextPath, lengthData, source='r')
    if pcOverlap[1] > 0.1:
        report.set_fail(log.OVERLAP)
        return None
    
    if report.success: report.set_fail(log.UNKNOWN)   
        
    print("unhandled case:")
    print("start" + ("* " if idx == 0 else " ") + str(path[0]))
    print("end" + ("* " if idx == -1 else " ") + str(path[-1]))

    if flip:
        print("nxstartf" + ("* " if nextIdx < 0  else " ") + str(nextPath[-1].flip_strands(lengthData, makeCopy=True)))
        print("nxendf" + ("* " if nextIdx >= 0 else " ") + str(nextPath[0].flip_strands(lengthData, makeCopy=True)))
        print("nxstart "  + str(nextPath[0]))
        print("nxend " + str(nextPath[-1]))

    else:
        print("nxstart" + ("* " if nextIdx >= 0 else " ") + str(nextPath[0]))
        print("nxend" + ("* " if nextIdx < 0 else " ") + str(nextPath[-1]))
        print("nxstartf " + str(nextPath[-1].flip_strands(lengthData, makeCopy=True)))
        print("nxendf " + str(nextPath[0].flip_strands(lengthData, makeCopy=True)))

    print( "overlap: " + str(round(pcOverlap[0],3)) + "  -  " + str(round(pcOverlap[1],3)) )
        
    #import time
    #break
    #time.sleep(5)
    #leftovers.append(nextPath)
        
    return None

def scaffold_all(path, validPaths, tigId, lengthData, normalize, param, linkedReads, arks, threshold):
    leftovers = []
    while len(validPaths) > 0:
        
        print("--------------:")

        pos1 = normalize(path[0])
        pos2 = normalize(path[-1])
        
        i = -1
        smallIndex = -1
        smallDist = 1e9
        side = -1
        nextSide = -1

        #find next closest validPath to either end
        for nextPath in validPaths:      
            i = i+1
            
            if arks:
                if not checkLinkedReads(path, nextPath, linkedReads, threshold):
                    #print("not in linked reads")
                    #input()
                    continue
            
            nextPos1 = normalize(nextPath[0])
            nextPos2 = normalize(nextPath[-1])
                        
            if pos1 is not None:
                if nextPos1 is not None:
                    if abs(pos1 - nextPos1) < smallDist:
                        smallIndex = i
                        smallDist =  abs(pos1 - nextPos1)
                        side, nextSide = 1,1
                if nextPos2 is not None:
                    if abs(pos1 - nextPos2) < smallDist:
                        smallIndex = i
                        smallDist =  abs(pos1 - nextPos2)
                        side, nextSide = 1,2
            if pos2 is not None:
                if nextPos1 is not None:
                    if abs(pos2 - nextPos1) < smallDist:
                        smallIndex = i
                        smallDist =  abs(pos2 - nextPos1)
                        side, nextSide = 2,1
                if nextPos2 is not None:
                    if abs(pos2 - nextPos2) < smallDist:
                        smallIndex = i
                        smallDist =  abs(pos2 - nextPos2)
                        side, nextSide = 2,2
                     
        nextPath = validPaths.pop(smallIndex)
        
        newPath = scaffold_pair(path, nextPath, side, nextSide, lengthData, param)
        if newPath is None: 
            leftovers.append(nextPath)
        else: 
            path = newPath
        
    return (path, leftovers)



def filter_paths(paths, tigId, lengthData, param):
                
    validPaths = []
    invalidPaths = []

    bigIndex = -1
    bigSize = -1

    #find all scaffolds that start or end with tigId
    #keep track of the biggest one
    i=-1
    for path in paths:
        i += 1
        if len(path) < 1:
            continue

        if not path[0].has_id(tigId) and not path[-1].has_id(tigId):
            invalidPaths.append(path)
            continue
        
        validPaths.append(path)
        size = path_length(path)

        if size > bigSize:
            bigSize = size
            bigIndex = len(validPaths)-1
          
    return (validPaths, invalidPaths, bigIndex)


def normalize_pos(fork, tigId, length):
    if fork.is_Nfork(): return None

    pos = fork.get_pos_by_id(tigId)
    if pos is None: return None
    if fork.get_strand_by_id(tigId) == -1:
        pos = length - pos
    return pos

def scaffold(paths, tigId, lengthData, param, linkedReads, arks, threshold, startPath=None):
    
    log.out("Scaffolding to contig: " + str(tigId), 1, param)
    normalize = lambda fork: normalize_pos(fork, tigId, lengthData[tigId])
    
    #1) find paths that can scaffold by tigId
    validPaths, invalidPaths, startHint = filter_paths(paths, tigId, lengthData, param)
    
    if startPath is None:
        if len(validPaths) <= 0:
            return (None, invalidPaths, [])
        path = validPaths.pop(startHint)
    else:
        if len(validPaths) <= 0:
            return (startPath, invalidPaths, [])
        path = startPath

    #2) try to scaffold all validPaths to path
    path, leftovers = scaffold_all(path, validPaths, tigId, lengthData, normalize, param, linkedReads, arks, threshold)
    
    return (path, invalidPaths, leftovers)
    
    

def remove_overlap(masterPath, slavePath, lengthData, source=None):
                                       #'r' for ref, 'q' for query, None = both
    def normalized_pos(fork, tigId):
        pos = fork.get_pos_by_id(tigId)
        if pos is None: return None
        if fork.get_strand_by_id(tigId) == -1:
            pos = lengthData[str(tigId)] - pos
        return pos
    
    def fill_dict(path):
        starts = dict()
        ends = dict()
        for fork in path:
            if source is None :
                tigIds = [fork.before_id(), fork.after_id()]
            else:
                tigIds = ([fork.rid] if source == 'r' else []) + \
                 ([fork.qid] if source == 'q' else [])
                 
            for tigId in tigIds:
                pos = normalized_pos(fork, tigId)
                tigId = str(tigId)
                if tigId not in starts:
                    starts[tigId] = pos
                    ends[tigId] = pos
                else:
                    if pos < starts[tigId]: starts[tigId] = pos
                    if pos > ends[tigId]: ends[tigId] = pos
        return (starts,ends)
               
    starts1, ends1 = fill_dict(masterPath)
    starts2, ends2 = fill_dict(slavePath)

    for tigId in set(starts1.keys()).union(set(starts2.keys())):
        if tigId in starts1 and tigId in starts2:
            start = max(starts1[tigId], starts2[tigId])
            end = min(ends1[tigId], ends2[tigId])

            while len(slavePath) > 0:
                pos = normalized_pos(slavePath[0], tigId)
                if pos is not None:
                    if pos >= start and pos <= end:
                        slavePath.pop(0)
                        continue
                pos = normalized_pos(slavePath[-1], tigId)
                if pos is not None:
                    if pos >= start and pos <= end:
                        slavePath.pop()
                        continue
                break
            

