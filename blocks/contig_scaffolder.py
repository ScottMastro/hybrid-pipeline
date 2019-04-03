import re
from path import can_join_forks
from path import Path

def scaffold_nonintersecting(paths, tigId, lengthData, param, verbose=False):
    
    parts, ret = [], []

    def get_pos(fork):
        if fork.before_id() == tigId:
            if fork.before_strand() == 1:
                return fork.before_pos()
            else:
                return lengthData[tigId] - fork.before_pos()
        
        if fork.after_id() == tigId:
            if fork.after_strand() == 1:
                return fork.after_pos()
            else:
                return lengthData[tigId] - fork.after_pos()
        return None

    def path_tuple(path):
        
        startPos = get_pos(path.first_fork())
        endPos = get_pos(path.last_fork())
        
        if startPos is None and endPos is None:
            return None
        elif startPos is None:
            t = (endPos, None, endPos, path)
        elif endPos is None:
            t = (startPos, startPos, None, path)
        else:
            t = (min(startPos,endPos), startPos, endPos, path)            
        return t
        
    for path in paths:
        if path.id_on_end(tigId):
            parts.append(path_tuple(path))
        else:
            ret.append(path)
            
    parts.sort(key=lambda parts: -parts[0])

    while len(parts) > 1:
        p1 = parts.pop()
        p2 = parts.pop()
        
        flip1 = p1[2] is None or p1[0] != p1[1]
        if flip1:
            f1 = p1[-1].first_fork().flip_strands(lengthData, makeCopy=True)
        else:
            f1 = p1[-1].last_fork()
        
        flip2 = p2[1] is None or p2[0] != p2[1]
        if flip2:
            f2 = p2[-1].last_fork().flip_strands(lengthData, makeCopy=True)
        else:
            f2 = p2[-1].first_fork()
        
        trim1 = f1.after_id() != tigId
        trim2 = f2.before_id() != tigId

        pos1 = f1.before_pos() if trim1 else f1.after_pos()
        pos2 = f2.after_pos() if trim2 else f2.before_pos()
        strand1 = f1.before_strand() if trim1 else f1.after_strand()
        strand2 = f2.after_strand() if trim2 else f2.before_strand()
    
        if strand1 == strand2:
            if (strand1 == 1 and pos1 <= pos2) or \
                (strand1 == -1 and pos2 <= pos1):
                    
                    path = p1[-1]
                    path2 = p2[-1]
                    if flip1: path.flip_strands(lengthData)
                    if flip2:  path2.flip_strands(lengthData)

                    if trim1: path.pop()
                    if trim2: path2.pop(0)
                    
                    path.add_path(path2)
                    t = path_tuple(path)
                    if t is None:
                        ret.append(path)
                        continue
                    
                    parts.append(t)
                    continue
                    
        parts.append(p2)
        ret.append(p1[-1])
    
    if len(parts) == 1:
        ret.append(parts[0][-1])

    return ret
        
        
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
        
    
def scaffold_pair(mainPath, extendPath, tigId, lengthData, param, verbose=False):
    #assumes mainPath is normalized
    
    extendPath.normalize(tigId, lengthData)       
    mainFork = mainPath.last_fork()
    extraFork = None
    if mainFork.after_id() != tigId:
        extraFork = mainPath.pop()
        mainFork = mainPath.last_fork()
        
    if mainFork.after_id() != tigId:
        print("error while scaffolding contig " + str(tigId))
        return (mainPath, extendPath)
    
    extended=False
    for i in range(len(extendPath)):
        if can_join_forks(mainFork, extendPath[i]):
            mainPath.extend(extendPath[i:])
            extended=True
            break
    if not extended: i = i+1
    
    if i > 0:
        overlap = Path()
        overlap.extend(extendPath[:i])
        startFork = overlap.first_fork()
        if startFork.before_id() != tigId:
            extraFork = overlap.pop(0)
            startFork = overlap.first_fork()
            
        if startFork.before_id() != tigId:
            print("giving up")

        startIndex, endIndex = None, None
        for j in range(len(extendPath)-i, len(mainPath)):
            if can_join_forks(mainPath[-j], startFork):
                startIndex = len(mainPath)-j
                break
            
        endFork = overlap.last_fork()
        if endFork.after_id() != tigId:
            extraFork = overlap.pop()
            endFork = overlap.last_fork()
            
        if endFork.after_id() != tigId:
            print("giving up")
            
        for k in range(len(mainPath)-j, len(mainPath)):
            if can_join_forks(endFork, mainPath[k]):
                endIndex = k
                break
            
        if startIndex is None and endIndex is None:
             return (mainPath, extendPath)
         
        mainPath.add_supplementary(overlap, startIndex, endIndex)       
        
    return  (mainPath, None)
    
    
    
def scaffold(paths, tigId, lengthData, param, verbose=False):
    
    ret = []
    scaff = []

    def path_tuple(path):
        minPos, maxPos = path.get_interval(tigId, lengthData)
        return (minPos, maxPos, path)

    for path in paths:
        if path.id_on_end(tigId):
            scaff.append(path_tuple(path))
        else:
            ret.append(path)

    scaff.sort(key=lambda scaff: -scaff[0])

    mainPath = scaff.pop()[2]
    mainPath.normalize(tigId, lengthData)
    while len(scaff) > 0:
        extendPath = scaff.pop()[2]
        mainPath, leftover = scaffold_pair(mainPath, extendPath, tigId, lengthData, param, verbose)
        if leftover is not None:
            ret.append(leftover)


    ret.append(mainPath)
    return ret
        
