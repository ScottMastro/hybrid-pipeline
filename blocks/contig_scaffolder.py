import re
import log
from path import can_join_forks
from path import path_overlap

from path import Path
from path import Node
from path import Graph

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
        
    
def scaffold_pair(mainPath, extendPath, tigId, lengthData, param):
    #assumes mainPath is normalized
    
    extendPath.normalize(tigId, lengthData)       

    log.out("Scaffolding pair: " + str(tigId), 3, param)
    log.out(str(mainPath.__repr__()), 3, param)
    log.out(str(extendPath.__repr__()), 3, param)

    mainFork = mainPath.last_fork()
    extraFork = None
    if mainFork.after_id() != tigId:
        extraFork = mainPath.pop()
        mainFork = mainPath.last_fork()
        
    if mainFork.after_id() != tigId:
        #this should not happen
        log.out("Error while scaffolding contig - ID does not match. Skipping.", 1, param)
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

def add_to_graph(graph, path, tigId, lengthData, param):
    
    #graph = Graph(paths[0])
    
    length = lengthData[tigId]
    #step 1. grab the first fork on the path
    #----------------------
    i=0
    pathPos = path[i].get_pos_by_id(tigId)
    
    #tigId not in path
    if pathPos is None:
        return (graph, path)
    
    if path[i].get_strand_by_id(tigId) == -1:
        pathPos = length - pathPos
    
    #step 2. grab the first node on the graph
    #----------------------
    currentNode = graph.get_start()
 
    graphPos = currentNode.fork.get_pos_by_id(tigId)
    if currentNode.fork.get_strand_by_id(tigId) == -1:
        graphPos = length - graphPos     

    closestNode = currentNode
    closestNodePos = graphPos

    #step 3. iterate the graph, find the closest node to the first path fork
    #----------------------
    while True:
        
        graphPos = currentNode.fork.get_pos_by_id(tigId)
        if currentNode.fork.get_strand_by_id(tigId) == -1:
            graphPos = length - graphPos     
            
        if abs(closestNodePos - pathPos) > abs(graphPos - pathPos):
            closestNode = currentNode
            closestNodePos = graphPos
            
        if len(currentNode.nextNode) < 1:
            print("all the way to the end!")
            break
        currentNode = currentNode.nextNode[0]
        print(currentNode)


    #step 4. try to anchor the first path fork to the graph
    #----------------------
    failToAnchor = False    
    while True:
        flip = False

        graphPos = closestNode.fork.get_pos_by_id(tigId)
        fork = path[i]
        
        #if path and graph are on opposite strands
        # path should merge back into graph
        if closestNode.fork.get_strand_by_id(tigId) != fork.get_strand_by_id(tigId):
            fork = fork.flip_strands(lengthData, makeCopy=True)
            flip = True
            forkPos = fork.get_pos_by_id(tigId)  

            if graphPos < forkPos:
                closestNode = closestNode.nextNode[0]
                continue
            
            if closestNode.fork.before_id() != tigId:
                if closestNode.prevNode[0].fork.get_pos_by_id(tigId) > forkPos:
                    closestNode = closestNode.prevNode[0]
                else:
                    closestNode = closestNode.nextNode[0]
                continue
            
            if fork.after_id() != tigId:
                i = i+1
                continue

        #if path and graph are on same strands
        # graph should merge into path
        else:
            forkPos = fork.get_pos_by_id(tigId)  
            if forkPos < graphPos:
                closestNode = closestNode.prevNode[0]
                continue
            
            if closestNode.fork.after_id() != tigId:
                if closestNode.nextNode[0].fork.get_pos_by_id(tigId) < forkPos:
                    closestNode = closestNode.nextNode[0]
                else:
                    closestNode = closestNode.prevNode[0]
                continue
            if fork.before_id() != tigId:
                i = i+1
                continue
        
        break
            
    #step 4a. anchored start of path, need to anchor end of path
    if not failToAnchor and not flip:
        
        closestStartNode = closestNode
        closestEndNode = closestNode
    
        j=i     
        lastAnchor = (None, None)
                
        while j < len(path):
        
            if path[j].after_id() != tigId:
                j = j +1
                continue
                
            if closestEndNode.fork.before_id() != tigId:
                if len(closestEndNode.nextNode) > 0:
                    closestEndNode = closestEndNode.nextNode[0]
                    continue
                else: break
                
            if closestEndNode.fork.get_pos_by_id(tigId) > path[j].get_pos_by_id(tigId):
                lastAnchor = (j, closestEndNode)
                j = j +1
                continue
            else:
                if len(closestEndNode.nextNode) > 0:
                    closestEndNode = closestEndNode.nextNode[0]
                    continue
                else: break

        j, closestEndNode = lastAnchor     
        if j is None: j = len(path)-1
                
    #step 5. merge path and graph
    
    node_i = None
    node_j = None
    
    k=0
    prevNode = None
    while k < len(path):
        node = Node(path[k])
        if prevNode is not None:
            prevNode.add_next(node)
            node.add_prev(prevNode)
        prevNode = node
        if k == i: node_i = node
        if k == j: node_j = node
        k = k + 1
   
    lastPathPos = path[j].get_pos_by_id(tigId)
    nd = closestStartNode
    switchMain = True
    
    while True:
        pos = nd.fork.get_pos_by_id(tigId)
        if pos is None:
            break
        if nd.fork.get_strand_by_id(tigId) == path[j].get_strand_by_id(tigId):
            if pos > lastPathPos:
                switchMain = False
                break
        
        if len(nd.nextNode) < 1: break
        else: nd = nd.nextNode[0]
    
    node_i.add_prev_main(closestStartNode)
    if switchMain:
        closestStartNode.add_next_main(node_i)
        if closestEndNode is not None:
            node_j.add_next(closestEndNode)
            closestEndNode.add_prev_main(node_j)
    else:
        closestStartNode.add_next(node_i)
        if closestEndNode is not None:
            node_j.add_next_main(closestEndNode)
            closestEndNode.add_prev(node_j)
        
    return (graph, None)
    
def scaffold(paths, tigId, lengthData, param):
    
    log.out("Scaffolding to contig: " + str(tigId), 1, param)
    
    graph = Graph(paths[4])
    ret = []
    for path in paths[5:]:  
    
        graph, ph = add_to_graph(graph, path, tigId, lengthData, param)

        if ph is not None:
            ret.append(ph)

def scaffold2(paths, tigId, lengthData, param, startPath=None):

    log.out("Scaffolding to contig: " + str(tigId), 1, param)
    length = lengthData[tigId]

    def normalizedPos(fork):
        pos = fork.get_pos_by_id(tigId)
        if pos is None: return None
        if fork.get_strand_by_id(tigId) == -1:
            pos = length - pos
        return pos
    
    validPaths = []
    invalidPaths = []

    bigIndex = -1
    bigSize = -1

    for path in paths:
        if len(path) < 1:
            invalidPaths.append(path)
            continue

        pos1 = normalizedPos(path[0])
        pos2 = normalizedPos(path[-1])
        if pos1 is None and pos2 is None:
            invalidPaths.append(path)
            continue
        validPaths.append(path)
        
        if pos1 is not None and pos2 is not None and \
            abs(pos1 - pos2) > bigSize:
                bigSize = abs(pos1 - pos2)
                bigIndex = len(validPaths)-1
          
    if startPath is None:
        if len(validPaths) <= 0:
            return (None, invalidPaths, [])
        path = validPaths.pop(bigIndex)
    else:
        if len(validPaths) <= 0:
            return (startPath, invalidPaths, [])
        path = startPath
        
    leftovers = []
    while len(validPaths) > 0:
        
        print("--------------:")

        
        pos1 = normalizedPos(path[0])
        pos2 = normalizedPos(path[-1])
        
        i = -1
        smallIndex = -1
        smallDist = 1e9
        side = -1
        nextSide = -1

        for nextPath in validPaths:      
            i = i+1
            nextPos1 = normalizedPos(nextPath[0])
            nextPos2 = normalizedPos(nextPath[-1])
                        
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

        okay = False
        for shift in [0,2]:
                
            idx = 0 if side == 1 else -1
            nextIdx = 0+shift if nextSide == 1 else -1-shift
            flip = side == nextSide
            
            try:
                fork = path[idx]
                nextFork = nextPath[nextIdx]
            except IndexError:
                okay = True
                print("shifted whole path: " + str(shift))
                break            
            
            if flip:
                nextFork = nextPath[nextIdx].flip_strands(lengthData, makeCopy=True)
    

    
            if side == 2: 
    
                if fork.after_id() != fork.rid:
                    okay = True
                    print("end anchoring failed for path" )
                    break            
                if nextFork.before_id() != nextFork.rid:
                    okay = True
                    print("end anchoring failed for nextpath" )
                    break       
                if fork.after_id() != nextFork.before_id():
                    okay = True
                    print("contigs don't match" )
                    break

                if can_join_forks(fork, nextFork):
                    if flip:
                        nextPath.flip_strands(lengthData)
                    path.add_path(nextPath)
                    print("success" + ( ", shift=" + str(shift) if shift > 0 else ""))
                    okay = True
                    break
                            
            elif side == 1: 
                
                if fork.before_id() != fork.rid:
                    okay = True
                    print("end anchoring failed for path" )
                    break            
                if nextFork.after_id() != nextFork.rid:
                    okay = True
                    print("end anchoring failed for nextpath" )
                    break
                
                if fork.before_id() != nextFork.after_id():
                    okay = True
                    print("contigs don't match" )
                    break
                

                if can_join_forks(nextFork, fork):
                    if flip:
                        nextPath.flip_strands(lengthData)
                    nextPath.add_path(path)
                    path = nextPath
                    print("success" + ( ", shift=" + str(shift) if shift > 0 else ""))
                    okay = True
                    break
                
            else:
                okay = True
                print("contigs don't match" )
                break
                
            pcOverlap = path_overlap(path, nextPath, lengthData, source='r')
            if pcOverlap[1] > 0.1:
                print( "overlap case: " + str(round(pcOverlap[0],3)) + "  -  " + str(round(pcOverlap[1],3)) )
                okay = True
                break
            
            if fork.rstrand != nextFork.rstrand:
                print("strand mismatch: ")
                print(str(fork))
                print(str(nextFork))
                print("nextPath=" + str(nextPath.__repr__()))
                okay = True
                break

            
        if okay: continue
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

        leftovers.append(nextPath)
        
        #import time
        #break
        #time.sleep(5)
        
    return (path, invalidPaths, leftovers)

def etc(paths, tigId, lengthData, param):

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

    if len(scaff) > 1:
    
        scaff.sort(key=lambda scaff: -scaff[0])

        mainPath = scaff.pop()[2]
        mainPath.normalize(tigId, lengthData)
        while len(scaff) > 0:
            extendPath = scaff.pop()[2]
            mainPath, leftover = scaffold_pair(mainPath, extendPath, tigId, lengthData, param)
            if leftover is not None:
                ret.append(leftover)

        ret.append(mainPath)
        
        
    return ret
        
    
