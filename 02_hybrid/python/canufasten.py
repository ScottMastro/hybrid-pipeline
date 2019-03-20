import aligner
from segment import Segment
from segment import join_paths
from segment import Path
from segment import try_join
from segment import pp

import copy
from contig_matrix import construct_matrix

def scaffold_left(novaId, canuId, novaSeq, canuSeq, canuPos, dir, align_buffer=1000, base_buffer=200):

    novaSeq_ = novaSeq[: min(align_buffer, len(novaSeq)-1)]

    canuLeft = max(0, int(canuPos-align_buffer/2))
    canuRight = min(int(canuPos+align_buffer/2), len(canuSeq)-1)
    canuSeq_ = canuSeq[canuLeft:canuRight]
    if dir == -1:
        canuSeq_ = aligner.reverse_complement(canuSeq_)

    canuAln, canuStart, canuEnd, aln, novaAln, novaStart, novaEnd = \
        aligner.align(canuSeq_, novaSeq_, "canu", "nova", reverse=False, printResults=False)

    buffer = base_buffer
    while buffer < len(aln):
        if aln[buffer-1] == "|":
            break
        else:
            buffer = buffer + 1

    ngapsNova = novaAln[:buffer].count("-")
    ngapsCanu = canuAln[:buffer].count("-")
    
    if ngapsNova + ngapsCanu > (base_buffer*2) * 0.08:
        print("WARNING: poor alignment at left end of supernova contig")
        return Path()

    cl = 0
    cr = canuLeft + canuStart + buffer - ngapsCanu
    nl = min(novaStart + buffer - ngapsNova, len(novaSeq)-1)
    nr = len(novaSeq)
    
    if dir == -1:
        cr = cr - canuLeft + len(canuSeq)-canuRight
    
    canuSeg = Segment(canuId, cl, cr, dir, len(canuSeq))
    novaSeg = Segment(novaId, nl, nr, 1, len(novaSeq))

    return Path([canuSeg, novaSeg])  


def scaffold_right(novaId, canuId, novaSeq, canuSeq, canuPos, dir, align_buffer=1000, base_buffer=200):
    novaLeft = len(novaSeq) - min(align_buffer, len(novaSeq)-1)
    novaSeq_ = novaSeq[novaLeft:]
    canuLeft = max(0, int(canuPos-align_buffer/2))
    canuRight = min(int(canuPos+align_buffer/2), len(canuSeq)-1)
    canuSeq_ = canuSeq[canuLeft:canuRight]
    if dir == -1:
        canuSeq_ = aligner.reverse_complement(canuSeq_)

    canuAln, canuStart, canuEnd, aln, novaAln, novaStart, novaEnd = \
        aligner.align(canuSeq_, novaSeq_, "canu", "nova", reverse=False, printResults=False)

    #verify we cut at a non-gap position
    #buffer bases will be included as canu
    buffer = base_buffer
    while buffer < len(aln):
        if aln[-buffer] == "|":
            break
        else:
            buffer = buffer + 1

    ngapsNova = novaAln[-buffer:].count("-")
    ngapsCanu = canuAln[-buffer:].count("-")
    
    if ngapsNova + ngapsCanu > (buffer*2) * 0.08:
        print("WARNING: poor alignment at right end of supernova contig")
        return Path()

    nl = 0
    nr = novaLeft + novaEnd - buffer + ngapsNova
    cl = canuLeft + canuEnd - buffer + ngapsCanu
    cr = len(canuSeq)
    
    if dir == -1:
        cl = cl - canuLeft + len(canuSeq) - canuRight
    
    '''
    size=20
    size2=10
    print(novaSeq[nr-size:nr] + "-"*size2 + " " + novaId)
    if dir == 1:
        print(canuSeq[cl-size:cl + size2] + " " + canuId)
    else:
        print(aligner.reverse_complement(canuSeq)[cl-size:cl + size2] + " " + canuId)
    '''
    
    novaSeg = Segment(novaId, nl, nr, 1, len(novaSeq))
    canuSeg = Segment(canuId, cl, cr, dir, len(canuSeq))
    
    return Path([novaSeg, canuSeg])  


def find_overlap(blockdf, novaId, canuId, chunkSize=1000, overlapThreshold=6000):
    ''' 
    find overlap between the left and right ends of the nova contig
    and the canu contig. The end of the matching block must be at least
    overlapThreshold bp away from the end(s) of the nova tig. Returns 3 values:
    (leftmost matching base, rightmost matching base, direction of alignment)
    returns None if a match is not found
    '''
    
    blocks = blockdf[(blockdf["contig"] == novaId) & (blockdf["chrom"] == canuId)]
    blocks = blocks.sort_values(by=["chunk_start"])

    leftBlock = blocks.iloc[0]
    rightBlock = blocks.iloc[-1]
    
    dir = leftBlock["direction"]
    if rightBlock["direction"] != dir:
        print("WARNING: blocks align in two different directions")
        return (None, None, None)
        
    size = (rightBlock["total_chunk"] -1) * chunkSize

    novaLeft = (leftBlock["chunk_start"]-1)*chunkSize
    novaRight = (rightBlock["chunk_end"]-1)*chunkSize

    if dir == -1:
        canuLeft = leftBlock["right"]
        canuRight = rightBlock["left"]
    elif dir == 1:
        canuLeft = leftBlock["left"]
        canuRight = rightBlock["right"]
    
    leftOverlap = novaLeft <= overlapThreshold
    rightOverlap = (size - novaRight) <= overlapThreshold

    return (canuLeft if leftOverlap else None, canuRight if rightOverlap else None, dir)


def build_network(novaTigs, canuTigs, novaData, canuData, blockdf, 
                  blockThreshold=5, chunkSize=1000,
                  overlapThreshold=6000, alignBuffer=8000, baseBuffer=200):
    ''' 
    Creates a node object for each contig and connects them if they can be
    scaffolded together.
    '''

    #build matrix
    matchDict = construct_matrix(novaTigs, canuTigs, blockdf, threshold=blockThreshold, asDict=True)
    
    pathDict = dict()

    for novaId in novaTigs:
        novaSeq = str(novaData[novaId])
        matches = matchDict[novaId]
    
        for canuId in matches:
            #check if nova, canu pair overlaps at the end of the nova contig
            left, right, dir = find_overlap(blockdf, novaId, canuId, chunkSize, overlapThreshold) 
            
            if left is None and right is None:
                continue
            
            canuSeq = str(canuData[canuId])
    
            leftPath, rightPath = Path(), Path()
            
            if left is not None:
                leftPath = scaffold_left(novaId, canuId, novaSeq, canuSeq, left, dir, alignBuffer, baseBuffer)
            if right is not None:
                rightPath = scaffold_right(novaId, canuId, novaSeq, canuSeq, right, dir, alignBuffer, baseBuffer)

            #combines left and right paths if canu spans all of nova
            path = join_paths(leftPath, rightPath)
            if not path.is_empty():
                
                pathDict.setdefault(path.get_left_id(), []).append(path)

                if path.get_left_id() != path.get_right_id():
                    pathDict.setdefault(path.get_right_id(), []).append(path)
    
    return pathDict
        

def collapse_internal(contigs, pathDict, preferContig):
    
    leftPathDict = copy.deepcopy(pathDict)
    

    
    for tigId in contigs:
        
        if tigId in leftPathDict:
            paths = leftPathDict[tigId]
            
            if(len(paths) == 2):
                
                p0 = copy.deepcopy(paths[0])
                p1 = copy.deepcopy(paths[1])
                
                #case 1: 
                #(tig00000569_pilon_pilon 0 91449 -1.0) (0 196 5455970 1) (tig00000569_pilon_pilon 5520142 29667941 -1.0) 
                #(tig00000569_pilon_pilon 0 5998502 -1.0) (6 200 39335996 1) 
                result = try_join(paths[0], paths[1])
                
                if result is not None:
                    p0.print_path()
                    p1.print_path()
                    print("result:")
                    result.print_path()
                    leftPathDict[tigId] = [result]






            
            if len(leftPathDict[tigId]) > 1:
                for pth in leftPathDict[tigId]:
                    pth.print_path()
                    
            print("-----------------")

    
def get_joinable(pathDict, path, pid):
    
    joinSet = []
    
    if pid in pathDict: 
        joinSet.extend(pathDict[pid])
    
    result = []
    for join in joinSet:
        
        if try_join(path, join, getResult=False):
            result.append(join)
    
    return result

def join(pathDict, pathA, pid):
    joinSetA = get_joinable(pathDict, pathA, pid)
    if len(joinSetA) == 1:
        pathB = joinSetA[0]
        joinSetB = get_joinable(pathDict, pathB, pid)
        
        if len(joinSetB) == 1:
            if joinSetA[0].equals(pathB) and joinSetB[0].equals(pathA):
                return(pathA, pathB)
                
    return None

def remove_from_dict(pathDict, path):
    
    removed = []
    
    for leftPath in pathDict[path.get_left_id()]:
        if not path.equals(leftPath):
            removed.append(leftPath)
    
    pathDict[path.get_left_id()] = removed
            
    if path.get_left_id() != path.get_right_id():
        removed = []

        for rightPath in pathDict[path.get_right_id()]:
            if not path.equals(rightPath):
                removed.append(rightPath)
        
        pathDict[path.get_right_id()] = removed

    return pathDict

def update_dict(pathDict, pathA, pathB, joined):
    
    pathDict = remove_from_dict(pathDict, pathA)
    pathDict = remove_from_dict(pathDict, pathB)

    pathDict.setdefault(joined.get_left_id(), []).append(joined)
    
    if joined.get_left_id() != joined.get_right_id():
        pathDict.setdefault(joined.get_right_id(), []).append(joined)

    return pathDict

def collapse_paths(contigs, pathDict, preferContig):
    
    copyDict = copy.deepcopy(pathDict)
    
    
    for tigId in contigs:
        
        if tigId in copyDict:
            
            i=0
            while True:
                paths = copyDict[tigId]   
                
                if i >= len(paths):
                    break
                
                path=paths[i]
                
                leftJoin = join(copyDict, path, path.get_left_id())
                
                if leftJoin is not None:
                    joined = try_join(leftJoin[0], leftJoin[1])
                    update_dict(copyDict, leftJoin[0], leftJoin[1], joined)
                    i=0
                    continue
                    
                if path.get_left_id() != path.get_right_id():
                    rightJoin = join(copyDict, path, path.get_right_id())
                    if rightJoin is not None:
                        joined = try_join(rightJoin[0], rightJoin[1])
                        update_dict(copyDict, rightJoin[0], rightJoin[1], joined)
                        i=0
                        continue

                i=i+1

            

    #oldstd=sys.stdout
    import sys
    sys.stdout = open("before.txt", "w")

    for tigId in pathDict:
        for path in pathDict[tigId]:
            path.print_path()
        
    sys.stdout = open("after.txt", "w")

    for tigId in copyDict:
        for path in copyDict[tigId]:
            path.print_path()
    
    sys.stdout = oldstd
    

    
    for tigId in contigs:   
         
        if tigId not in copyDict:
            continue
        
        if len(copyDict[tigId])>1:
            pp(copyDict[tigId])
            print ("=====")

            