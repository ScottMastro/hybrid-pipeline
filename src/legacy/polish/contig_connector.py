import supernova_helper as supernova

LEFT='L'
RIGHT='R'
def flip_direction(currentDirection):
    return RIGHT if currentDirection == LEFT else LEFT
    
    
def backtrack(currentUnitigId, traversed, refGraph, direction=RIGHT):
    
    leftTigs, rightTigs = refGraph.get_unitig_neighbours(currentUnitigId)
    backTrackTigs = leftTigs if direction == RIGHT else rightTigs

    if len(backTrackTigs) > 1:
        print(backTrackTigs)
        untraversedTigIds = [i for i in backTrackTigs if i not in traversed]
        
        for untraversedTigId in untraversedTigIds:
    
            l, r = refGraph.get_unitig_neighbours(untraversedTigId)
            possibleBranches = set(l + r)
            possibleBranches.remove(currentUnitigId)
            
            #the untraversed tig has no neighbours other than the current tig
            if len(possibleBranches) == 0:
                #todo:resolve this (get contigs, run structural solver??)
                traversed[untraversedTigId] = True
                print("structural unitig:", untraversedTigId, "continuing with", currentUnitigId)
                print("Backtracked and found structural unitig, moving on...")
            else:
                print("unsure what to do code:backtrack.....")
                break
    
    return traversed

def zero_path(currentContigId, aligndf, traversed, refGraph, queryGraph, direction=RIGHT):
    
    print("No unitig path to follow. Checking supernova for hints....")
    
    qAlignments = aligndf[(aligndf["rid"] == currentContigId)]
    checkSize=2000

    if direction == RIGHT:
        qAlignments = qAlignments.sort_values(by="rend", ascending=False)
        otherAlignments = qAlignments[qAlignments["rend"] > int(qAlignments.iloc[0]["rend"]) - checkSize]
    elif direction == LEFT:
        qAlignments = qAlignments.sort_values(by="rstart", ascending=True)
        otherAlignments = qAlignments[qAlignments["rstart"] < int(qAlignments.iloc[0]["rend"]) - checkSize]

    qTig = qAlignments.iloc[0]
    qTigId = qTig["qid"]
    otherTigIds = list(set(otherAlignments["qid"]))

    # only one contig or a pair of haplotigs are at the end of our contig
    if len(otherTigIds) == 1 or (len(otherTigIds) == 2 and \
               supernova.find_haplotig(otherTigIds[0], queryGraph) == otherTigIds[1]):
    
        qTigAlignments = aligndf[(aligndf["qid"] == qTigId) & (aligndf["rid"] != currentContigId)]
        rContigsNext = list(set(qTigAlignments["rid"]))

        # only one contig or a pair of haplotigs are at the end of our contig
        if len(rContigsNext) == 1:
            print("ohhhhh yeahhhh.")
            
            nextAlignment = qTigAlignments.iloc[0]
            nextContigId = nextAlignment["rid"]
            #todo: attempt to scaffold

            sameStrand = nextAlignment["strand"] == qTig["strand"]
            newDirection = direction if sameStrand else flip_direction(direction)
               
            if newDirection == RIGHT:
                newUnitigId = refGraph.get_closest_unitig(nextContigId, nextAlignment["rstart"])
            elif newDirection == LEFT:
                newUnitigId = refGraph.get_closest_unitig(nextContigId, nextAlignment["rend"])

            print("scaffolding using supernova:", qTigId)
            print("next unitig:", newUnitigId, "on contig:", nextContigId)
            print("next direction:", newDirection)

            return (nextContigId, newUnitigId, newDirection)
            
        # we can do better....... some cases here will be valid
        else:
            print("unsure what to do code:zeropath_2.....")
            #more advanced logic here....
            input()
    

    else:
        print("unsure what to do code:zeropath_1.....")
        #more advanced logic here....
        input()

    
    print("Done.")
    print("(Can we keep going though?)")


def one_path(currentUnitigId, nextUnitigId):
    print("next unitig:", nextUnitigId)
    print("Only one unitig path available.")
    return nextUnitigId

def two_path(currentUnitigId, nextTigs, aligndf, traversed, refGraph, direction=RIGHT):
    
    #check if the two unitigs represent large structural differences:
    nextLeftTigs1, nextRightTigs1 = refGraph.get_unitig_neighbours(nextTigs[0])
    nextLeftTigs2, nextRightTigs2 = refGraph.get_unitig_neighbours(nextTigs[1])

    if len(nextLeftTigs1) > 1 or len(nextRightTigs1) > 1 or len(nextLeftTigs2) > 1 or len(nextRightTigs2) > 1:
        print("unsure what to do code:two_path1.....")
        input()

    #flip if reverse complement:
    if currentUnitigId in nextRightTigs1 and currentUnitigId not in nextLeftTigs1:
        nextTigs1 = nextLeftTigs1
    elif currentUnitigId in nextLeftTigs1 and currentUnitigId not in nextRightTigs1:
        nextTigs1 = nextRightTigs1
    else:
        print("unsure what to do code:two_path2a.....")

    if currentUnitigId in nextRightTigs2 and currentUnitigId not in nextLeftTigs2:
        nextTigs2 = nextLeftTigs2
    elif currentUnitigId in nextLeftTigs2 and currentUnitigId not in nextRightTigs2:
        nextTigs2 = nextRightTigs2
    else:
        print("unsure what to do code:two_path2b.....")


    if len(nextTigs1) == 1 and len(nextTigs2) == 1 and nextTigs1[0] == nextTigs2[0]:
        print("we good dude.")
        print("branching unitigs:", nextTigs, "next unitig", nextTigs1[0])
        print("Structral variants, branching back single unitig.")
        #todo:resolve this
        
        traversed[nextTigs[0]] = True
        traversed[nextTigs[1]] = True

        nextUnitigId = nextTigs1[0]
        return (nextUnitigId, traversed)
    else:
        print("unsure what to do code:two_path3.....")
        input()

def three_path(currentContigId, currentUnitigId, nextTigs, aligndf, traversed, refGraph, queryGraph, direction=RIGHT):
    u2cDict = {i : refGraph.unitig_to_contig(i) for i in nextTigs}
    c2uDict = {}
    
    for key in u2cDict:
        for value in u2cDict[key]:
            if value in c2uDict:
                c2uDict[value].append(key)
            else:
                c2uDict[value] = [key]
    if currentContigId in c2uDict:
        
        if len(c2uDict[currentContigId]) == 1:
            unitigCandidateId = c2uDict[currentContigId][0]
            
            end   = max(refGraph.unitig_pos(currentContigId, currentUnitigId))
            start = min(refGraph.unitig_pos(currentContigId, unitigCandidateId))

            #is this supported by a supernova contig??
            
            qContigsStartIds = list(aligndf[(aligndf["rid"] == currentContigId) & \
                                    (aligndf["rstart"] <= min(end,start))]["qid"])
                
            qContigsEndIds = list(aligndf[(aligndf["rid"] == currentContigId) & \
                                    (aligndf["rstart"] >= max(end,start))]["qid"])

            qContigIds = list(set(qContigsStartIds).intersection(set(qContigsEndIds)))
            
            # only one contig or a pair of haplotigs support this connection
            if len(qContigIds) == 1 or (len(qContigIds) == 2 and \
                   supernova.find_haplotig(qContigIds[0], queryGraph) == qContigIds[1]):
                print("we good son.")
                
            else:
                print("unsure what to do code:three_path1.....")
                input()

        else:
            print("unsure what to do code:three_path2.....")
            input()
        
        
    print("unsure what to do code:three_path3.....")
    input()



#===================================================================

#chr2:22,299,575-22,602,231      'utg00000635' ->  ['tig00000636', 'tig00000670']
def resolve_path(rid, refGraph, queryGraph, aligndf, seqData, lengthData, param):
    currentContigId = rid
    currentUnitigId = refGraph.get_closest_unitig(currentContigId, 0)
    traversed = dict()
    direction=RIGHT
    
    #BIG todo: have more accurate unitig-contig mapping AFTER pilon polishing

    while(True):
        print("------------------------------\n", "current unitig:", currentUnitigId)
        
        leftTigs, rightTigs = refGraph.get_unitig_neighbours(currentUnitigId)
        nextTigs = leftTigs if direction == LEFT else rightTigs
        traversed[currentUnitigId] = True
        #todo: update contigId

        traversed = backtrack(currentUnitigId, traversed, refGraph, direction)

        if len(nextTigs) == 0:
            results = zero_path(currentContigId, aligndf, traversed, refGraph, queryGraph, direction)
            currentContigId, currentUnitigId, direction = results

        elif len(nextTigs) == 1:
            currentUnitigId = one_path(currentUnitigId, rightTigs[0])
            
        elif len(nextTigs) == 2:
            results = two_path(currentUnitigId, nextTigs, aligndf, traversed, refGraph, direction)
            currentUnitigId, traversed = results
            
        else:
            three_path(currentContigId, currentUnitigId, nextTigs, aligndf, traversed, refGraph, queryGraph, direction)





    
    
    '''
    ustart, uend = refGraph.unitig_pos(rid, unitigId)
    unitigRegion = SimpleRegion(rid, ustart, uend)
    
    #todo: utig coordinates are a bit off.....

    while(True):
        leftTigs, rightTigs = refGraph.get_unitig_neighbours(unitigId)
    
        #go right
        if len(rightTigs) == 1:
            unitigId = rightTigs[0]
        else:
            print(unitigId)
            print(rightTigs)
            for utig in rightTigs:
                print(utig + ":")
                print("\t" + ",".join(canu.unitig_to_contig(refGraph, utig)) )
                print("\t" + ",".join(canu.unitig_pos(refGraph, rid, utig)) )

            input()

    
    hapid = supernova.find_haplotig(qid, queryGraph)
    
    #find easiest examples, 2 haplotigs + 1 to 1 alignments
    if hapid is None: 
        print(qid, " no haplotig found")
        return

    alignments = aligndf[(aligndf["qid"] == qid) | (aligndf["qid"] == hapid)]
    rmatches = set(alignments["rid"])
    
    if len(rmatches) > 1:
        print(qid, " matches multiple ref contigs")
        return
    
    print(qid, " is good")

    #todo: watch sort order!!
    for idx, row in (alignments).iterrows():
        #todo: pair haplotig alignments!!!
        alnA = alignments[alignments["qid"] == qid].iloc[0]
        alnB = alignments[alignments["qid"] == hapid].iloc[0]
    '''