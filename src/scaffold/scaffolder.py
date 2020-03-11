import weld.path_helper as path_helper
import scaffold.build_scaffold as helper
import scaffold.critical_forks as critical_forks
from scaffold.scaffold_codes import fail_enum as f


def remove_overlapping_paths(paths, lengthData):
        
    paths.sort(key=lambda path: path_helper.path_length(path))
    lengths = [path_helper.path_length(path) for path in paths]
    
    keep = []
    
    for i,path1 in enumerate(paths):
        shouldKeep=True
        print(i)
        #if (path_helper.path_length(scaffolds[i]) < 50000): continue
        #if(i<340): continue
        for j,path2 in enumerate(reversed(paths)):
            
            if lengths[j] >= lengths[i] and i != j:
                o1, o2 = path_helper.path_overlap(path1, path2, lengthData, source='r')
                if max(o1, o2) > 0.99:
                    shouldKeep=False
                    break

        if shouldKeep: keep.append(path1)

def filter_by_size(paths, size):
    
    longPaths = [ path for path in paths if path_helper.path_length(path) > size ]
    shortPaths = [ path for path in paths if path_helper.path_length(path) <= size ]
    return(longPaths, shortPaths)


def scaffold(paths, tigId, unitigs, lengthData, param):
    '''
    Takes a tigId and a list of paths. Attempts to scaffold every Path that
    contains a Fork with tigId on the end.
    Returns a tuple of:
    1) a list of scaffolded Paths or None
    2) a list of Paths that do not have tigId on the end
    3) a list of the leftovers Paths that failed to scaffold
    '''
    relevantPaths, irrelevantPaths = path_helper.filter_paths(paths, tigId, ends=False)    
    
    if len(relevantPaths) < 2: 
        print("Nothing to scaffold")
        return (paths, [])

    #build critical forks and sort left to right
    criticalForkList, excludePaths = [], []
    for path in relevantPaths: 
        criticalForkList.extend( critical_forks.build(path, tigId, lengthData) )

    criticalForkList.sort(key=lambda cf: cf.minPos)

    #plot critical segments
    #plotter.plot_segments(criticalForkList, tigId, lengthData)
    #return (paths, [])

    i=0
    results = []
    while i < len(criticalForkList) -1:
        print(criticalForkList[i].minTigId, criticalForkList[i+1].minTigId) 
        result = helper.can_scaffold(criticalForkList, i, unitigs, lengthData, param)
        resultCode, resultHint = result
        
        #one segment fully overlaps another, remove the smaller one
        if resultCode == f.FULL_OVERLAP:
            excludePaths.append(criticalForkList[resultHint].path)
            criticalForkList.pop(resultHint)
            if resultHint == i and i > 0:
                results.pop()
                i=i-1
            continue
            
        #segments are coming from the same path with nothing inbetween, merge
        if resultCode == f.SAME_PATH:
            criticalForkList[i].merge(criticalForkList[i+1])
            criticalForkList.pop(i+1)
            continue
        
        if resultCode == f.UNSURE:
            print("!!!!!")
            #todo: does this happen?
            input()
        
        results.append(result)
        i = i + 1
        
    scaffoldCandidates = [cf.path for cf in criticalForkList]

    print([f.failDict[resultCode] for resultCode,_ in results])
    print([cf.minTigId for cf in criticalForkList])

    #need to create this so that a path is not duplicated if it has
    #multiple critical forks
    redundantDict = dict()
    for i,pth in enumerate(scaffoldCandidates):
        for j,otherpth in enumerate(scaffoldCandidates):
            if i < j and pth == otherpth: 
                redundantDict[j] = i 
                print("WARNING: WATCH OUT FOR REDUNDANCY!")
    
    for i in redundantDict.keys(): scaffoldCandidates[i] = None
    #scaffoldCopy = copy.deepcopy(scaffoldCandidates)

    i=0
    redoChance = True
    while i < len(scaffoldCandidates) -1:
        redoChance = True

        if scaffoldCandidates[i] is None:
            scaffoldCandidates[i] = scaffoldCandidates[redundantDict[i]]
            scaffoldCandidates[redundantDict[i]] = None
        if scaffoldCandidates[i+1] is None:
            scaffoldCandidates[i+1] = scaffoldCandidates[redundantDict[i+1]]
            scaffoldCandidates[redundantDict[i+1]] = None

        resultCode, resultHint = results[i]
        if resultCode == f.PASS:
            fork1, fork2 = resultHint
            path1, path2 = scaffoldCandidates[i], scaffoldCandidates[i+1]
            if path1 is not None and path2 is not None:
                scaffold = helper.scaffold_pair(path1, path2, fork1, fork2, lengthData, param)
    
                #try again
                if scaffold is None and redoChance:
                    cf1 = critical_forks.build(scaffoldCandidates[i], tigId, lengthData)
                    cf2 = critical_forks.build(scaffoldCandidates[i+1], tigId, lengthData)
                    result = helper.can_scaffold([cf1[0],cf2[0]], 0, unitigs, lengthData, param)
                    results[i] = result
                    redoChance = False
                    continue
                if scaffold is not None:
                    scaffoldCandidates[i] = None
                    scaffoldCandidates[i+1] = scaffold
            
        i = i+1
    
    scaffolds = [s for s in scaffoldCandidates if s is not None]
    print(len(scaffolds), " scaffolds returned.")
    

    paths = irrelevantPaths + scaffolds
    
    return (paths, excludePaths)


    #NOTES:
    # tigId = rids[4]
    # tig00031394_pilon_pilon has missassembly at end w/small contigs
    # require min contig length for scaffolding?

    # tigId = rids[8]
    # tig00031671_pilon_pilon same story on chr1, small contigs at end

    # tigId = rids[7]
    # tig00000130_pilon_pilon
    # ends in ambiguous region on chr10, small contig 1051 matches
    # we probably don't want to include it

    # tigId = rids[15]
    # tig00031818_pilon_pilon
    # seems to be a very subtle misassembly in canu
    # tig00031818_pilon_pilon, 588157 matches contig 1077
    # tig00031818_pilon_pilon, 610644 matches contig 318 (22487 bp dist)
    # but in HG38, aligns to (107916331, 107918330) and (106226108, 106228105)
    # (1688226 bp dist), checked 318 at the potential scaffold region looks good

    # tigId = rids[19]
    # tig00031382_pilon_pilon
    # very small overlap with 248 but it seems to produce a productive scaffold!

    # tigId = rids[20]
    # tig00031577_pilon_pilon
    # beautiful job at recognizing same-chrom misassembly in canu

    # tigId = rids[24]
    # tig00031412_pilon_pilon
    # **********two ugly supernovas at the end (small), keep canu and remove sn???

    # tigId = rids[45]
    # tig00002045_pilon_pilon
    #same as above, we should throw out small sn contigs at end!
    
    
    # tigId = rids[50]
    # tig00031958_pilon_pilon
    # continue peel example
    #same as above, we should throw out small sn contigs at end!

    # tigId = rids[109]
    # tig01067927_pilon_pilon
    # incorrect scaffold since there is no excess supernova...
    # not much can be done? polish to fix???

    # tigId = rids[118]
    # tig00031616_pilon_pilon
    # peel back example

    # tigId = rids[122]
    # tig00032199_pilon_pilon
    # weird one, 3 misassemblies???

    # tig00031438_pilon_pilon misassembly in the middle???