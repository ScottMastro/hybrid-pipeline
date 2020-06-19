import sys
sys.path.append("../..")

import utils.log as logger
import structures.path as pathtools
import hybridtools.scaffold.build_scaffold as builder
import hybridtools.scaffold.excluded_regions as salvager

'''
 The scaffolder takes welded paths, whose walk represents
 hybrid contigs, and attempts to connect paths together at their ends.
 
 The salvager component attempts to pick up any remaining paths
 that possibly represent unique regions that may have been lost
 during the scaffolding process.
'''

#==================================================
# Logging helpers
#==================================================

OUTPUT_SCAFFOLD_NAME = "scaffolds.txt"
OUTPUT_LEFTOVERS_NAME = "leftovers.bed"

def write_scaffold(sid, scaffold, fileName, lengthData):
    for i,fork in enumerate(scaffold[:-1]):
        frm = fork.after_pos_norm(lengthData)
        to = scaffold[i+1].before_pos_norm(lengthData)
        chrPos = str(fork.after_id()) + ":" + \
                 str(frm) + "-" +  str(to) 
        strand = fork.after_strand()
        if strand > 0: strand = "+" 
        elif strand < 0: strand = "-"
        info = [sid, chrPos, strand, scaffold.pid]
        logger.FileLogger().write_cols(fileName, info)

#==================================================
# Main scaffold function
#==================================================

def scaffold(paths, tigId, unitigs, lengthData, param):
    '''
    Takes a tigId and a list of paths. Attempts to scaffold every Path that
    contains a Fork with tigId on the end.
    Returns a tuple of:
    1) a list of scaffolded Paths
    2) a list of Paths that were discarded (ex. redundantly overlaps with scaffold)
    '''
    
    logger.log("Scaffolding with contig: " + tigId, logger.LOG_DETAILS)
    relevantPaths, irrelevantPaths = pathtools.filter_paths(paths, tigId, endsOnly=False)    
    scaffolds, excludePaths = builder.scaffold_tigs(relevantPaths, tigId, unitigs, lengthData, param)
    
    return (irrelevantPaths + scaffolds, excludePaths)

#==================================================
# Main salvage function
#==================================================

def salvage(scaffolds, qids, rids, lengthData, minSize=10000):
        
    logger.log("Salvaging unused segments.", logger.LOG_DETAILS)

    unusedQuery = salvager.get_unused_regions(scaffolds, qids, lengthData, minSize)
    unusedRef = salvager.get_unused_regions(scaffolds, rids, lengthData, minSize)
    unusedDict = unusedQuery ; unusedDict.update(unusedRef)
    
    logger.log("Salvaging the ends of scaffolds.", logger.LOG_DETAILS)
    scaffolds, unusedDict = salvager.extend_path_ends(scaffolds, unusedDict, lengthData)

    scaffoldId, scaffoldPrefix = 1, "scaff"
    for scaffold in scaffolds:
        sid = scaffoldPrefix + str(scaffoldId)
        write_scaffold(sid, scaffold, OUTPUT_SCAFFOLD_NAME, lengthData)
        scaffold.pid=sid
        scaffoldId += 1

    for tigId in unusedDict:
        lid=1
        for region in unusedDict[tigId]:
            s,e = region.start, region.end
            info = [tigId, s, e, str(tigId) + "_" + str(lid), 
                    round(100*abs(e-s)/lengthData[tigId],2)]
            lid += 1
            
            logger.log("Salvaged: " + str(region), logger.LOG_DEBUG, indent=1)

            logger.FileLogger().write_cols(OUTPUT_LEFTOVERS_NAME, info)

    logger.FileLogger().flush_all()

    return scaffolds, unusedDict