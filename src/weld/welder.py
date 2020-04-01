import sys
sys.path.append("..")
import utils.log as logger

import weld.build_path as builder

# LOGGING HELPERS
# -----------------

def write_path(pid, megaPath, terminatingForks, scaffoldForks, fileName, lengthData):
    for i,fork in enumerate(megaPath[:-1]):
        chrPos = str(fork.after_id()) + ":" + \
                 str(fork.after_pos_norm(lengthData)) + "-" + \
                 str(megaPath[i+1].before_pos_norm(lengthData))
        strand = fork.after_strand()
        if strand > 0: strand = "+" 
        elif strand < 0: strand = "-"
        regionType = "gap" if fork.is_switch_reference() else ""
        if fork in terminatingForks:
            regionType = "megagap"
        if fork in scaffoldForks:
            regionType = "scaffold"
        info = [pid, chrPos, strand, regionType]
        logger.FileLogger().write_cols(fileName, info)

OUTPUT_PATH_NAME = "paths.txt"
pathId = 1 

# -----------------

def weld(contig, seqData, lengthData, param):
    '''
    Welds together a contig into a Path. Blocks in megablocks have their ends
    anchored between query and reference. NNNs in query sequence are replaced 
    using reference sequence. A single Path is returned.
    '''
    if contig is None or len(contig) == 0 : return None

    pathPrefix = "path"
    global pathId

    megaPaths = []
    terminatingForks = []
    for megablock in reversed(contig):
        #path = welder.weld_megablock(megablock, seqData, lengthData, param)
        
        blockPaths = builder.weld_megablock(megablock, seqData, param)
        terminatingForks.extend([path[-1] for path in blockPaths])
        
        megaPath = builder.join_blockpaths(blockPaths, lengthData, param)
        megaPaths.append(megaPath)

    scaffoldForks = [mpath[-2] for mpath in megaPaths[:-1]]
    path = builder.join_mblockpaths(megaPaths, lengthData, param)
    
    if len(path) < 1: return None
    
    pid = pathPrefix + str(pathId)
    path.set_path_id(pid)
    write_path(pid, path, terminatingForks, scaffoldForks, OUTPUT_PATH_NAME, lengthData)
    pathId += 1
      
    return path
