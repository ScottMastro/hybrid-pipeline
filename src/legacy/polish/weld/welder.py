import weld.build_path as builder
import weld.path_helper as path_helper
from weld.paths import Path
 
def weld(contig, seqData, lengthData, param):
    '''
    Welds together a Contig into a Path. Blocks in Megablocks have their ends
    anchored between query and reference. NNNs in query sequence are replaced 
    using reference sequence. A single Path is returned.
    '''
    if contig is None: return Path()
    
    plotDir = param.DOT_PLOT
    contigPaths = []
      
    for megablock in reversed(contig.mblocks):
        blockPaths = builder.weld_megablock(megablock, seqData, param)
        if plotDir is not None: 
            import weld.dotplot as dotplot
            dotplot.gaps_dotplot(blockPaths, seqData, lengthData, plotDir)
        
        megaPath = builder.join_blockpaths(blockPaths, lengthData, param)
        
        if megaPath is not None and len(megaPath) > 0:
            contigPaths.append(megaPath)
     
    path = builder.join_megablockpaths(contigPaths, lengthData, param)
        
    ok = path_helper.check_path(path)
    
    if not ok:
        print("Not ok")
        input()
        
    return path