import sys
sys.path.append("./analysis")
import log as logger
import parameters
import file_handler as io
import weld.build_path as welder
import weld.path_helper as path_helper
import stitch.stitcher as stitcher
import scaffold.scaffolder as scaffolder
import excluded_regions as salvager

OUTPUT_BLOCK_NAME = "blocks.bed"
OUTPUT_MBLOCK_ALL_NAME = "mblocks_all.bed"
OUTPUT_MBLOCK_NAME = "mblocks.bed"
OUTPUT_PATH_NAME = "paths.txt"
OUTPUT_SCAFFOLD_NAME = "scaffolds.txt"
OUTPUT_LEFTOVERS_NAME = "leftovers.bed"


def main():

    #--------------------------------------
    # READ IN DATA
    #--------------------------------------

    print("Parsing parameters...")
    param = parameters.get_parameters()

    print("Reading Canu fasta...")
    refData = io.read_fasta(param.REF_FA)
    print("Reading Supernova fasta...")
    queryData = io.read_fasta(param.QUERY_FA)
    
    seqData = dict()
    seqData.update(refData)
    seqData.update(queryData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    
    print("Reading alignment data...")
    aligndf = io.parse_alignments(param.SUMMARY)
    alignDict = {str(tigId):rows for (tigId,rows) in tuple(aligndf.groupby('qid'))}

    rids, qids = list(refData.keys()), list(queryData.keys())
    rids.sort(key=lambda x: -lengthData[x])
    qids.sort(key=lambda x: -lengthData[x])

    #if not io.validate_ids(rids, qids): return
    #todo: validate output files here first!

    confBed = io.parse_confident_regions(param.REF_BED, rids, param)   
    
    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)

    #--------------------------------------
    # STITCH BLOCKS
    #--------------------------------------

    try:
        contigs = io.unpickle(param.OUTPUT_DIR + "/contigs.pickle")
        print("Loading stitched data...")

    except:
        print("Stitching contigs...")

        blockId, mblockId = 1,1
        blockPrefix, mblockPrefix = "block", "mblock"

        def add_id(intervals, prefix, i):
            for interval in intervals:
                interval.iid = prefix + str(i)
                i += 1
            return i
        
        def write_bed(intervals, fileName):
            logger.interval2bed(intervals, "q_" + fileName, q=True)
            logger.interval2bed(intervals, "r_" + fileName, q=False)

        contigs = []
        for tigId in qids:
            logger.out("Doing contig: " + tigId, 1, param)
    
            blocks = stitcher.stitch_blocks(tigId, alignDict, param)
            blockId = add_id(blocks, blockPrefix, blockId)
            write_bed(blocks, OUTPUT_BLOCK_NAME)
            
            mblocks = stitcher.stitch_megablocks(tigId, blocks, lengthData, param)
            mblockIdwrite_scaffold = add_id(mblocks, mblockPrefix, mblockId)
            write_bed(mblocks, OUTPUT_MBLOCK_ALL_NAME)
            
            contig = stitcher.stitch_contig(tigId, mblocks, lengthData, param)
            if contig is not None:
                write_bed(contig.mblocks, OUTPUT_MBLOCK_NAME)
                contigs.append(contig)
        
        logger.FileLogger().flush_all()
        #io.pickle(contigs, param.OUTPUT_DIR + "/contigs.pickle")

    #--------------------------------------
    # WELD BLOCKS
    #--------------------------------------

    try:
        paths = io.unpickle(param.OUTPUT_DIR + "/paths.pickle")
        print("Loading welded data...")
    except:
        print("Welding contigs...")
        
        pathId, pathPrefix = 1, "path"

        def write_path(pid, megaPath, terminatingForks, scaffoldForks, fileName):
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

        paths = []
        for contig in contigs:
            if contig is None: continue
            megaPaths = []
            pid = pathPrefix + str(pathId)
            terminatingForks = []
            for megablock in reversed(contig.mblocks):
                #path = welder.weld_megablock(megablock, seqData, lengthData, param)
                
                blockPaths = welder.weld_megablock(megablock, seqData, param)
                terminatingForks.extend([path[-1] for path in blockPaths])
                megaPath = welder.join_blockpaths(blockPaths, lengthData, param)
                megaPaths.append(megaPath)

            scaffoldForks = [mpath[-2] for mpath in megaPaths[:-1]]
            path = welder.join_megablockpaths(megaPaths, lengthData, param)
            path = path_helper.clean_NNNs(path)
            if len(path) < 1: continue
            path.set_path_id(pid)    
            paths.append(path)
            
            write_path(pid, path, terminatingForks, scaffoldForks, OUTPUT_PATH_NAME)
            pathId += 1

        logger.FileLogger().flush_all()
        #io.pickle(paths, param.OUTPUT_DIR + "/paths.pickle")        


    #--------------------------------------
    # SCAFFOLD BLOCKS
    #--------------------------------------
    try:
        scaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds_raw.pickle")
        print("Loading scaffolded data...")
    except:
        print("Scaffolding contigs...")

        ridSet = set()
        for path in paths:
            for fork in path: ridSet.add(fork.rid)
    
        excludeList = []
    
        for tigId in rids:
            if tigId not in ridSet: continue
            logger.out("Scaffolding with contig: " + tigId, 1, param)
    
            paths, exclude = scaffolder.scaffold(paths, tigId, confBed, lengthData, param)
            excludeList.extend(exclude)
            #if stop: input()
            print("npaths:",len(paths))

        scaffolds = paths
                    
        io.pickle(scaffolds, param.OUTPUT_DIR + "/scaffolds_raw.pickle")        

    #--------------------------------------
    # SALVAGE LEFTOVERS
    #--------------------------------------
    
    try:
        scaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds.pickle")
        print("Loading salvaged data...")
    except:
        print("Salvaging contigs...")

        unusedQuery = salvager.get_unused_regions(scaffolds, qids, lengthData)
        unusedRef = salvager.get_unused_regions(scaffolds, rids, lengthData)
        unusedDict = unusedQuery ; unusedDict.update(unusedRef)
        
        scaffolds, unusedDict = salvager.extend_path_ends(scaffolds, unusedDict, lengthData)
    
        scaffoldId, scaffoldPrefix = 1, "scaff"
        def write_scaffold(sid, scaffold, fileName):
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
        for scaffold in scaffolds:
            sid = scaffoldPrefix + str(scaffoldId)
            write_scaffold(sid, scaffold, OUTPUT_SCAFFOLD_NAME)
            scaffold.pid=sid
            scaffoldId += 1
    
        logger.FileLogger().flush_all()

    
        for tigId in unusedDict:
            lid=1
            for region in unusedDict[tigId]:
                s,e = region.start, region.end
                info = [tigId, s, e, str(tigId) + "_" + str(lid), 
                        round(100*abs(e-s)/lengthData[tigId],2)]
                lid += 1
                logger.FileLogger().write_cols(OUTPUT_LEFTOVERS_NAME, info)

        logger.FileLogger().flush_all()
        
        io.pickle(scaffolds, param.OUTPUT_DIR + "/scaffolds.pickle")     


    #--------------------------------------
    # OUTPUT
    #--------------------------------------

    print("Writing sequence...")

    logger.FileLogger().flush_all()
    io.write_hybrid(scaffolds, seqData, param)
    
    leftoverRegions = []
    for key in unusedDict:
        leftoverRegions.extend(unusedDict[key])
    io.write_leftover(leftoverRegions, seqData, param,
                      minSize=10000, compressNs=10, spacerNs=500)






    '''
            
    mblockId=1
    mblockPrefix="mblock"

    mid = mblockPrefix+str(mblockId)
    output_path(megaPath, mid, seqData, lengthData, param)    
    
    mblockId+=1
    '''
            
            #analyzer.analyze_blockpaths(blockPaths, seqData, lengthData, param)
            
            #for i,blockPath in enumerate(blockPaths):
            #    print(blockPath)
                
                #polisher.polish_block(blockPath, seqData, lengthData, param)
                #if i+1 < len(blockPaths):
                #    blockPathBefore, blockPathAfter = blockPath, blockPaths[i+1]
                #    polisher.polish_interblock(blockPathBefore, blockPathAfter, seqData, lengthData, param)

            #input()
            #try:
            #    polisher.polish_block(megaPath, seqData, lengthData, param)
            #except:
            #    pass
            
            

    


#debugging help
def p(tigId, pos, rc=False, dist=10000, length=2000):
    print("length = " + str(lengthData[tigId]))
    print("before:\n")
    if not rc:
        print(seqData[tigId][pos-dist-int(length/2):pos-dist+int(length/2)])
    else:
        print(seqData[tigId][lengthData[tigId]-pos-dist-int(length/2):lengthData[tigId]-pos-dist+int(length/2)])

    print("after:\n")
    if not rc:
        print(seqData[tigId][pos+dist-int(length/2):pos+dist+int(length/2)])
    else:
        print(seqData[tigId][lengthData[tigId]-pos+dist-int(length/2):lengthData[tigId]-pos+dist+int(length/2)])

def pp(tigId, pos, rc=False, length=2000):
    print("length = " + str(lengthData[tigId]))
    if not rc:
        print(seqData[tigId][pos-int(length/2):pos+int(length/2)])
    else:
        print(seqData[tigId][lengthData[tigId]-pos-int(length/2):lengthData[tigId]-pos+int(length/2)])



if __name__== "__main__":
  main()
  print("done")
  exit()