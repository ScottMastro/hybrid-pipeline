import sys
sys.path.append("./analysis")
import log.log as logger
import parameters
import file_handler as io
import weld.build_path as welder
import weld.path_helper as path_helper
import stitch.stitcher as stitcher
import scaffold.scaffolder as scaffolder
import excluded_regions as salvager
import analysis.analyze_sequence as analyzer

    
def main():

    print("Parsing parameters...")
    param = parameters.get_parameters()
    writer = logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)

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

    confBed = io.parse_confident_regions(param.REF_BED)    
    
    #--------------------------------------
    #--------------------------------------
    #--------------------------------------

    writer = logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)

    try:
        contigs = io.unpickle(param.OUTPUT_DIR + "/contigs.pickle")
        print("Loading stitched contigs...")

    except:
        print("Stitching contigs...")
        
        contigs = []
        for tigId in qids:
            logger.out("Doing contig: " + tigId, 1, param)
    
            blocks = stitcher.stitch_blocks(tigId, alignDict, param)
            mblocks = stitcher.stitch_megablocks(tigId, blocks, lengthData, param)
            contig = stitcher.stitch_contig(tigId, mblocks, lengthData, param)
        
            if contig is not None:  
                contigs.append(contig)
        
        io.pickle(contigs, param.OUTPUT_DIR + "/contigs.pickle")
    
    print("Welding contigs...")
    paths, emptyIds = [], []
    for contig in contigs:
        megapaths = []       
        for megablock in reversed(contig.mblocks):
            #path = welder.weld_megablock(megablock, seqData, lengthData, param)
            
            blockPaths = welder.weld_megablock(megablock, seqData, param)
            analyzer.analyze_blockpaths(blockPaths, seqData, param)

            megaPath = welder.join_blockpaths(blockPaths, lengthData, param)
            megapaths.append(megaPath)
            
        paths.append(megapaths)
      
    #emptyIds.append(tigId) if len(path) < 1 else       
            
    unusedQuery = salvager.get_unused_regions(paths, qids, lengthData)


    #todo:
    #output megapath bridges
    #output N filling + flank seqence
    


    paths = [ path_helper.clean_NNNs(path) for path in paths ]
    paths = [ path for path in paths if len(path) > 0 ]

    io.pickle(paths, "paths.pickle")

    '''
    paths = io.unpickle("CF062_paths.pickle")
    paths = io.unpickle("hg002_paths10.pickle")
    paths = io.unpickle("hg002_scaffolds.pickle")

    '''

    #--------------------------------------
    #--------------------------------------
    #--------------------------------------
    
    #paths = io.unpickle("hg002_paths10.pickle")
    #paths = io.unpickle("hg002_scaffolds.pickle")

    paths = [ path_helper.clean_NNNs(path) for path in paths ]
    paths = [ path for path in paths if len(path) > 0 ]

    ridSet = set()
    for path in paths:
        for fork in path: ridSet.add(fork.rid)

    print("Iterating over reference contigs...")
    scaffolds, excludeList = [], []

    for tigId in rids:
        if tigId not in ridSet: continue
        logger.out("Scaffolding with contig: " + tigId, 1, param)

        paths, exclude = scaffolder.scaffold(paths, tigId, confBed, lengthData, param)
        excludeList.extend(exclude)
        print("npaths:",len(paths))
        #if stop: input()
        
    scaffolds = paths
    #io.pickle(scaffolds, "scaffolds.pickle")

    #plt.hist([ np.log10(path.path_length()) for path in paths ])
    #plt.show()
    #print(sum([ path.path_length() for path in paths ]))

    #--------------------------------------
    #--------------------------------------
    #--------------------------------------
    
    unusedRef = salvager.get_unused_regions(paths, rids, lengthData)
    unusedQuery = salvager.get_unused_regions(paths, qids, lengthData)

        
    
    writer.flush_all()
    io.write_hybrid(scaffolds, seqData, param)


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