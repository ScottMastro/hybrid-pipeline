import os
import sys
import parameters as param
import utils.log as logger

import utils.file_handler as io
import utils.fasta_handler as fasta

import stitch.stitcher as stitcher
import weld.welder as welder
import scaffold.scaffolder as scaffolder

def main(param):

    parser = param.set_hybrid_parameters()
    parameters = param.get_hybrid_parameters(parser)

    if not param.validate_hybrid_parameters(parameters):
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(param.OUTPUT_DIR): os.mkdir(param.OUTPUT_DIR)

    logger.Logger(clean=True, level=param.VERBOSE, wait=param.WAIT)
    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)

    #==================================================
    # Read in data
    #==================================================
    print("ref=", param.REF_FA)
    print("query=", param.QUERY_FA)
    print("summary=", param.SUMMARY)

    logger.log("Reading reference fasta...")
    refData = fasta.read_fasta(param.REF_FA)
    logger.log("Reading query fasta...")
    queryData = fasta.read_fasta(param.QUERY_FA)
    
    seqData = dict()
    seqData.update(refData) ; seqData.update(queryData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    
    logger.log("Reading alignment data...")
    alignDict = io.parse_alignments(param.SUMMARY)

    rids, qids = list(refData.keys()), list(queryData.keys())
    if not fasta.validate_ids(rids, qids): sys.exit()
    
    rids.sort(key=lambda x: -lengthData[x])
    qids.sort(key=lambda x: -lengthData[x])

    confBed = io.parse_bed(param.REF_BED)   

    #==================================================
    # Stitch alignments
    #==================================================

    contigs = io.unpickle(param.OUTPUT_DIR + "/contigs.pickle")
    if contigs is not None:
        logger.log("Found existing checkpoint pickle. Loading stitched data...")
    else:
        logger.log("Stitching contigs...")
        
        contigs = []
        for tigId in qids:
            contig = stitcher.stitch(tigId, alignDict, lengthData, param)
            if contig is not None: contigs.append(contig)

        logger.FileLogger().flush_all()
        io.pickle(contigs, param.OUTPUT_DIR + "/contigs.pickle")

    #todo: output information about BLOCK ORDER! (missassembly detection)

    #==================================================
    # Weld blocks
    #==================================================
    
    paths = io.unpickle(param.OUTPUT_DIR + "/paths.pickle")
    if paths is not None:    
        logger.log("Found existing checkpoint pickle. Loading welded data...")
    else:
        logger.log("Welding contigs...")

        paths = []
        for contig in contigs:
            path = welder.weld(contig, seqData, lengthData, param)
            if path is not None: paths.append(path)

        logger.FileLogger().flush_all()
        io.pickle(paths, param.OUTPUT_DIR + "/paths.pickle")        
        
    #==================================================
    # Scaffold paths
    #==================================================

    rawScaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds_raw.pickle")
    if rawScaffolds is not None:    
        logger.log("Found existing checkpoint pickle. Loading scaffolded data...")
    else:
        logger.log("Scaffolding contigs...")

        ridSet = set()
        for path in paths:
            for fork in path: ridSet.add(fork.rid)
        
        for tigId in rids:
            if tigId not in ridSet: continue
            paths, exclude = scaffolder.scaffold(paths, tigId, confBed, lengthData, param)

        rawScaffolds = paths
        logger.FileLogger().flush_all()
        io.pickle(rawScaffolds, param.OUTPUT_DIR + "/scaffolds_raw.pickle")        

    #==================================================
    # Salvage leftovers
    #==================================================

    scaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds.pickle")
    if scaffolds is not None:    
        logger.log("Found existing checkpoint pickle. Loading salvaged data...")
    else:
        logger.log("Salvaging contigs...")

        scaffolds, unusedDict = scaffolder.salvage(rawScaffolds, qids, rids, lengthData, minSize=10000)   
        io.pickle(scaffolds, param.OUTPUT_DIR + "/scaffolds.pickle")

    #==================================================
    # Write output
    #==================================================

    logger.log("Writing sequence...")

    fasta.write_hybrid(scaffolds, seqData, param)
    
    leftoverRegions = []
    for key in unusedDict:
        leftoverRegions.extend(unusedDict[key])
    fasta.write_leftover(leftoverRegions, seqData, param,
                      minSize=10000, compressNs=10, spacerNs=500)
    
    logger.FileLogger().flush_all()
