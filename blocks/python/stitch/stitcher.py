import sys
sys.path.append("..")
import log.log as logger

import stitch.plot_intervals as plotter
import stitch.build_intervals as builder
import copy

OUTPUT_BLOCK_NAME = "blocks.bed"
OUTPUT_MBLOCK_NAME = "megablocks.bed"

def interval2bed(intervals, fileName):
    
    for interval in intervals:
        rgb =  logger.rbg_generator(interval.rid)
        rInfo = str(interval.rid) + ":" + str(interval.left(q=False)) + "-" + str(interval.right(q=False))
        qInfo = [interval.qid, interval.left(q=True), interval.right(q=True), 
                 rInfo, interval.percent_identity(), interval.get_dir(q=False,string=True),
                interval.left(q=True), interval.right(q=True), rgb]
        
        logger.FileLogger().write_cols(fileName, qInfo)


def stitch(tigId, alignDict, lengthData, param, q=True):    
    """Performs all stitching tasks: creating chunks, blocks and megablocks.
    tigId is the ID of the assembled contig to stitch.
    alignDict is the dictionary containing chunk alignments between assemblies.
    lengthData is a dictionary containing the lengths of all assembled contigs.
    param is the Parameters object.
    q is True if tigId is a query contig, False if reference contig.
    
    q=True implies chunks, blocks, etc will be formed along query contig.
    q=False is only intended to plot from the perspective of the reference.

    If param.PLOT_BLOCKS is not None, Contig object will be plotted.

    Returns a Contig object that contains all constructed megablocks.
    Returns None if no chunks alignments exist."""

    plotDir = param.PLOT_BLOCKS
    
    if tigId not in alignDict:
        logger.out("Contig " + str(tigId) + " has no chunk alignments.", 1, param)
        return None
    
    rows = alignDict[tigId]

    #save Intervals for plotting
    cList, bList, mList = [], [], []

    mblockList=[]
    length = lengthData[tigId]
    
    for idx, row in rows.iterrows():
        #if row[2] == 'tig00007745_pilon_pilon': break

        '''
        idx=0
        for i in range(idx): row = rows.iterrows().send(None)        
        '''       
        #chunks
        chunks = builder.construct_chunks(row, param)
        if plotDir is not None: cList.extend(copy.deepcopy(chunks))
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        if plotDir is not None: bList.extend(copy.deepcopy(blocks))
        
        #megablock
        megablocks = builder.construct_megablocks(blocks, length, param)
            
        mblockList.extend(megablocks)
        if plotDir is not None: mList.extend(copy.deepcopy(megablocks))

    #contig
    contig = builder.construct_contig(mblockList, tigId, length, param)    
    
    #plot
    if len(cList) > 0 and plotDir is not None:
        
        tList = [] if contig is None else contig.mblocks
        if not q: tList = []
        intervalList=[cList, bList, mList, tList]
        plotter.plot_levels(intervalList, length, q=q, outputPath=plotDir + tigId)

    return contig


def stitch_blocks(tigId, alignDict, param):    
    if tigId not in alignDict:
        logger.out("Contig " + str(tigId) + " has no chunk alignments.", 1, param)
        return []
    
    rows = alignDict[tigId]
    blockList = []
    for idx, row in rows.iterrows():
        
        #chunks
        chunks = builder.construct_chunks(row, param)
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        blockList.extend(blocks)

    interval2bed(blockList, OUTPUT_BLOCK_NAME)


    return blockList


def stitch_megablocks(tigId, blocks, lengthData, param):    
    
    mblockList=[]
    length = lengthData[tigId]
    rids = set([block.rid for block in blocks])
    
    for rid in rids:
        rblocks = [block for block in blocks if block.rid == rid]
        megablocks = builder.construct_megablocks(rblocks, length, param)
        mblockList.extend(megablocks)

    return mblockList


def stitch_contig(tigId, mblocks, lengthData, param):    
    
    length = lengthData[tigId]
    contig = builder.construct_contig(mblocks, tigId, length, param)    

    if contig is not None:
        interval2bed(contig.mblocks, OUTPUT_MBLOCK_NAME)


    return contig


        