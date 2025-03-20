import sys
sys.path.append("../..")

import utils.log as logger
import stitch.build_intervals as builder

'''
 The stitcher takes BLAST alignments between two assemblies
 and creates a hierarchy of data structures representing
 these alignments at different scales of the genome.
 
 The main structures are defined as Interval objects and 
 have the following hierarchy:
    
      chunks < blocks < megablocks < contigs
'''

#==================================================
# Stitching subfunctions
#==================================================

def stitch_blocks(tigId, alignDict, param):
    """
    Extracts in alignments (chunks) from alignDict. Contiguous chunks produce
    block intervals. Returns a list of blocks for contig tigId.
    """
    if tigId not in alignDict:
        logger.log("Contig " + str(tigId) + " has no alignments.", logger.LOG_DETAILS)
        return []
    
    rows = alignDict[tigId]
    blockList = []
    for idx, row in rows.iterrows():
        
        #chunks
        chunks = builder.construct_chunks(row, param)
        
        #blocks
        blocks, trash = builder.construct_blocks(chunks, param)
        blockList.extend(blocks)

    return blockList


def stitch_megablocks(tigId, blocks, lengthData, param):    
    """
    Takes blocks formed by contiguous chunks and creates larger blocks
    with more tolerance of divergence (megablocks).
    Returns a list of megablocks.
    """

    mblockList=[]
    length = lengthData[tigId]
    rids = set([block.rid for block in blocks])
    
    for rid in rids:
        rblocks = [block for block in blocks if block.rid == rid]
        megablocks = builder.construct_megablocks(rblocks, length, param)
        mblockList.extend(megablocks)

    return mblockList


def stitch_contig(tigId, mblocks, lengthData, param):    
    """
    Takes megablocks and creates a non-overlapping set spanning a contig.
    Returns a list of megablocks.
    """

    length = lengthData[tigId]
    return builder.construct_contig(mblocks, tigId, length, param)

#==================================================
# Logging helpers
#==================================================
    
def interval2bed(intervals, fileName, q=True):
    for interval in intervals:
        tigId = str(interval.qid if q else interval.rid)
        rgb = logger.rbg_generator(tigId)
        name = str(interval.iid) + "_id=" + str(round(interval.percent_identity(), 3)) + \
                "_cov=" +  str(interval.total_coverage()) 
        info = [tigId, interval.left(q), interval.right(q), name,
                interval.percent_identity(), interval.get_dir(q, string=True),
                interval.left(q), interval.right(q), rgb]
                
        logger.FileLogger().write_cols(fileName, info)

def write_bed(intervals, fileName):
    interval2bed(intervals, "q_" + fileName, q=True)
    interval2bed(intervals, "r_" + fileName, q=False)

OUTPUT_BLOCK_NAME = "blocks.bed"
OUTPUT_MBLOCK_ALL_NAME = "mblocks_all.bed"
OUTPUT_MBLOCK_NAME = "mblocks.bed"
blockId, mblockId = 1,1
blockPrefix, mblockPrefix = "block", "mblock"

#==================================================
# Main stitch function
#==================================================

def stitch(tigId, alignDict, lengthData, param):
    """
    Performs all stitching operations. Outputs block information as BED files.
    Returns a list of contigs (contig = list of megablocks).
    """
    tigInfo = tigId + " (" + str(lengthData[tigId]) + " bp)"
    logger.log("Stitching contig: " + tigInfo, logger.LOG_DETAILS)

    global blockId, mblockId
    
    def add_id(intervals, prefix, i):
        for interval in intervals:
            interval.iid = prefix + str(i)
            i += 1
        return i

    # chunks and blocks
    blocks = stitch_blocks(tigId, alignDict, param)
    write_bed(blocks, OUTPUT_BLOCK_NAME)
    
    # megablocks
    mblocks = stitch_megablocks(tigId, blocks, lengthData, param)
    blockId = add_id(blocks, blockPrefix, blockId)

    mblockId = add_id(mblocks, mblockPrefix, mblockId)
    write_bed(mblocks, OUTPUT_MBLOCK_ALL_NAME)
    
    # contig
    contig = stitch_contig(tigId, mblocks, lengthData, param)
    if contig is None or len(contig) < 1: return None
    write_bed(contig, OUTPUT_MBLOCK_NAME)
    
    logger.log(str(len(contig)) + " Megablocks formed: ", logger.LOG_DEBUG)
    for mblock in contig: 
        logger.log(str(mblock.iid) + "\t" + mblock.log_string(), logger.LOG_DEBUG, indent=1)
    
    return contig
