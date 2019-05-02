import numpy as np
import log
import block_builder as blocker
import plot_blocks as plotter
from gfa_handler import GFA

def stitch(aligndf, qname, param):    
    
    #QUERY ONLY TODO
    #dataframe containing only alignments for one contig
    df = aligndf.loc[aligndf[str(aligndf.columns[0])] == qname]
    if len(df) == 0:
        log.out("Contig " + str(qname) + " has no chunk alignments!", 1, param)
        return None
    
    mblockList=[]
    
    #row=next(df.iterrows())[1]
    for idx, row in df.iterrows():
        
        nchunks = int(row[1])
        length = nchunks*param.CHUNK_SIZE
        
        chunks = blocker.construct_chunks(row, param)
        blocks, trash = blocker.construct_blocks(chunks, param)
        blocks = blocker.remove_redundancy(blocks, param)            
        mblocks = blocker.construct_megablocks(blocks, length, param)
        mblockList.extend(mblocks)

    contig = blocker.construct_contig(qname, length, mblockList, param, q=True)
    return contig

def plot(aligndf, qname, param):

    df = aligndf.loc[aligndf[str(aligndf.columns[0])] == qname]
    
    blockList=[]
    trashList=[]
    mblockList=[]
    
    #row=next(qdf.iterrows())[1]
    for idx, row in df.iterrows():
        
        nchunks = int(row[1])
        length = nchunks*param.CHUNK_SIZE
        
        chunks = blocker.construct_chunks(row, param)
        blocks, trash = blocker.construct_blocks(chunks, param)
        
        blockList.append(blocks)
        trashList.append(trash)

        blocks = blocker.remove_redundancy(blocks, param)            
        mblocks = blocker.construct_megablocks(blocks, length, param)

        mblockList.append(mblocks)


    mblocks = []
    for mblockSublist in mblockList:
        mblocks.extend(mblockSublist)
    contig = blocker.construct_contig(qname, length, mblocks, param, q=True)

    contigList = []
    for mblockSublist in mblockList:
        newList = []
        for mblock in mblockSublist:
            if mblock in contig.mblocks:
                newList.append(mblock)
        contigList.append(newList)

    plotter.plot_levels(blockList, trashList, mblockList, contigList, nchunks,
                        step=500, outputPath="./plots/" + str(qname) )

def stitch_contigs(aligndf, param):

    #gfa = GFA()

    qnames=np.unique(aligndf[str(aligndf.columns[0])])
    contigs=dict()
    for qname in qnames:
        #print("Plotting " + str(qname))
        #plot(aligndf, qname, param)
        
        contigs[qname] = stitch(aligndf, qname, param)
        
        #gfa = add_contig(contigs[-1], gfa, isQuery=True)

    #gfa.write("hello.gfa")

    return contigs

    
'''
import sys
import file_handler as fh

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
else:
    #prefix="C:/Users/scott/Dropbox/hybrid-pipeline/blocks"
    prefix="/home/scott/Dropbox/hybrid-pipeline/blocks"
    summary_file = prefix + "/summary.txt"
    block_file = prefix + "/blocks_new.txt"

index = 0
aligndf = fh.parse_alignments(summary_file)   
'''