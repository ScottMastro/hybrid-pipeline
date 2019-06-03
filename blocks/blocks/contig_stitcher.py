import numpy as np
import log
import block_builder as blocker
import plot_blocks as plotter
import copy
from intervals import Contig

def stitch(aligndf, tigId, param):    
    
    #QUERY ONLY
    #dataframe containing only alignments for one contig
    df = aligndf.loc[aligndf[str(aligndf.columns[0])] == tigId]
    if len(df) == 0:
        log.out("Contig " + str(tigId) + " has no chunk alignments!", 1, param)
        return None
    
    mblockList=[]
    
    #row=next(df.iterrows())[1]
    for idx, row in df.iterrows():
        
        #if row[2] == 'tig00032236_pilon_pilon':
        #    break
        
        nchunks = int(row[1])
        length = nchunks*param.CHUNK_SIZE
        
        chunks = blocker.construct_chunks(row, param)
        blocks, trash = blocker.construct_blocks(chunks, param)
        blocks = blocker.remove_overlap(blocks, param)            
        mblocks = blocker.construct_megablocks2(blocks, length, param)
        mblockList.extend(mblocks)

    mblockList = blocker.clean_megablocks(mblockList, param)

    contig = Contig(tigId, length, mblockList)

    return contig

def stitch_r(aligndf, rname, param, mblockOnly=False):    
    
    #REF ONLY
    #dataframe containing only alignments for one contig
    df = aligndf.loc[aligndf[str(aligndf.columns[2])] == rname]
    if len(df) == 0:
        log.out("Contig " + str(rname) + " has no chunk alignments!", 1, param)
        return None
    
    mblockList=[]
    
    #row=next(df.iterrows())[1]
    for idx, row in df.iterrows():
        length = row[3]
        
        chunks = blocker.construct_chunks(row, param)
        blocks, trash = blocker.construct_blocks(chunks, param)
        mblocks = blocker.construct_megablocks2(blocks, length, param)
        mblockList.extend(mblocks)

    if mblockOnly:
        return(mblockList)

    contig = blocker.construct_contig(rname, length, mblockList, param, q=False)
    return contig

def plot(aligndf, tigId, param, q=True):

    if q:
        df = aligndf.loc[aligndf[str(aligndf.columns[0])] == tigId]
    else:
        df = aligndf.loc[aligndf[str(aligndf.columns[2])] == tigId]

    blockList=[]
    trashList=[]
    mblockList=[]
    mblockListCopy=[]
    
    #row=next(qdf.iterrows())[1]
    for idx, row in df.iterrows():
        
        nchunks = int(row[1])
        length = nchunks*param.CHUNK_SIZE if q else int(row[3])
        
        chunks = blocker.construct_chunks(row, param)
        blocks, trash = blocker.construct_blocks(chunks, param)

        blockList.append(copy.deepcopy(blocks))
        trashList.append(copy.deepcopy(trash))
        
        blocks = blocker.remove_overlap(blocks, param)            

        mblocks = blocker.construct_megablocks2(blocks, length, param)
        mblockList.append(copy.deepcopy(mblocks))
        mblockListCopy.extend(mblocks)

    mblockListCopy = blocker.clean_megablocks(mblockListCopy, param)

    contigList = []
    for mblocks in mblockList:
        if len(mblocks) < 1:
            contigList.append([])
            continue
        tig = mblocks[0].rid if q else mblocks[0].qid
        lst = []
        for mblocksCopy in mblockListCopy:
            tigCopy = mblocksCopy[0].rid if q else mblocksCopy[0].qid
            if tig == tigCopy:
                lst.append(mblocksCopy)
                
        contigList.append(lst)
        
    print('Generating plot')
    if q:
        plotter.plot_levels(blockList, trashList, mblockList, contigList, nchunks,
                            step=500,  \
                            #outputPath="/home/scott/Dropbox/hybrid-pipeline/blocks/plots/"  \
                            outputPath="/media/scott/Rotom/assembly_data/plots/" \
                            #outputPath="/media/scott/HDD/sickkids/NA24385/plots/" \
                            + str(tigId))
    else:
        plotter.plot_levels_r(blockList, trashList, mblockList, contigList, length,
                            step=500, 
                            #outputPath="/home/scott/Dropbox/hybrid-pipeline/blocks/plotscanu/" \
                            outputPath="/media/scott/Rotom/assembly_data/plots/" \
                            #outputPath="/media/scott/HDD/sickkids/NA24385/plots/" \
                            + str(tigId))

def plot_contigs(aligndf, param):

    qnames=np.unique(aligndf[str(aligndf.columns[0])])
    rnames=np.unique(aligndf[str(aligndf.columns[2])])
    q=True
        
    for tigId in qnames if q else rnames:

        print("Plotting " + str(tigId))
        plot(aligndf, tigId, param, q)
        