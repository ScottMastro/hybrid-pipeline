import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import os

def makeRectangle(left, right, stack_factor, height, gap, chunk=True,):

    if(chunk):
        y=-(stack_factor*(gap + height) + gap + height)
    else:
        y=(stack_factor*(gap + height) + gap)
    
    return Rectangle( (left, y), right-left, height)

def makeAnnotation(left, i, stack_factor, height, gap, chunk=True,):
   
    if(chunk):
        y=-(stack_factor*(gap + height) + gap + height) + height/8
    else:
        y=(stack_factor*(gap + height) + gap) - height/8
    
    return [i, left + height/12, y]


def plotBlocks(row_df, block_summary, height=100, gap=20):
    
    row_df=pd.DataFrame.reset_index(row_df, drop=True)
    #block_summary=block_summary.sort_values(by=['left'])
    block_summary=pd.DataFrame.reset_index(block_summary,drop=True)

    global prev_right
    prev_right=0
    global stack_factor
    stack_factor=1
    global max_stack_factor
    max_stack_factor=1

    def stackRect(left, right):
        global prev_right
        global stack_factor
        global max_stack_factor
        
        if(left < prev_right):
            stack_factor=stack_factor+1
            if(stack_factor > max_stack_factor):
                max_stack_factor = stack_factor
        else:
            stack_factor=1
            
        prev_right=right
        return stack_factor
        
    #chunks-----------------------------------
    left_chunk=row_df["left"]
    right_chunk=row_df["right"]
    chunk_id=row_df["chunk"]

    prev_right=0
    stack_factor=1
    chunks = [makeRectangle(l, r, stackRect(l, r), height, gap)
                for l,r in zip(left_chunk,right_chunk)]
    prev_right=0
    stack_factor=1
    chunk_label = [makeAnnotation(l, i, stackRect(l, r), height, gap)
                for l,r,i in zip(left_chunk,right_chunk,chunk_id)]
    
    #blocks-----------------------------------

    left_block=block_summary["left"]
    right_block=block_summary["right"]

    prev_right=0
    stack_factor=1
    blocks = [makeRectangle(l, r, stackRect(l, r), height, gap, False)
                for l,r in zip(left_block,right_block)]
    
    
    visible_range=20050
    
    #plot-----------------------------------

    fig, ax = plt.subplots(1)
    
    pc_block = PatchCollection(blocks, facecolor='b', alpha=0.5, edgecolor='b')
    ax.add_collection(pc_block)

    pc_chunk = PatchCollection(chunks, facecolor='r', alpha=0.5, edgecolor='r')
    ax.add_collection(pc_chunk)
    
    plt.ylim(-1600, 800)
    fig.set_size_inches(12, 5)

    for index,block in block_summary.iterrows():
        
        left=block["left"]
        chrom=block["chrom"]
        contig=block["contig"]
                
        #left
        x_limit=(left - visible_range/2, left + visible_range/2)
        plt.xlim(x_limit)

        for label, x, y in chunk_label:
            if x >= x_limit[0] and x <= x_limit[1]:
                plt.annotate(str(label), xy=(x, y), clip_on=True)
        
        plt.title(str(chrom) + " - " + str(contig) + " - block " + str(index))

        fig.savefig("/media/scott/Rotom/plots/" + str(chrom) + "_" + str(contig) + "_" + str(index) + "_left.png")

        
    for index,block in block_summary.iterrows():
        
        right=block["right"]
        chrom=block["chrom"]
        contig=block["contig"]
        
        #right
        x_limit=(right - visible_range/2, right + visible_range/2)
        plt.xlim(x_limit)

        for label, x, y in chunk_label:
            if x >= x_limit[0] and x <= x_limit[1]:
                plt.annotate(str(label), xy=(x, y), clip_on=True)
        
        plt.title(str(chrom) + " - " + str(contig) + " - block " + str(index))
        fig.savefig("/media/scott/Rotom/plots/" + str(chrom) + "_" + str(contig) + "_" + str(index) + "_right.png")
    plt.close('all')
    
    
def plotQuery(block_summary, query, show_label=False, height=100, gap=20):
    
    summary=block_summary.loc[block_summary["chrom"]==query]
    summary=summary.sort_values(by=['left'])

    global prev_right
    prev_right=0
    global stack_factor
    stack_factor=1
    global max_stack_factor
    max_stack_factor=1

    def stackRect(left, right):
        global prev_right
        global stack_factor
        global max_stack_factor
        
        if(left < prev_right):
            stack_factor=stack_factor+1
            if(stack_factor > max_stack_factor):
                max_stack_factor = stack_factor
        else:
            stack_factor=1
            
        prev_right=right
        return stack_factor
        
    #query-----------------------------------
    left=0
    right=summary["span"].iloc[0]

    query_rect = [makeRectangle(left, right, 0, height, gap)]
    
    #blocks-----------------------------------

    left_block=summary["left"]
    right_block=summary["right"]

    prev_right=0
    stack_factor=1
    blocks = [makeRectangle(l, r, stackRect(l, r), height, gap, False)
                for l,r in zip(left_block,right_block)]
    
    if(show_label):
        block_id=summary["block"]
        prev_right=0
        stack_factor=1
        block_label = [makeAnnotation(l, i, stackRect(l, r), height, gap)
                    for l,r,i in zip(left_block,right_block,block_id)]
            
    #plot-----------------------------------

    fig, ax = plt.subplots(1)
    
    pc_block = PatchCollection(query_rect, facecolor='b', alpha=0.5, edgecolor='b')
    ax.add_collection(pc_block)

    pc_chunk = PatchCollection(blocks, facecolor='r', alpha=0.5, edgecolor='r')
    ax.add_collection(pc_chunk)
    
    plt.ylim(-200, 1800)
    plt.xlim(left, right)

    fig.set_size_inches(12, 5)

    if(show_label):
        for label, x, y in block_label:
            plt.annotate(str(label), xy=(x, y), clip_on=True)

    plt.title(str(query) + " - blocks")
    fig.savefig("/media/scott/Rotom/plots/" + str(query) + ".png")
    plt.show()
    plt.close('all')

















import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

def draw_rect(intervals, edge):
    
    if len(intervals) < 1:
        return None
        
    intervals.sort(key=lambda interval: interval.left_id())


    y=1
    ycoord = [-1 for i in intervals]
    lastInterval = None
    done = False

    while not done:
        done = True
        for i in range(len(intervals)):
            if(ycoord[i] < 0):
                done = False
                if(lastInterval is None or intervals[i].left_id() > lastInterval.right_id()):
                    lastInterval = intervals[i]
                    ycoord[i] = y
                    
        y = y + 1
        lastInterval = None      
        
    xcoordLeft = [interval.left_id() for interval in intervals]
    xcoordRight = [interval.right_id() for interval in intervals]
    fills = [interval.fill for interval in intervals]
    
    rects = [Rectangle( (xl-0.5, y), xr-xl+1, 1) for xl,xr,y in zip(xcoordLeft, xcoordRight, ycoord)]
    
    return PatchCollection(rects, facecolor=np.array(fills), edgecolor=edge)

def get_ymin(collection):
    return [path.get_extents().min[1] for path in collection.get_paths()]
def get_xmin(collection):
    return [path.get_extents().min[0] for path in collection.get_paths()]
def get_ymax(collection):
    return [path.get_extents().max[1] for path in collection.get_paths()]
def get_xmax(collection):
    return [path.get_extents().max[0] for path in collection.get_paths()]
def get_coords(path):
    x1=path.get_extents().min[0] 
    x2=path.get_extents().max[0] 
    y1=path.get_extents().min[1] 
    y2=path.get_extents().max[1] 
    return (x1, x2, y1, y2)

def plot_levels(blockList, trashList, mblockList, contigList, \
                nchunks, step=None, outputPath=None, annotate=True):
            
    #remove extension
    if outputPath is not None:
        outputPath = os.path.splitext(outputPath)[0]
    if step is None:
        step = nchunks+1
    
    chunkPlotWidth = 0.05
    chunkPlotHeight = chunkPlotWidth*5
    
    chunks = []
    blocks = []
    mblocks = []
    contig = []

    np.random.seed(nchunks)
    
    for blockSublist, trashSublist, mblockSublist, contigSubList in \
        zip(blockList, trashList, mblockList, contigList):
            
        col = [np.random.rand(), np.random.rand(), np.random.rand()]
        
        for block in blockSublist:
            block.fill = col
            blocks.append(block)
            for chunk in block.components:
                chunk.fill = col
                chunks.append(chunk)
        for trash in trashSublist:
            for chunk in trash.components:
                chunk.fill = col
                chunks.append(chunk)
        for mblock in mblockSublist:
            mblock.fill = col
            mblocks.append(mblock)
        for mblock in contigSubList:
            mblock.fill = col
            contig.append(mblock)
    
    chunkRects = draw_rect(chunks, "black")
    blockRects = draw_rect(blocks, "black")
    mblockRects = draw_rect(mblocks, "black")
    contigRects = draw_rect(contig, "black")

    ymaxChunk = 0 if chunkRects is None else max(get_ymax(chunkRects)) + 1
    ymaxBlock = 0 if blockRects is None else max(get_ymax(blockRects)) + 1
    ymaxMblock = 0 if mblockRects is None else max(get_ymax(mblockRects)) + 1
    ymaxContig = 0 if contigRects is None else  max(get_ymax(contigRects)) + 1

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(ncols=1, nrows=4, figure=fig, \
                           height_ratios=[ymaxChunk+0.5, 
                                          ymaxBlock+0.5, 
                                          ymaxMblock+0.5,
                                          ymaxContig+0.5])
    
    ax0 = fig.add_subplot(gs[0, 0])
    if chunkRects is not None:
        ax0.add_collection(chunkRects)
        
    ax0.set_ylim(0.5, ymaxChunk)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['top'].set_visible(False)
    
    #---------------------------------------
    
    ax1 = fig.add_subplot(gs[1, 0])
    if blockRects is not None:
        ax1.add_collection(blockRects)
    
    ax1.set_ylim(0.5, ymaxBlock)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)    
    
    #---------------------------------------
    
    ax2 = fig.add_subplot(gs[2, 0])
    if mblockRects is not None:
        ax2.add_collection(mblockRects)
    
    ax2.set_ylim(0.5, ymaxMblock)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    
    #---------------------------------------
    
    ax3 = fig.add_subplot(gs[3, 0])
    if contigRects is not None:
        ax3.add_collection(contigRects)
    
    ax3.set_ylim(0.5, ymaxContig)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    
    if annotate:
        if blockRects is not None:
            for block, path in zip(blocks, blockRects.get_paths()):
                x1, x2, y1, y2 = get_coords(path)
                
                label = str(round(block.start(q=False)/1000.0))
                ax1.annotate(label, xy=(x1, y1), clip_on=True, ha='center',
                             path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
                    
                label = str(round(block.end(q=False)/1000.0))
                ax1.annotate(label, xy=(x2, y2), clip_on=True, ha='center',
                             path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
            
        if mblockRects is not None:
            for mblock, path in zip(mblocks, mblockRects.get_paths()):
                x1, x2, y1, y2 = get_coords(path)
                
                space=150
        
                if x2-x1 < space*2:
                    label = str(mblock.rid)
                    ax2.annotate(label, xy=((x1+x2)/2, (y1+y2)/2), clip_on=True, va='center', ha='center',
                                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
                else:
                    for i in np.arange(x1, x2, space):
                        label = str(mblock.rid)
                        ax2.annotate(label, xy=((2*i+space)/2, (y1+y2)/2), clip_on=True, va='center', ha='center',
                                     path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])

    for i in range(0, nchunks, step):

        part = int(i/step)
        minX = i   
        maxX = min(nchunks, i+step)
        
        ax0.set_xlim(minX, maxX)
        ax1.set_xlim(minX, maxX)
        ax2.set_xlim(minX, maxX)
        ax3.set_xlim(minX, maxX)

        plt.xticks(np.arange(minX, maxX, 10))
        fig.set_size_inches(max(3, (maxX - minX)*chunkPlotWidth), (ymaxChunk+ymaxBlock+ymaxMblock)*chunkPlotHeight)
        
        if outputPath is not None:
            fig.savefig(outputPath + "_" + str(part) + ".png")
        
    plt.close('all')
    
