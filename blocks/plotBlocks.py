import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
    
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
