import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

def construct_rect(intervals, col_generator, q=True):
    
    if len(intervals) < 1:
        return None
        
    intervals.sort(key=lambda interval: interval.left(q))

    y=1
    ycoord = [-1 for i in intervals]
    lastInterval = None
    done = False

    while not done:
        done = True
        for i in range(len(intervals)):
            if(ycoord[i] < 0):
                done = False
                if(lastInterval is None or lastInterval.right(q) <= intervals[i].left(q)):
                    lastInterval = intervals[i]
                    ycoord[i] = y
        if y==20:
            for i in range(len(intervals)):
                if(ycoord[i] < 0): ycoord[i] = 21
                
        y = y + 1
        lastInterval = None      
                
    xcoordLeft = [interval.left(q)/1000 for interval in intervals]
    xcoordRight = [interval.right(q)/1000 for interval in intervals]
    fills = [col_generator(interval.rid if q else interval.qid) \
                           for interval in intervals]
    
    rects = [Rectangle( (xl, y), xr-xl, 1) for xl,xr,y in zip(xcoordLeft, xcoordRight, ycoord)]
    
    return PatchCollection(rects, facecolor=np.array(fills), edgecolor=col_generator(None))


def get_ymax(collection):
    return [path.get_extents().max[1] for path in collection.get_paths()]
def get_coords(path):
    x1=path.get_extents().min[0] 
    x2=path.get_extents().max[0] 
    y1=path.get_extents().min[1] 
    y2=path.get_extents().max[1] 
    return (x1, x2, y1, y2)


def plot_levels(intervalList, length, q=True, outputPath=None, step=500):
            
    #remove extension
    if outputPath is not None:
        outputPath = os.path.splitext(outputPath)[0]
    if step is None: step = length+1
    
    chunkPlotWidth = 0.05
    chunkPlotHeight = chunkPlotWidth*5
    maxLevels = 15
    
    def col_generator(tigId): 
        if tigId is None: return "black"
        np.random.seed((length + hash(tigId) % 2**30))
        return [np.random.rand(), np.random.rand(), np.random.rand()]

    rectangles = [construct_rect(intervals, col_generator, q) for intervals in intervalList]

    ymax = [min(maxLevels, 0 if rects is None else max(get_ymax(rects)) + 1) \
            for rects in rectangles]

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(ncols=1, nrows=4, figure=fig, height_ratios=[y + 0.5 for y in ymax])
    
    axs = []
    for i in range(len(intervalList)):
        
        ax = fig.add_subplot(gs[i, 0])
        if rectangles[i] is not None:
            ax.add_collection(rectangles[i])
        
        ax.set_ylim(0.5, ymax[i])
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)    
        axs.append(ax)
        
    blockidx = 1
    if rectangles[blockidx] is not None:
        for block, path in zip(intervalList[blockidx], rectangles[blockidx].get_paths()):
            x1, x2, y1, y2 = get_coords(path)
            h = abs(y1-y2)/4

            label = str(round(block.start(not q)/1000.0))
            axs[blockidx].annotate(label, xy=(x1, y1 + h), clip_on=True, ha='left', va='center',
                         path_effects=[PathEffects.withStroke(linewidth=2.5,foreground="w")])
                
            label = str(round(block.end(not q)/1000.0))
            axs[blockidx].annotate(label, xy=(x2, y2 - h), clip_on=True, ha='right', va='center',
                         path_effects=[PathEffects.withStroke(linewidth=2.5,foreground="w")])
        
    mblockidx = 2
    if rectangles[mblockidx]  is not None:
        for mblock, path in zip(intervalList[mblockidx], rectangles[mblockidx].get_paths()):
            x1, x2, y1, y2 = get_coords(path)
            
            space=150
            label = "(" + (mblock.rid if q else mblock.qid) + ")"

            if x2-x1 < space*2:
                axs[mblockidx].annotate(label, xy=((x1+x2)/2, (y1+y2)/2), clip_on=True, va='center', ha='center',
                             path_effects=[PathEffects.withStroke(linewidth=2.5,foreground="w")])
            else:
                for i in np.arange(x1, x2, space):
                    x = (2*i+space)/2
                    if x > x2: continue
                    axs[mblockidx].annotate(label, xy=(x, (y1+y2)/2), clip_on=True, va='center', ha='center',
                                 path_effects=[PathEffects.withStroke(linewidth=2.5,foreground="w")])


    for i in range(0, int(length/1000), step):

        part = int(i/step)
        minX = i
        maxX = min(length/1000, i+step)
        
        for ax in axs:
            ax.set_xlim(minX, maxX)

        plt.xticks(np.arange(minX, maxX, 10))
        fig.set_size_inches(max(3, (maxX - minX)*chunkPlotWidth), sum(ymax)*chunkPlotHeight)
        
        if outputPath is not None:
            fig.savefig(outputPath + "_" + str(part) + ".png")
                
    plt.close('all')