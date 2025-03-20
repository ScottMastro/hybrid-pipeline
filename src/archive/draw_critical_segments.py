import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def plot_segments(criticalForkList, tigId, lengthData, q=False):
    plt.switch_backend('Qt5Agg')
    plt.close("all")
    # Figure and Axes
    fig, ax = plt.subplots(1,1,figsize=(14,14), facecolor='#f7f7f7', dpi= 80)
    tigLength = lengthData[tigId]

    indexDict, i = {}, 0
    tigIds = []
    for cf in criticalForkList:
        tig = cf.path[0].get_id(not q)
        tigIds.append(tig)
        if tig in indexDict: continue
        else:
            indexDict[tig] = i
            ax.hlines(y=i, xmin=0, xmax=tigLength, color='black', alpha=0.05, linewidth=1)
            i = i + 1
                           
    height = i+1
    # Vertical Lines
    #ax.vlines(x=.2*tigLength, ymin=0, ymax=height, color='black', alpha=0.4, linewidth=1, linestyles='dotted')
    #ax.vlines(x=.4*tigLength, ymin=0, ymax=height, color='black', alpha=0.4, linewidth=1, linestyles='dotted')
    #ax.vlines(x=.6*tigLength, ymin=0, ymax=height, color='black', alpha=0.4, linewidth=1, linestyles='dotted')
    #ax.vlines(x=.8*tigLength, ymin=0, ymax=height, color='black', alpha=0.4, linewidth=1, linestyles='dotted')
    
    # Func to draw line segment
    def newline(p1, p2, color='black'):
        #ax = plt.gca()
        l = mlines.Line2D([p1[0],p2[0]], [p1[1],p2[1]], color='skyblue', linewidth=3)
        ax.add_line(l)
        return l
                

    starts = [cf.minPos for cf in criticalForkList]
    ends = [cf.maxPos for cf in criticalForkList]
    index = [indexDict[tig] for tig in tigIds]
    
    for x in starts + ends:
        ax.vlines(x=x, ymin=0, ymax=height, color='black', alpha=0.6, linewidth=1, linestyles='dotted')
    
    # Line Segments
    for p1, p2, tig in zip(starts, ends, tigIds):
        i = indexDict[tig]
        ax.annotate(tig, xy=((p1+p2)/2,i-0.1),  xycoords='data',
                    ha='center', va='top')
        newline([p1, i], [p2, i])

    # Points
    ax.scatter(y=index, x=starts, s=15, color='#0e668b', alpha=0.7)
    ax.scatter(y=index, x=ends, s=15, color='#a3c4dc', alpha=0.7)
    
    # Decoration
    ax.set_facecolor('#f7f7f7')
    ax.set_title("Alignments to " + tigId, fontdict={'size':22})
    ax.set(xlim=(0, tigLength), ylim=(-0.5, height), ylabel='index')
    
    fig.set_size_inches(10, 3)

    plt.show()
    plt.savefig("./scaffold_plots/" + tigId + ".png")
    
    
    