import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import random
import matplotlib.patches as patches

def get_spline_nodes(xStart, xEnd, length, direction=0, yLimit=(10,1500), chaosFactor=0.5, xStretch=2):
    
    yMid = 0
    xGap = xEnd- xStart
    yMin, yMax = yLimit
    
    if length / xGap < 3:
        xStretch = 0
    #length = max(length, xGap*2)
    
    degree = int(np.log2(length)*chaosFactor)
    if degree < 1: degree = 1
        
    def make_nodes():
         
        yHeight = min(yMax*direction, length*direction, key=abs)
        if direction == 0: yHeight = 2*min(yMax, length, key=abs)
        
        nodes = [ [x, yHeight*random.random() - (yHeight/2 if direction == 0 else 0)] \
                   for x in np.linspace(xStart, xEnd, degree + 1, endpoint=False) ]
    
        nodes = nodes[1:]
        xPos = [node[0] for node in nodes]
    
        for node in nodes:
            node[0] = (1-chaosFactor)*node[0] + \
                (chaosFactor)*(random.uniform(xStart - xGap*xStretch, xEnd + xGap*xStretch))
            
        nodes[0][1] = nodes[0][1] + (0 - nodes[0][1])*(1-chaosFactor)
        nodes[0][-1] = nodes[0][-1] + (0 - nodes[0][-1])*(1-chaosFactor)

        for i in range(len(nodes)-1):
            nodes[i][1] = nodes[i][1] + (nodes[i+1][1] -  nodes[i][1])*(1-chaosFactor)
            nodes[i][0] = nodes[i][0] + (nodes[i+1][0] -  nodes[i][0])*(1-chaosFactor)
        

        return nodes, xPos

    
    nodes, xPos = make_nodes()
    lengthCheck = False
    lastIteration = None
    while True:
        
        if direction == 0:
            x = [xStart, xStart+xGap/1000] + [x[0] for x in nodes] + [xEnd - +xGap/1000, xEnd]
            y = [yMid, yMid] + [y[1] for y in nodes] + [yMid, yMid]
        else:
            x = [xStart, xStart] + [x[0] for x in nodes] + [xEnd, xEnd]
            y = [yMid, yMid + direction] + [y[1] for y in nodes] + [yMid + direction, yMid]

        tck,u     = interpolate.splprep( [x,y] ,s = 0 )
        xSpline, ySpline = interpolate.splev( np.linspace( 0, 1, 8*degree ), tck,der = 0)
    
        total = 0
        for i in range(len(xSpline[:-1])):
            total = total + np.sqrt( (xSpline[i+1] - xSpline[i])**2 + (ySpline[i+1] - ySpline[i])**2 )
            
            #ySum = ySum + abs(abs(ySpline[i+1]) - abs(xSpline[i]))
            #xSum = xSum + abs(abs(xSpline[i+1]) - abs(xSpline[i]))

        #print( str(total - length))
        
        lenDiff = total - length 

        if not lengthCheck and lenDiff < -max(xGap, length/100):
            degree = int(degree*1.5)
            nodes, xPos = make_nodes()
            #print("adding nodes")
            continue
        else:
            lengthCheck = True
        
        #print(lenDiff)
        
        if lastIteration is not None and abs(lastIteration - lenDiff) <= 1:
            break
        if lenDiff < length/100: break
    
        lastIteration = total - length
        
        for i, node in enumerate(nodes):
            node[0] = node[0] + (xPos[i] - node[0])/16
            node[1] = max(yMin * direction, node[1] - (node[1])/16, key=abs)

    if direction == 0:
        nodes = np.array( [ [xStart, yMid], [xStart+xGap/1000, yMid] ] + \
                              nodes + \
                          [ [xEnd - xGap/1000, yMid], [xEnd, yMid] ] )
    else:
        nodes = np.array( [ [xStart, yMid], [xStart, yMid + direction] ] + \
                              nodes + \
                          [ [xEnd, yMid + direction], [xEnd, yMid] ] )
       
    return nodes

def draw_spline(nodes, colour, ax, alpha=1, precision=1000):
    x = nodes[:,0]
    y = nodes[:,1]
    
    tck,u = interpolate.splprep([x,y], s=0)
    xSpline, ySpline = interpolate.splev(np.linspace(0, 1, precision), tck, der=0)
    
    total = 0
    for i in range(len(xSpline)-1):
        total = total + np.sqrt( (xSpline[i+1] - xSpline[i])**2 + (ySpline[i+1] - ySpline[i])**2 )

    print("length=" + str(total))
    ax.plot(xSpline, ySpline, color=colour, alpha=alpha, linewidth=0.5)
    #ax.plot( x, y, 'o', color="orange")

def col_generator(tigId, lengthData): 
    if tigId is None: return "black"
    np.random.seed((lengthData[tigId] + hash(tigId) % 2**30))
    return [np.random.rand(), np.random.rand(), np.random.rand()]


    
def plot_block_path2(plotPaths, lengthData):

    def normalized_pos(fork, q=True):

        pos = fork.qpos if q else fork.rpos
        if pos is None: return None
        if (fork.qstrand if q else fork.rstrand) == -1:
            pos = lengthData[(fork.qid if q else fork.rid)] - pos
        return pos
    
    xGap=2000
    alpha = 0.8
    
    fig, ax = plt.subplots()
    xPos = 0
    y = 10


    for megaPath in plotPaths[:6]:       
        
        for i, blockPath in enumerate(megaPath):
            rColour = col_generator(megaPath[0][0].rid, lengthData)
            qColour = col_generator(megaPath[0][0].qid, lengthData)
            #ax.plot([blockPath[0].qpos,blockPath[-1].qpos], [y,y], \
             #color=qColour, alpha=alpha)
    
            ax.plot([blockPath[0].qpos,blockPath[-1].qpos], [0,0], \
             color="black", alpha=alpha, lw=0.1)
    
            xMin = min(normalized_pos(blockPath[0]), normalized_pos(blockPath[-1]))
            xMax = max(normalized_pos(blockPath[0]), normalized_pos(blockPath[-1]))
            
            
            if blockPath[0].qstrand == -1: x1, x2 = xMax, xMin
            else: x1, x2 = xMin, xMax
            qArrow = patches.FancyArrowPatch((x1, y), (x2, y), color=qColour, alpha=0.6, mutation_scale=20, arrowstyle="wedge")
            qArrow2 = patches.FancyArrowPatch((x1, 0), (x2, 0), color=qColour, alpha=0.6, mutation_scale=20, arrowstyle="wedge")

            if blockPath[0].rstrand == -1: x1, x2 = xMax, xMin
            else: x1, x2 = xMin, xMax
            rArrow = patches.FancyArrowPatch((x1, -y), (x2, -y), color=rColour, alpha=0.6, mutation_scale=20, arrowstyle="wedge")
            
            if i+1 < len(megaPath):
                xMin = xMax
                xMax = min(normalized_pos(megaPath[i+1][-1]), normalized_pos(megaPath[i+1][0]))
                
                if blockPath[-1].rstrand == -1: x1, x2 = xMax, xMin
                else: x1, x2 = xMin, xMax
                rArrow2 = patches.FancyArrowPatch((x1, 0), (x2, 0), color=rColour, alpha=0.6, mutation_scale=20, arrowstyle="wedge")
                
                ax.add_patch(rArrow2)
                
            ax.add_patch(qArrow2)
            ax.add_patch(qArrow)
            ax.add_patch(rArrow)
    
    
        xMin = min(normalized_pos(megaPath[0][0]), normalized_pos(megaPath[-1][-1]))
        xMax = max(normalized_pos(megaPath[0][0]), normalized_pos(megaPath[-1][-1]))
    
        rect = patches.Rectangle( (xMin, -y*4), xMax-xMin, y*8,
            fill=False, color="black", clip_on=False)
        ax.add_patch(rect)
    
    plt.axis( [ None, None, -100, 100] )
    fig.set_size_inches(28.5, 8.5)

    plt.show()

        
    '''
        xPairs = []
        tigIds = []
        
        xPairsAlt = []
        tigIdsAlt = []
    
        for i, path in enumerate(mpath[:-1]):
            xPairs.append([path.after_pos(), mpath[i+1].before_pos()])
            tigIds.append(path.after_id())
            xPairsAlt.append([path.before_pos(), mpath[i+1].after_pos()])
            tigIdsAlt.append(path.before_id())
    
        
    
        for i, ((x1, x2), (x1Alt, x2Alt)) in enumerate(zip(xPairs, xPairsAlt)):
            colour = col_generator(tigIds[i], lengthData)
            colourAlt = col_generator(tigIdsAlt[i], lengthData)
    
            xDist = abs(x2 - x1)
            xDistAlt = abs(x2Alt - x1Alt)
            print("expected length=" + str(xDistAlt))
            ax.plot([xPos, xPos+xDist], [yPos, yPos], color=colour, alpha=alpha, lw=1)
            if xDist < xGap:
                ax.scatter([(xPos+xPos+xDist)/2], [yPos], color=colour, alpha=alpha)
            else:
                nodes = get_spline_nodes(xPos, xPos + xDist, xDistAlt, \
                                         direction=1, yLimit=(4000,10000), chaosFactor=0.7, xStretch=0)
                draw_spline(nodes, colourAlt, ax, precision=1000)
    
    
            xPos = xPos + xDist + xGap
        
    plt.axis( [ -1000, None, -10000, 10000] )
  
    plt.show()

'''
    
    
def plot_block_path3(plotPaths, lengthData):

    def normalized_pos(fork, q=True):
        pos = fork.qpos if q else fork.rpos
        if pos is None: return None
        if (fork.qstrand if q else fork.rstrand) == -1:
            pos = lengthData[(fork.qid if q else fork.rid)] - pos
        return pos

    forkGap=10000
    smallForkGap=10000
    bigForkGap=50000
    
    alpha = 0.8
    fig, ax = plt.subplots()
    nextForkPos = 0
    yLimit = (10, 1500)
    for megaPath in plotPaths[:8]:       
    
        for i, blockPath in enumerate(megaPath):
            
            # draw alternative path
            
            lengthAlt = abs(blockPath[-1].after_pos() - blockPath[0].before_pos())
            length = abs(blockPath[-1].before_pos() - blockPath[0].after_pos())
            
            fg = forkGap
            if length > bigForkGap:
                fg = bigForkGap
            elif length < forkGap:
                fg = smallForkGap
            
            col = col_generator(blockPath[0].before_id(), lengthData)
            nodes = get_spline_nodes(nextForkPos, nextForkPos+fg, lengthAlt, \
                                     direction=1, yLimit=(500, 3500), chaosFactor=0.9, xStretch=0)
            draw_spline(nodes, col, ax, alpha=0.3, precision=1000)

            # draw block
            col = col_generator(blockPath[0].after_id(), lengthData)
            nodes = get_spline_nodes(nextForkPos, nextForkPos+fg, length, \
                                     direction=0, yLimit=yLimit, chaosFactor=0.9, xStretch=0)
            draw_spline(nodes, col, ax, precision=1000)
            plt.scatter(nextForkPos, 0, marker="D", color="black", zorder=10, s=5)


            nextForkPos = nextForkPos + fg

            # draw between block
            if i+1 < len(megaPath):

                lengthAlt = abs(megaPath[i+1][0].after_pos() - blockPath[-1].before_pos())
                length = abs(megaPath[i+1][0].before_pos() - blockPath[-1].after_pos())
                
                fg = forkGap
                if length > bigForkGap:
                    fg = bigForkGap
                elif length < forkGap:
                    fg = smallForkGap
                
                # draw alternative path
                col = col_generator(blockPath[-1].before_id(), lengthData)
                nodes = get_spline_nodes(nextForkPos, nextForkPos+fg, length, \
                                         direction=-1, yLimit=(500, 3500), chaosFactor=0.9, xStretch=0)
                draw_spline(nodes, col, ax, alpha=0.3, precision=1000)
            
    
                col = col_generator(blockPath[-1].after_id(), lengthData)
    
                nodes = get_spline_nodes(nextForkPos, nextForkPos+fg, length, \
                                         direction=0, yLimit=yLimit, chaosFactor=0.5, xStretch=0)
                draw_spline(nodes, col, ax, precision=1000)
    
                plt.scatter(nextForkPos, 0, marker="D", color="black", zorder=10, s=5)
                nextForkPos = nextForkPos + fg


                    
    #plt.axis( [ None, None, -100, 100] )
    fig.set_size_inches(48.5, 8.5)
    plt.show()
  
        
        
        

    
    
    
    
    
    
    
    
    
    
    
    
    

        
    
    
    
nodes = get_spline_nodes(100, 300, 30000, direction=0, yLimit=(10,1500), chaosFactor=0.8, xStretch=0)
        


x = nodes[:,0]
y = nodes[:,1]

tck,u     = interpolate.splprep( [x,y] ,s = 0 )
xSpline, ySpline = interpolate.splev( np.linspace( 0, 1, 100000 ), tck,der = 0)

total = 0
for i in range(len(xSpline)-1):
    total = total + np.sqrt( (xSpline[i+1] - xSpline[i])**2 + (ySpline[i+1] - ySpline[i])**2 )


fig, ax = plt.subplots()
fig.set_size_inches(8.5, 8.5)

plt.plot( x, y, 'o', xSpline, ySpline, linewidth=0.4)
plt.legend( [ 'data' , 'spline'] )
plt.axis( [ -500, 1500, -2000, 2000] )
plt.show()
print(total)



