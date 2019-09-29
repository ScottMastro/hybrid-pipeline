import gfapy
import re
import uuid

t="\t"

class GFA2:
    def __init__(self):
        self.gfa = gfapy.Gfa(version='gfa2')
        self.paths=[]
        self.segments = dict()
            
    def add_segment(self, segID, length, depth=1):
        
        if segID in self.segments:
            return
        
        self.segments[segID] = True
        
        dpString=""
        if(depth > 0):
            dpString=t+"RC:i:"+str(depth*length)
        
        self.gfa.add_line("S" +t+ str(segID) +t+ str(length) +t+ "*" + dpString)
                    
    def add_edge(self, edge1, edge2, beg1, end1, beg2, end2, dir1="+", dir2="+"):
        self.gfa.add_line("E" +t+ "*" +t+ str(edge1)+dir1 +t+ \
                          str(edge2)+dir2 +t+ \
                          str(beg1) +t+ str(end1) +t+ \
                          str(beg2) +t+ str(end2) +t+ "*")
            
    def write(self, filename):
        self.gfa.to_file(filename)
        #f=open(filename,"a")
        #f.writelines("\n".join(self.paths))
        #f.close()
        
    def __repr__(self):
        return self.gfa.__repr__()
    def __str__(self):
        return self.gfa.__str__()

class GFA:
    def __init__(self):
        self.gfa = gfapy.Gfa(version='gfa1')
        self.paths=[]
        self.segments = dict()
        self.links = dict()

    def add_segments(self, ids, lengths):
        '''
        Adds a set of segments to gfa. 
        ids and lengths are lists that are zipped to create a dict (unique ids)
        '''
        d = dict(zip(ids, lengths))
        for id, length in d.items():  
            self.gfa.add_line("S"+t+str(id)+t+"*"+t+"LN:i:"+str(length))
            
    def add_segment(self, id, length, depth=0):
        if id not in self.segments:
            self.segments[id] = True
            dpString=""
            
            if(depth > 0):
                dpString=t+"RC:i:"+str(depth*length)
            
            self.gfa.add_line("S"+t+str(id)+t+"*"+t+"LN:i:"+str(length) +dpString)
        
    def add_containment(self, container, contained, dir1, dir2, pos1, pos2):
        self.gfa.add_line("C" +t+ str(container) +t+ dir1 +t+ \
                          str(contained) +t+ dir2 +t+ \
                          str(pos1) +t+ str(pos2-pos1) + "M")

    def add_path(self, name, segIDs):
    
        self.paths.append("P" +t+ str(name) +t+ "+,".join(segIDs) +"+" +t+ \
                          ",".join(["*" for x in segIDs]))
            

    def add_link(self, fromId, toId, overlap="*", fromDir="+", toDir="+"):
        if fromId+toId not in self.links:
            self.links[fromId+toId] = True

            self.gfa.add_line("L" +t+ str(fromId) +t+ fromDir +t+ \
                              str(toId) +t+ toDir +t+ overlap)
            
            
    def write(self, filename):
        self.gfa.to_file(filename)
        #f=open(filename,"a")
        #f.writelines("\n".join(self.paths))
        #f.close()
        
    def __repr__(self):
        return self.gfa.__repr__()
    def __str__(self):
        return self.gfa.__str__()

    
def get_seg(nodeA, nodeB, lengthData):
    
    id = str(nodeA.fork.after_id())
    start = nodeA.fork.after_pos()
    if nodeA.fork.after_strand() == -1:
        start = lengthData[id] - start              
    
    if id != str(nodeB.fork.before_id()):
        return (None, None)
    
    end = nodeB.fork.before_pos()
    if nodeB.fork.before_strand() == -1:
        end = lengthData[id] - end     
        
    segID = re.sub("\D", "", id) + "[" + str(start).zfill(10)  + "-" + str(end).zfill(10) + "]" + \
        re.sub("\D", "", str(nodeA.fork.before_id()))
        
    length = abs(end - start)
    return (segID, length)

def get_seg_head(node, lengthData, switch=False):
    
    if switch:
        id = str(node.fork.after_id())
        end = node.fork.after_pos()
    if not switch:
        id = str(node.fork.before_id())
        end = node.fork.before_pos()
        
    start = 0
        
    segID = re.sub("\D", "", id) + "[" + str(start).zfill(10)  + "-" + str(end).zfill(10) + "]"
    length = abs(end - start)
    return (segID, length)

def get_seg_tail(node, lengthData, switch=False):
    
    if not switch:
        id = str(node.fork.after_id())
        start = node.fork.after_pos()
    if switch:
        id = str(node.fork.before_id())
        start = node.fork.before_pos()
        
    end = lengthData[id]
        
    segID = re.sub("\D", "", id) + "[" + str(start).zfill(10)  + "-" + str(end).zfill(10) + "]"
    length = abs(end - start)
    return (segID, length)


def add_segs(gfa, currentNode, lengthData, prevSegIDs = [], skipNextMain=False):
    
    skip = True
    for prevNode in currentNode.prevNode:
        if skip:
            skip = False
            continue
        prevIDs = add_segs(gfa, prevNode, lengthData, prevSegIDs=[], skipNextMain=True)
        prevSegIDs = prevSegIDs + prevIDs
        
    currentSegIDs = []
    #tailSegIDs = []
    #headSegIDs = []

    for nextNode in currentNode.nextNode:
        
        segID, length = get_seg(currentNode, nextNode, lengthData)
        if segID is not None:
            gfa.add_segment(segID, length, depth=10)
            currentSegIDs.append(segID)

        segID, length = get_seg(nextNode, currentNode, lengthData)
        if segID is not None:
            gfa.add_segment(segID, length, depth=1)
            currentSegIDs.append(segID)
       # else:
            #segID, length = get_seg_tail(currentNode, lengthData, switch=True)
       #     gfa.add_segment(segID, length, depth=1)
        #    tailSegIDs.append(segID)
            
       #     segID, length = get_seg_head(nextNode, lengthData, switch=True)
      #      gfa.add_segment(segID, length, depth=1)
       #     headSegIDs.append(segID)

   # if len(currentNode.nextNode) == 0:
        
 #       segID, length = get_seg_tail(currentNode, lengthData)
 #       gfa.add_segment(segID, length, depth=1)
 #       currentSegIDs.append(segID)

 #       segID, length = get_seg_tail(currentNode, lengthData, switch=True)
 #       gfa.add_segment(segID, length, depth=10)
 #       currentSegIDs.append(segID)

 #   if len(prevSegIDs) == 0:
        
   #     segID, length = get_seg_head(currentNode, lengthData)
 #       gfa.add_segment(segID, length, depth=10)
  #     prevSegIDs.append(segID)
#
  #      segID, length = get_seg_head(currentNode, lengthData, switch=True)
   #     gfa.add_segment(segID, length, depth=1)
 #       prevSegIDs.append(segID)

    for segID in prevSegIDs:
        for nextSegID in currentSegIDs: # + tailSegIDs:
            gfa.add_link(segID, nextSegID)
                    
    skip = skipNextMain
    for nextNode in currentNode.nextNode:
        if skip:
            skip = False
            continue
        add_segs(gfa, nextNode, lengthData, currentSegIDs ) #+ headSegIDs)
    
    return currentSegIDs 

def add_segs2(gfa, currentNode, lengthData, visited=dict(), combine=[], link=None):
    
    if currentNode in visited:
        return []
    visited[currentNode] = True
    
    if len(currentNode.nextNode) == 1:
        nextNode = currentNode.nextNode[0] 
        if nextNode not in visited:
            ids = set()
            ids.add(currentNode.fork.before_id())
            ids.add(currentNode.fork.after_id())
            ids.add(nextNode.fork.before_id())
            ids.add(nextNode.fork.after_id())
            if len(ids) == 2:
                return add_segs2(gfa, nextNode, lengthData, visited, combine + [currentNode], link)
                
    id1 = currentNode.fork.before_id()
    id2 = currentNode.fork.after_id()
    start1 = currentNode.fork.get_pos_by_id(id1)
    if currentNode.fork.get_strand_by_id(id1) == -1: start1 = lengthData[id1] - start1
    start2 = currentNode.fork.get_pos_by_id(id2)
    if currentNode.fork.get_strand_by_id(id2) == -1: start2 = lengthData[id2] - start2
    end1 = start1
    end2 = start2

    for node in combine:
        p1 = node.fork.get_pos_by_id(id1)
        if node.fork.get_strand_by_id(id1) == -1: p1 = lengthData[str(id1)] - p1
        p2 = node.fork.get_pos_by_id(id2)
        if node.fork.get_strand_by_id(id2) == -1: p2 = lengthData[str(id2)] - p2

        if p1 < start1: start1 = p1
        if p2 < start2: start2 = p2
        if p1 > end1: end1 = p1
        if p2 > end2: end2 = p2

    segId1 = clean_id(id1) + "[" + str(start1) + "," + str(end1) + "]"
    segId2 = clean_id(id2) + "[" + str(start2) + "," + str(end2) + "]"
    #segIds.append(segId1)
    #segIds.append(segId2)
   
    gfa.add_segment(segId1, end1-start1, depth=1000)
    gfa.add_segment(segId2, end1-start1, depth=1000)

    if link is None:
        s1 = "_"+str(uuid.uuid1())
        gfa.add_segment(s1, 1, depth=1)
        link = s1
        
    gfa.add_link(link, segId1)
    gfa.add_link(link, segId2)

    s2 = "_"+str(uuid.uuid1())
    gfa.add_segment(s2, 1, depth=1)
    gfa.add_link(segId1, s2)
    gfa.add_link(segId2, s2)
    
    for nextNode in currentNode.nextNode:
        add_segs2(gfa, nextNode, lengthData, visited, [], s2)    



def get_containment(gfa, node1, node2, lengthData, queryData):
    
    if node1.fork.before_id() != node2.fork.after_id():
        return
    
    id1 = str(node1.fork.before_id())
    id2 = str(node2.fork.before_id())

    beg1 = node1.fork.before_pos()
    if node1.fork.before_strand() == -1:
        beg1 = lengthData[id1] - beg1
    end1 = node2.fork.after_pos()
    if node2.fork.after_strand()== -1:
        end1 = lengthData[id1] - end1  
        
    beg2 = node1.fork.after_pos()
    if node1.fork.after_strand()== -1:
        beg2 = lengthData[id2] - beg2
    end2 = node2.fork.before_pos()
    if node2.fork.before_strand()== -1:
        end2 = lengthData[id2] - end2  
        
    if end1 < beg1 : (beg1, end1) = (end1, beg1)
    if end2 < beg2 : (beg2, end2) = (end2, beg2)

    if id1 in queryData:
        container = id1
        contained = id2
        pos1 = beg1
        pos2 = end1
        x = beg2
        y = end2
        
    if id2 in queryData:
        container = id2
        contained = id1
        pos1 = beg2
        pos2 = end2
        x = beg1
        y = end1
        
    gfa.add_containment(clean_id(container), clean_id(contained) + str(x) + str(y),
                 "+", "+", pos1, pos2)

def add_segs3(gfa, currentNode, lengthData, visited, queryData):
    
    if currentNode in visited:
        return
    
    visited[currentNode] = True
    
    id1 = str(currentNode.fork.after_id())
    if id1 in queryData:
        gfa.add_segment(clean_id(id1), lengthData[id1])
    id2 = str(currentNode.fork.before_id())
    if id2 in queryData:
        gfa.add_segment(clean_id(id2), lengthData[id2])

    for prevNode in currentNode.prevNode:
        if prevNode in visited:
            continue
        
        get_containment(gfa, prevNode, currentNode, lengthData, queryData)


    for nextNode in currentNode.nextNode:
        if nextNode in visited:
            continue

        get_containment(gfa, currentNode, nextNode, lengthData, queryData)

    
    for node in currentNode.prevNode + currentNode.nextNode:
        if node not in visited:
            add_segs3(gfa, node, lengthData, visited, queryData)


def clean_id(id): return re.sub("\D", "", str(id))

def graph_to_gfa(graph, lengthData, queryData):
    
    
    gfa = GFA()
    for key in refData:
        gfa.add_segment(key, lengthData[key], depth=1)


    #add_segs(gfa, graph.get_start(), lengthData)
    add_segs2(gfa, graph.start, lengthData)

    return gfa


def get_edge(gfa, node1, node2, lengthData):
    
    if node1.fork.before_id() != node2.fork.after_id():
        return
    
    id1 = str(node1.fork.before_id())
    id2 = str(node2.fork.before_id())

    beg1 = node1.fork.before_pos()
    if node1.fork.before_strand() == -1:
        beg1 = lengthData[id1] - beg1
    end1 = node2.fork.after_pos()
    if node2.fork.after_strand()== -1:
        end1 = lengthData[id1] - end1  
        
    beg2 = node1.fork.after_pos()
    if node1.fork.after_strand()== -1:
        beg2 = lengthData[id2] - beg2
    end2 = node2.fork.before_pos()
    if node2.fork.before_strand()== -1:
        end2 = lengthData[id2] - end2  
        
    if end1 < beg1 : (beg1, end1) = (end1, beg1)
    if end2 < beg2 : (beg2, end2) = (end2, beg2)

    gfa.add_edge(clean_id(id1), clean_id(id2), beg1, end1, beg2, end2, dir1="+", dir2="+")


def add_edges(gfa, currentNode, lengthData, visited):
    
    if currentNode in visited:
        return
    visited[currentNode] = True

    id1 = str(currentNode.fork.after_id())
    gfa.add_segment(clean_id(id1), lengthData[id1])
    id2 = str(currentNode.fork.before_id())
    gfa.add_segment(clean_id(id2), lengthData[id2])

    for prevNode in currentNode.prevNode:
        if prevNode in visited:
            continue
        get_edge(gfa, prevNode, currentNode, lengthData)
        
    for nextNode in currentNode.nextNode:
        if nextNode in visited:
            continue
        get_edge(gfa, currentNode, nextNode, lengthData)

    for node in currentNode.prevNode + currentNode.nextNode:
        add_edges(gfa, node, lengthData, visited)


def graph_to_gfa2(graph, lengthData):
    
    gfa = GFA2()

    #add_segs(gfa, graph.get_start(), lengthData)
    visited = dict()    
    add_edges(gfa, graph.start, lengthData, visited)
        
    return gfa