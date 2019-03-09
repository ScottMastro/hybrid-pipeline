import gfapy

t="\t"

class GFA:

    gfa=0

    def __init__(self):
        self.gfa = gfapy.Gfa(version='gfa1')

    def add_segments(self, ids, lengths):
        '''
        Adds a set of segments to gfa. 
        ids and lengths are lists that are zipped to create a dict (unique ids)
        '''
        d = dict(zip(ids, lengths))
        for id, length in d.items():  
            self.gfa.add_line("S"+t+str(id)+t+"*"+t+"LN:i:"+str(length))
            
    def add_segment(self, id, length, depth=0):
                
        dpString=""
        
        if(depth > 0):
            dpString=t+"RC:i:"+str(depth*length)
        
        self.gfa.add_line("S"+t+str(id)+t+"*"+t+"LN:i:"+str(length) +dpString)
            
    def add_containment(self, containerId, containedId, startPos, length, containerDir="+", containedDir="+"):
        self.gfa.add_line("C" +t+ str(containerId) +t+ containerDir +t+ \
                          str(containedId) +t+ containedDir +t+ \
                          str(startPos) +t+ str(length)+"M")
        
    def add_link(self, fromId, toId, overlap="*", fromDir="+", toDir="+"):
        self.gfa.add_line("L" +t+ str(fromId) +t+ fromDir +t+ \
                          str(toId) +t+ toDir +t+ overlap)
            
            
    def write(self, filename):
        self.gfa.to_file(filename)

    def __repr__(self):
        return self.gfa.__repr__()
    def __str__(self):
        return self.gfa.__str__()

    
def add_contig(contig, gfa, isQuery=True):
    
    if(contig.is_empty()):
        gfa.add_segment(contig.id, contig.size, depth=1)
        return gfa

    index = 1
    
    for i in range(len(contig.mblocks)):

        if isQuery:
            
            if i == 0:
                size=contig.mblocks[i].left_q() 
                gfa.add_segment(contig.id + "_" + str(index), size, depth=1)
                index = index + 1
                
            size=abs(contig.mblocks[i].right_q() - contig.mblocks[i].left_q())
            
            gfa.add_segment(contig.id + "_" + str(index), size, depth=1)
            gfa.add_segment(contig.mblocks[i].rid + "_" + contig.id, size, depth=10)
            gfa.add_link(contig.id + "_" + str(index-1), contig.id + "_" + str(index))
            gfa.add_link(contig.id + "_" + str(index-1), contig.mblocks[i].rid + "_" + contig.id)

            index = index + 1            
            
            size= (contig.size if i == (len(contig.mblocks)-1) else contig.mblocks[i+1].left_q()) - contig.mblocks[i].right_q()
            gfa.add_segment(contig.id + "_" + str(index), size, depth=1)

            gfa.add_link(contig.id + "_" + str(index-1), contig.id + "_" + str(index))
            gfa.add_link(contig.mblocks[i].rid + "_" + contig.id, contig.id + "_" + str(index))
            index = index + 1            

        else:
            #todo
            pass
            
    return gfa

