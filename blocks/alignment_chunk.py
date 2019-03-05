class Chunk:
    id=-1
    
    rstart=-1
    rend=-1
    qstart=-1
    qend=-1

    def __init__(self, cid, rst, red, qst, qed):
        self.id=cid
        self.rstart=rst
        self.rend=red
        self.qstart=qst
        self.qend=qed
    
    def get_dir(self, string=False):
        if(self.rstart < self.rend):
            return "+" if string else 1
        elif(self.rend < self.rstart):
            return "-" if string else -1
        else:
            return "*" if string else 0
    
    #left and right relative to query
    def left_q(self): return self.qstart
    def right_q(self): return self.qend
    def left_r(self): return self.rstart
    def right_r(self): return self.rend
    
    def size_q(self): return abs(self.qstart - self.qend)
    def size_r(self): return abs(self.rstart - self.rend)

    def __repr__(self):
        return str(self.id) + ":" + str(self.qstart) + "-" + str(self.qend)

    def __str__(self):
        return "chunk id=" + str(self.id) + \
                " query=" + str(self.qstart) + "-" + str(self.qend) + \
                " ref=" + str(self.rstart) + "-" + str(self.rend) + \
                " dir=" + self.get_dir(string=True)
    

def construct_chunks(row, chunkSize=1000):
    cid=[int(x) for x in row[4].split(',')]
    rst=[int(x) for x in row[5].split(',')]
    red=[int(x) for x in row[6].split(',')]
    qst=[int(x) for x in row[7].split(',')]
    qed=[int(x) for x in row[8].split(',')]
    
    #convert start/end position from chunk to contig
    qst=[(i-1)*chunkSize + s for i, s in zip(cid, qst)]
    qed=[(i)*chunkSize - (chunkSize-e) for i, e in zip(cid, qed)]

    return [Chunk(cid[i], rst[i], red[i], qst[i], qed[i]) for i in range(0, len(cid))]