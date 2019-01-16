import readCSV as reader
import plotBlocks as plot
import pandas as pd
import csv
import sys

csv.field_size_limit(999999999999)
 
class Chunk:
    id=-1
    start=-1
    end=-1
    qstart=-1
    qend=-1

    def __init__(self, cid, cst, ced, qst, qed):
        self.id=cid
        self.start=cst
        self.end=ced
        self.qstart=qst
        self.qend=qed
    
    def get_dir(self):
        if(self.start < self.end):
            return 1
        elif(self.end < self.start):
            return -1
        else:
            return "???"
    
    def left(self): return min(self.start, self.end)
    def right(self): return max(self.start, self.end)
    
    def print_info(self):
        print(str(self.id) + " : " + str(self.start) + " to " + str(self.end))

class Block:
    index=-1
    direction=0
    chunks=[]
    #chunks in the contig
    nchunk=0
    #size of the chr
    span=0
    chrom=""
    contig=""
  
    def __init__(self, chunk, chrom="xxx", contig="tig", nchunk=0, span=0):
        global block_index
        self.index=block_index
        block_index=block_index+1
        self.direction=0
        self.chunks=[]
        self.chrom=chrom
        self.contig=contig
        self.nchunk=nchunk
        self.span=span

        self.why=""
        self.chunks.append(chunk)

    def get_dir(self): return self.direction   
    def verify_direction(self, chunk): 
        
        d=1
        if(chunk.left() > self.chunks[-1].left()):
            d=-1
        
        if self.direction == 0:
            self.direction = d  
            return True
        
        return d == self.direction
    
    def add(self, chunk): self.chunks.append(chunk)
    def first_chunkid(self): return self.chunks[-1].id
    def last_chunkid(self): return self.chunks[0].id
    def first_chunk(self): return self.chunks[-1]
    def last_chunk(self): return self.chunks[0]
    def trim_end(self, n):
        if(n >= (len(self.chunks) - theta1)):
            return True
        self.chunks = self.chunks[n:]
        return False
    
    def trim_start(self, n): 
        if(n >= (len(self.chunks)-theta1)):
            return True
        self.chunks = self.chunks[:-n]
        return False
    
    def left(self): 
        return min(self.chunks[0].left(), self.chunks[-1].right())
    def right(self): 
        return max(self.chunks[0].right(), self.chunks[-1].left())
    
    def size(self): return len(self.chunks)
    def valid(self): return len(self.chunks) > 0

    def print_info(self):
        print(str(self.first_chunkid()) + " to " + str(self.last_chunkid()) + " : " + str(self.size()))

    def __lt__(self, other):
        return self.chunks[0].id > other.chunks[0].id

def distance(chunk1, chunk2):
    return min(abs(chunk1.left() - chunk2.right()) )

def construct_blocks(line):
    
    blocks = []
    trash = []
    
    cid=[int(x) for x in line.iloc[0,4].split(',')]
    st=[int(x) for x in line.iloc[0,5].split(',')]
    ed=[int(x) for x in line.iloc[0,6].split(',')]
    qst=[int(x) for x in line.iloc[0,7].split(',')]
    qed=[int(x) for x in line.iloc[0,8].split(',')]
    
    chunks=[Chunk(cid[i], st[i], ed[i], qst[i], qed[i]) for i in range(0, len(cid))]
    
    contig_ = line.iloc[0,0]
    nchunk_ = line.iloc[0,1]
    chrom_ = line.iloc[0,2]
    span_ = line.iloc[0,3]

    while len(chunks) > 0:
         
        x = chunks.pop()
        
        block = Block(x, chrom_, contig_, nchunk_, span_)
        i = 1
                
        while len(chunks) > 0:

            if(i > len(chunks)):
                if(block.size() > theta1):
                    blocks.append(block)
                    break
                else:
                    trash.append(block)    
                    break
                                    
            next = chunks[-i]
                 
            id_dist = x.id - next.id
            
           # print(chunks[-i].id)
           # if(chunks[-i].id == 67235):
          #      time.sleep(20)
                
            if id_dist == 1 :
                if(abs(x.start - next.start) - chunk_size < theta2):
                    if(block.verify_direction(next)):
                        block.add(next)
                        x = chunks.pop(-i)
                        continue
                    
                i = i + 1
                continue
                
            elif id_dist == 0 :
                i = i + 1
                continue
                 
            else:
                chunk_skip = id_dist-1

                if(chunk_skip <= theta3):
                    #print(str(abs(x.start - next.start)- chunk_size) + " - " + str(chunk_skip*theta4))
                    if(abs(x.start - next.start) - chunk_size < (chunk_skip*theta4)):
                        if(block.verify_direction(next)):
                            block.add(next)
                            x = chunks.pop(-i)
                            continue
                        
                    i = i + 1
                    continue 
                    
                 
            if(block.size() > theta1):
                blocks.append(block)
                #block.print_info()

            else:
                trash.append(block)
            break
        
    return [blocks, trash]

def chunk_split(c1, c2):
    size = max(c1.qend, c2.qend) - min(c1.qstart, c2.qstart) + 1
    overlap_start = max(c1.qstart, c2.qstart)
    overlap_end = min(c1.qend, c2.qend)
    return (overlap_end - overlap_start) / size < theta5

def remove_redundancy(blocks):
    i=0
    while ( i < len(blocks)-1 ):
        
        x = blocks[i].first_chunkid() - blocks[i+1].last_chunkid()
        
        if x <= 0:
            #print(str(i) + "  " + str(x))

            if x == 0:
                if(chunk_split(blocks[i].first_chunk(), blocks[i+1].last_chunk())):
                    i=i+1
                    continue
            
            if(blocks[i].size() > blocks[i+1].size()):
                empty = blocks[i+1].trim_end(-x + 1)
                if(empty):
                    blocks.pop(i+1)
                    i=i-1
                    
            else:
                empty = blocks[i].trim_start(-x + 1)
                if(empty):
                    blocks.pop(i)
                    i=i-1

        i=i+1
        
    return blocks
        
        
def blocks2df(blocks):
    
    df=pd.DataFrame()
    df["block"] = [x.index for x in blocks]
    df["contig"] = [x.contig for x in blocks]
    df["chrom"] = [x.chrom for x in blocks]
    df["span"] = [x.span for x in blocks]
    df["size"] = [x.size() for x in blocks]
    df["chunk_start"] = [x.first_chunkid() for x in blocks]
    df["chunk_end"] = [x.last_chunkid() for x in blocks]
    df["total_chunk"]= [x.nchunk for x in blocks]
    df["left"] = [x.left() for x in blocks]
    df["right"] = [x.right() for x in blocks]
    df["direction"] = [x.direction for x in blocks]
    #df["why"] = [x.why for x in blocks]

    return df

#----globals-------------------------
block_index = 0
chunk_size = 1000

#min chunks per block
theta1 = 3
#max dist between chunks
theta2 = 5000

#max chunk skip allowed
theta3 = 2
#chunk skip size
theta4 = 2000

#max overlap for a repeated chunk
theta5 = 0.25

#----globals-------------------------

index = 0
#in_file = "../input/blastn.canu_vs_GRCh38.dist.II.concat.id95.cov40.filter.sort.summarize.txt"
in_file=sys.argv[1]
out_file = sys.argv[2]
#out_file = "../out/canu_hg38_blocks.txt"

data = reader.read(in_file)   
#data = data.head(48)
#data = data.tail(24)

prev_contig=""
contig=""

contig_blocks = []
all_blocks = []
to_plot=False

while(index < len(data)):
    prev_contig=contig
    
    if(index % 500 == 0):
        print(str(index) + " / " + str(len(data)))
    
    line = data[index:index+1]
    contig=line.iloc[0,0]
    
    blocks, trash = construct_blocks(line)
    index=index+1    

    if(to_plot):
        df = blocks2df(blocks)
        plot.plotBlocks(reader.line2df(line), df)

    if(prev_contig == contig):
        contig_blocks.extend(blocks)
        if(index < len(data)):
            continue

    elif(prev_contig == ""):
        contig_blocks = blocks
        if(index < len(data)):
            continue
    
    contig_blocks.sort()
    contig_blocks = remove_redundancy(contig_blocks)
    all_blocks.append(contig_blocks)
    contig_blocks = blocks

df = pd.concat([blocks2df(x) for x in all_blocks])

df.to_csv(out_file, index=False, sep=",")

print("done")
