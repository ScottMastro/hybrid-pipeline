import readCSV as reader
import plotBlocks as plot
import parser
from block import Chunk
from block import Block

import pandas as pd
import sys
 
#----globals-------------------------
chunkSize = 1000
chunkPerBlock = 3        #min chunks per block, theta1
maxDist = 5000           # #max dist between chunks, theta2
chunkSkip = 2            #max chunk skip allowed, theta3
bpSkip = 2000            #chunk skip size, theta4
repeatOverlap = 0.25     #max overlap for a repeated chunk, theta5

#----globals-------------------------

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
else:
    summary_file = "/media/scott/Rotom/assembly_data/CF062/summary.txt"
    block_file = "/media/scott/Rotom/assembly_data/CF062/blocks_new.txt"

index = 0
aligndf = parser.parse_alignments(summary_file)   
#data = data.head(48)

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


def main():
    
    prev_contig=""
    contig=""
    blockIndex = 0
    
    contig_blocks = []
    all_blocks = []
    to_plot=False
    
  while(index < len(aligndf)):
    prev_contig=contig
    
    if(index % 500 == 0):
        print(str(index) + " / " + str(len(aligndf)))
    
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
        if(index < len(aligndf)):
            continue
    
    contig_blocks.sort()
    contig_blocks = remove_redundancy(contig_blocks)
    all_blocks.append(contig_blocks)
    contig_blocks = blocks
    
    df = pd.concat([blocks2df(x) for x in all_blocks])
    
    df.to_csv(block_file, index=False, sep=",")
  
  
if __name__== "__main__":
  main()
  print("done")
