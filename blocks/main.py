import plotBlocks as plot
import file_handler as fh
#from chunk import Chunk
#from block import Block
#from block_builder import construct_blocks
#from block_builder import remove_redundancy
#from block_builder import evaluate_blocks
import pandas as pd
import sys
 
#----globals-------------------------
chunkSize = 1000
chunkPerBlock = 3        #min chunks per block, theta1
maxDist = 5000           #max dist between chunks, theta2
chunkSkip = 2            #max chunk skip allowed, theta3
bpSkip = 2000            #chunk skip size, theta4
repeatOverlap = 0.25     #max overlap for a repeated chunk, theta5


if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
else:
    #prefix="C:/Users/scott/Desktop/novafix/blocks"
    prefix="/home/scott/Dropbox/hybrid-pipeline/blocks"
    summary_file = prefix + "/summary.txt"
    block_file = prefix + "/blocks_new.txt"

index = 0
aligndf = fh.parse_alignments(summary_file)   
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
    
    nlines=len(aligndf)
    prev_contig=""
    contig=""
    index = 0
    
    contig_blocks = []
    all_blocks = []
    to_plot=False
    
    while(index < nlines):
        prev_contig=contig
        
        if(index % 500 == 0):
            print(str(index) + " / " + str(nlines))
        
        line = aligndf[index:index+1]
        contig=line.iloc[0,0]

        #print(index)
        blocks, trash = construct_blocks(line)
        #evaluate_blocks(line, blocks, trash, print_=False)
        #time.sleep(0.5)

        index=index+1    
    
        if(to_plot):
            df = blocks2df(blocks)
            plot.plotBlocks(fh.line_to_df(line), df)
    
        if(prev_contig == contig):
            contig_blocks.extend(blocks)
            if(index < nlines):
                continue
    
        elif(prev_contig == ""):
            contig_blocks = blocks
            if(index < nlines):
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
