import numpy as np
import alignment_chunk as chunker
import alignment_block as blocker
import alignment_megablock as mblocker
import alignment_contig as contiger

from gfa_handler import GFA

def distance(chunk1, chunk2):
    return min(abs(chunk1.left() - chunk2.right()) )
    



def count_chunk_alignments(qdf, qname):
    nchunk = qdf.iloc[0,1]
    countDict = {n+1 : 0 for n in range(nchunk)}
    
    for row in qdf[str(qdf.columns[4])]:
        chunkIds=[int(i) for i in row.split(',')]
        for n in chunkIds:
            countDict[n] = countDict[n] + 1     
            
    import matplotlib.pylab as plt
    lists = sorted(countDict.items()) # sorted by key, return a list of tuples  
    x, y = zip(*lists) # unpack a list of pairs into two tuples
    plt.figure(figsize=(30,8))
    plt.plot(x, y)
    plt.show()
    return countDict
    
tab="\t"


def add_containments(contig, gfa, isQuery=True):
    
    for mblock in contig.mblocks:
        if isQuery:
            gfa.add_containment(contig.id, mblock.rid, mblock.left_q(), mblock.size_q())
        else:
            gfa.add_containment(contig.id, mblock.qid, min(mblock.left_r(), mblock.right_r()), mblock.size_r())
            
    return gfa
    
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


def create_blocks(aligndf, chunkSize=1000): #, qname):

    gfa = GFA()
    
    #add reference segments
    #gfa.add_segments(aligndf[str(aligndf.columns[2])], 
    #                 aligndf[aligndf.columns[3]])
    #add query segments
    #gfa.add_segments(aligndf[str(aligndf.columns[0])], 
    #                 aligndf[aligndf.columns[1]]*1000)
        
    qnames=np.unique(aligndf[str(aligndf.columns[0])])
    
    megablocks = { str(name) : [] for name in qnames}
    contigs = []
    
    #qname=0
    for qname in qnames:
        #dataframe containing only alignments for one query contig
        qdf = aligndf.loc[aligndf[str(aligndf.columns[0])] == qname]
        
        #row=next(qdf.iterrows())[1]
        for idx, row in qdf.iterrows():
            
            qid = str(row[0]) 
            qlen = int(row[1])*chunkSize
            rid = str(row[2])

            chunks = chunker.construct_chunks(row, chunkSize)
            blocks, trash = blocker.construct_blocks(rid, qid, chunks)
            blocks = blocker.remove_redundancy(blocks, chunkPerBlock, repeatOverlap)            
            mblocks = mblocker.construct_megablocks(row, blocks, trash, chunkSize, 1.5)
            if len(mblocks) > 0:
                megablocks[qid].extend(mblocks)
        
        print(qid)
        contigs.append(contiger.construct_contig(qid, qlen, megablocks[qid], isQuery=True))
        gfa = add_contig(contigs[-1], gfa, isQuery=True)
    
    gfa.write("hello.gfa")
    

import sys
import file_handler as fh
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