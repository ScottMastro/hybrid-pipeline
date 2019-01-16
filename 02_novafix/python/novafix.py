import pandas as pd
from contigMatrix import constructMatrix

import numpy as np
import matplotlib.pyplot as plt
import re
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import subprocess
import sys

seq_buffer = 500
chunk_size = 1000

print("Starting")
print("Loading data")

block_file=sys.argv[1]
summary_file=sys.argv[2]
canufa=sys.argv[3]
novafa=sys.argv[4]
outdir=sys.argv[5]

blockdf = pd.read_csv(block_file) 
aligndf = pd.read_csv(summary_file, header = None, sep='\t')

output = outdir+"/nfixed.fa"
canu_part_output = outdir+"/nfixed_canu.fa"
nova_part_output = outdir+"/nfixed_nova.fa"

#with gzip.open(novafa, "rt") as handle:
#    records = list(SeqIO.parse(handle, "fasta"))

records = list(SeqIO.parse(novafa, "fasta"))

nova_reads = dict(zip([r.id for r in records], [r.seq for r in records]))

print("Creating matrix")

matrix = constructMatrix(blockdf, threshold=35)

def find(s, ch):
    nloc = dict()

    index = [i for i, ltr in enumerate(s) if ltr == ch]
    if(len(index) < 1): return nloc
    last_start = index[0]
    prev = last_start -1
    
    for i in index:
        if(abs(prev - i) != 1):
            nloc[last_start] = prev
            last_start = i
        
        prev = i
        
    return nloc
    

def trim(aln, base_buffer):
    canu_aln = aln[0][0].rstrip()
    nova_aln = aln[0][1].rstrip()
    
    while(nova_aln[0] == "-" or canu_aln[0] == "-"):
        canu_aln=canu_aln[1:]
        nova_aln=nova_aln[1:]

    while(nova_aln[-1] == "-" or canu_aln[-1] == "-"):
        canu_aln=canu_aln[:-1]
        nova_aln=nova_aln[:-1]
        
    nidx = nova_aln.index("N")
    counter = base_buffer
    left = nova_aln[:nidx]
    left_buffer = ""
    index = nidx-1
    while(index >=0 and counter > 0):
        if not canu_aln[index] == "-":
            left_buffer = canu_aln[index] + left_buffer
            counter=counter-1
            
        left=left[:-1]
        index=index-1
    
    left = left + left_buffer
    
    counter = base_buffer
    right = nova_aln[nidx:]
    right_buffer = ""
    index = nidx-1
    while(index < len(canu_aln) and counter > 0):
        if not canu_aln[index] == "-":
            right_buffer = right_buffer +  canu_aln[index]
            counter=counter-1
            
        right=right[1:]
        index=index+1
        
    right = right_buffer + right
    trimmed = (left + right).replace("N", "").replace("-", "")
    return trimmed
        
    
def stitch(aln, pre, post):
    base_buffer=7
    
    pre = pre[:-base_buffer]
    post = post[base_buffer:]

    _pre = pre
    _post = post
    canu_aln = aln[0][0].rstrip()
    nova_aln = aln[0][1].rstrip()

    while(len(pre) > 0):
        if pre[0] == nova_aln[0]:
            pre = pre[1:]
    
        canu_aln = canu_aln[1:]
        nova_aln = nova_aln[1:]
        
    while(len(post) > 0):
        if post[-1] == nova_aln[-1]:
            post = post[:-1]
    
        canu_aln = canu_aln[:-1]
        nova_aln = nova_aln[:-1]

    canu_aln = canu_aln.replace("-", "")
    
    if(len(canu_aln) <= base_buffer*2 ):
        return trim(aln, base_buffer)
    return _pre + canu_aln + _post


def create_N_fasta():

    fasta_file = "../out/ngaps.fa"
    f = open(fasta_file, "w")

    tigcounter = 1
    for novaTig in nova_reads:
        print(novaTig + "(" + str(tigcounter) + ")")
        tigcounter=tigcounter+1
        nloc = find(nova_reads[novaTig], 'N')
        if(len(nloc) < 1): continue
        
        counter = 1
        for start in nloc.keys():
    
            seq_start = max(0, start - seq_buffer)
            seq_end = min(nloc[start] + seq_buffer, len(nova_reads[novaTig]))
            novaSeq = str(nova_reads[novaTig][seq_start:seq_end])
            novaSeq = re.sub(r'N+', 'N', novaSeq)
            fa_id = ">" + novaTig + "_gap" + str(counter)
            
            f.write(fa_id + "\n")
            f.write(novaSeq + "\n")
            
            
            counter=counter+1

#STEPS:
# 1. Make new empty fasta file
# 2. Iterate through nova_reads contig-by-contig --> stop when you get to an N
# 3. Divide index of the N by 1000 then + 1--> This is the chunk number
# 4. Identify blocks with the same chunk_start and end as chunk # that align to the same canu tig --> if the chunks found don't correspond, keep +1/-1 until they do
# 5. Find left-->right of missing region
# 6. in /allen/novafix do bash sh_extractcanu *canu tig name minus pilon_pilon* *left* *right*
# 7. Replace Ns with the canu sequence in the new fasta file --> back to step 2

#   novaTig = "290"

f = open(output, "w")
fcanu = open(canu_part_output, "w")
fnova = open(nova_part_output, "w")
#hybrid = open(outdir+"/hybrid.fa")

for novaTig in nova_reads:
    nloc = find(nova_reads[novaTig], 'N')
    if(len(nloc) < 1): continue
    counter=0
    for start in nloc.keys():
        counter=counter+1
        seq_start = max(0, start - seq_buffer)
        seq_end = min(nloc[start] + seq_buffer, len(nova_reads[novaTig]))

        novaSeq = str(nova_reads[novaTig][seq_start:seq_end])
        
        chunk_start = seq_start/chunk_size + 1
        chunk_end = seq_end/chunk_size + 1
        
        blocks_start = blockdf.loc[(blockdf['contig'] == int(novaTig))]
        blocks_start = blocks_start[blocks_start['chunk_start'] <= chunk_start]
        blocks_start = blocks_start[blocks_start['chunk_end'] >= chunk_start]
        
        blocks_end = blockdf.loc[(blockdf['contig'] == int(novaTig))]
        blocks_end = blocks_end[blocks_end['chunk_start'] <= chunk_end]
        blocks_end = blocks_end[blocks_end['chunk_end'] >= chunk_end]

        potential_canu = set(blocks_start['chrom']) & set(blocks_end['chrom'])

        try:
            in_matrix = set(matrix[int(novaTig)].dropna().keys())
        except KeyError:
            "ignored"
            #fa_id = ">" + novaTig + "_unchanged" 
            #hybrid.write(fa_id + "\n")
            #hybrid.write(str(nova_reads[novaTig]) + "\n")
            continue
            
        canuTig = list(potential_canu & in_matrix)

        if(len(canuTig) == 1):
            try:
                canuTig = canuTig[0]
                alignments = aligndf[aligndf[2] == canuTig]
                alignments = alignments[alignments[0] == int(novaTig)]  
                
                chunks = alignments.iloc[0,4].split(",")
                chunk_start_index = chunks.index(str(chunk_start))
                chunk_end_index = chunks.index(str(chunk_end))
                
                canu_starts = alignments.iloc[0,5].split(",")
                canu_ends = alignments.iloc[0,6].split(",")
               
                canuStart = min( int(canu_starts[chunk_start_index]), int(canu_starts[chunk_end_index]), int(canu_ends[chunk_start_index]), int(canu_ends[chunk_end_index]) )
                canuEnd = max( int(canu_starts[chunk_start_index]), int(canu_starts[chunk_end_index]), int(canu_ends[chunk_start_index]), int(canu_ends[chunk_end_index]) )

                canuSeq = subprocess.check_output(['bash','extractCanu.sh', canufa, canuTig.split("_")[0], str(canuStart - seq_buffer), str(canuEnd + seq_buffer)]).rstrip()

                if(len(canuSeq) > 10000):
                    print("too long")
                    continue
                if(len(canuSeq) < 1):
                    print("too short???")
                    continue
                
                pre = novaSeq[0:novaSeq.index("N")]
                post = novaSeq[novaSeq.rindex("N")+1:len(novaSeq)]

                aln = pairwise2.align.globalms(canuSeq, novaSeq, 2, -3, -5, -2, one_alignment_only=True)
                canuRev=str(Seq(canuSeq).reverse_complement())
                alnRev = pairwise2.align.globalms(canuRev, novaSeq, 2, -3, -5, -2, one_alignment_only=True) 

                fa_id = ">" + novaTig + "_gap" + str(counter)

                fcanu.write(fa_id + "\n")
                if(aln[0][2] > alnRev[0][2]):
                    stitched = stitch(aln, pre, post)
                    fcanu.write(canuSeq + "\n")
                else:
                    stitched = stitch(alnRev, pre, post)
                    fcanu.write(canuRev + "\n")

                fnova.write(fa_id + "\n")
                fnova.write(novaSeq + "\n")

                f.write(fa_id + "\n")
                f.write(stitched + "\n")
                
                print("Ns were fixed :)")
            except ValueError:
                "ignored"
                continue
