#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import libraries
#----------------
import sys

import pandas as pd
import gzip
from Bio import SeqIO
import canufasten
import novafix
import time

#input parameters
#----------------
blockThreshold=10         # % overlap to establish connection in matrix
chunkSize=1000           # chunk length in bp
overlapThreshold=6000    # bp distance to end of nova tig before scaffolding
alignBuffer = 8000       # number of bp aligned when scaffolding
baseBuffer = 200         # number of extra canu bases to take when scaffolding

novaBuffer = 1000         # number of bp aligned when N filling (nova)
canuBuffer = 1600         # number of bp aligned when N filling (canu)

#input files
#----------------

if len(sys.argv) > 1 :
    block_file=sys.argv[1]
    summary_file=sys.argv[2]
    canufa=sys.argv[3]
    novafa=sys.argv[4]
    outdir=sys.argv[5]
else:
    prefix = "/media/scott/Rotom/assembly_data/CF062/"
    #prefix = "../"
    
    block_file = prefix + "blocks.txt"
    block_file = prefix + "blocks.txt"
    summary_file = prefix + "summary.txt"
    novafa = prefix + "OSK7121_03C_supernova.pseudohap2.2.fasta"
    canufa = prefix + "CF062B2D.contigs.fasta.PILON2.fasta"
    outdir = prefix + "canufasten"

output_base = outdir + "/fastened"
fa = ".fasta"

#load data
#----------------

print("Loading data")
blockdf = pd.read_csv(block_file, dtype={'contig': str, 'chrom': str})
aligndf = pd.read_csv(summary_file, header = None, sep='\t')

print("Reading Supernova fasta")
if(novafa[-2:] == "gz"):
    with gzip.open(novafa, "rb") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
else:
    records = list(SeqIO.parse(novafa, "fasta"))
novaData = dict(zip([r.id for r in records], [r.seq for r in records]))

print("Reading Canu fasta")
if(canufa[-2:] == "gz"):
    with gzip.open(canufa, "rb") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
else:
    records = list(SeqIO.parse(canufa, "fasta"))
canuData = dict(zip([r.id for r in records], [r.seq for r in records]))
records = None

#get all contigs
#----------------
novaTigs = [x for x in novaData]
canuTigs = [x for x in canuData]
contigs = []
contigs.extend(novaTigs)
contigs.extend(canuTigs)

#scaffold
#---------------
paths = canufasten.build_network(novaTigs, canuTigs, novaData, canuData, blockdf, 
                    blockThreshold, chunkSize, overlapThreshold, alignBuffer, baseBuffer)


nfixes = novafix.fix_Ns(novaTigs, canuTigs, novaData, canuData, blockdf, aligndf,
                    novaBuffer, canuBuffer, blockThreshold, chunkSize)


    
   
#output results
#----------------

#segment =('0', 196, 5456171, 1)
def fill_Ns(segment):
    newSegment = []
    if segment[0] in nfixes:
        fixes = nfixes[segment[0]]
        for f in fixes:
            l, m, r = f
            #todo: what if segment ends in the middle of Ns????
            if l[2] > segment[1] and r[1] < segment[2] and segment[3] == l[3]:
                sl = (segment[0], segment[1], l[2], segment[3])
                sr = (segment[0], r[1], segment[2], segment[3])
                newSegment.append
                
def collapse_path(path):
    segments = []
    collapsed = []
    
    for bridge in path:
        for segment in bridge:
            segments.append(segment)
    
    prevSeg = segments[0]
    for i in range(1, len(segments)):
        seg = segments[i]
        
        if seg[0] == prevSeg[0]:       #same contig
            if seg[3] != prevSeg[3]:   #different strands
                print("WARNING: unexpected bridge between reverse complements")
            else:
                prevSeg = (prevSeg[0], prevSeg[1], seg[2], prevSeg[3])
                continue
        
        collapsed.append(prevSeg)
        prevSeg = seg
        
    collapsed.append(prevSeg)
    return collapsed
            
for path in paths:
    print(collapse_path(path))
            

            
def output_left_path(path, novaSeq, canuSeq):
    with open('guru99.txt', 'w') as f:
    
        if path is not None:
            c, n = path
            printSize=5000
            dist=10000
            if(c[2] > 0):
    
                f.writelines(">canu_-" + str(dist) + "_" + c[0] + "\n")
                f.writelines(canuSeq[c[2]-printSize-dist:c[2]-dist]+ "\n")
    
                
                f.writelines(">canu_" + c[0] + "\n")
                f.writelines(canuSeq[c[2]-printSize:c[2]]+ "\n")
                f.writelines(">nova_" + n[0]  + "\n")
                f.writelines(novaSeq[n[1]:n[1]+printSize]+ "\n" )
                
                f.writelines(">nova_+" + str(dist) + "_" + n[0] + "\n")
                f.writelines(novaSeq[n[1]+dist:n[1]+dist+printSize]+ "\n" )
                
                f.writelines("\t\n")
                f.writelines("\t\n")
                f.writelines("\t\n")
    



    if len(matches < 1):
        full = [(novaId, 0, len(novaSeq)-1)]
        printPath.append(full)
#        continue
       




