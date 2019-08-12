import sys
sys.path.append("./analysis")
sys.path.append('./weld')
sys.path.append('./stitch')

from parameters import Parameters
import file_handler as io
import welder
import stitcher

import contig_scaffolder as scaffolder
import contig_output as output
import path_helper

import sys
import copy 
import re

import numpy as np
from Bio import SeqIO
import gzip
import pickle
import pandas as pd
from analysis import region_extract
import os
#----Parameters-------------------------

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
    novafa = sys.argv[3]      #supernova fasta
    canufa = sys.argv[4]      #canu fasta

#----Parameters-------------------------

def validate_ids(rids, qids):
    intersect = set(rids).intersection(qids) 
    if len(intersect) > 0:
        print("ERROR: IDs for input fasta files overlap.")
        print("\"" + str(intersect.pop()) + "\" is found in both fasta files.")
        print("Please ensure each ID is unique.")
        return False
    return True

def main():

    #----Load data-------------------------

    param = Parameters()

    print("Reading Canu fasta")
    refData = io.read_fasta(param.referenceFasta)
    print("Reading Supernova fasta")
    queryData = io.read_fasta(param.queryFasta)
    
    rids = list(refData.keys())
    qids = list(queryData.keys())

    #if not validate_ids(rids, qids): return
    
    seqData = dict()
    seqData.update(refData)
    seqData.update(queryData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    
    aligndf = io.parse_alignments(param.summaryFile)
    
    lrdict = scaffolder.linkedReadsDict(rids, param.arks_output)
    #lrdf = scaffolder.linkedReadsDF(rids, lengthData, param.arks_output)
    unitigs = scaffolder.unitigsDict(param.unitigsBed)

    #--------------------------------------
    
    paths = []
    emptyIds = []
    print("Iterating over query contigs...")

    outputPath="/Users/allen bao/Documents/assembly_data/plots/"
    outputPath=None

    for tigId in qids:
        
       # if(int(tigId) < 295): 
       #    continue
        print("Contig: " + tigId)
        contig = stitcher.stitch(tigId, aligndf, lengthData, param, \
                                 plotDir=outputPath, q=True)
        if contig is None: emptyIds.append(tigId)
        else:
            path = welder.weld(contig, seqData, lengthData, param)
            

            if len(path) < 1: emptyIds.append(tigId)
            else: 
                paths.append(path)
    
    
    with open('CF062_paths2.pickle', 'wb') as handle:
        pickle.dump(paths, handle)
        
    with open('CF062_paths.pickle', 'rb') as handle:
        paths = pickle.load(handle)
    
     
    '''
    counter=1
    counterStart=counter
    size=20000
    i=22720*1000 + counter*size
    while True:
        print(">sn"+str(counter))
        print(seqData["246"][i:i+size])
        counter = counter +1
            
        i = i + size
        if i > 23080*1000 or counter > counterStart+25:
            break
        
        
    counter=1
    counterStart=counter
    size=20000
    i=5177*1000 + counter*size
    while True:
        print(">cu"+str(counter))
        print(seqData["tig00001574_pilon_pilon"][i:i+size])
        counter = counter +1
            
        i = i + size
        if i > 5552*1000 or counter > counterStart+25:
            break
    '''
        

    
    '''
    with open('giabpath10k.pickle', 'wb') as handle:
        pickle.dump(paths, handle)
    with open('giabpathsmall.pickle', 'wb') as handle:
        pickle.dump(smallPaths, handle)

    '''
    '''
    with open('giabpath10k.pickle', 'rb') as handle:
        pths = pickle.load(handle)
    with open('giabpathsmall.pickle', 'rb') as handle:
        smallPaths = pickle.load(handle)
    
    paths = []
    for path in pths: 
        if len(path) > 1: paths.append(path)
    
    paths_ = copy.deepcopy(paths)
    paths = copy.deepcopy(paths_)

    '''
    
    with open('CF062_paths2.pickle', 'rb') as handle:
        p1 = pickle.load(handle)
    #with open('CF062_paths.pickle', 'rb') as handle:
        #p2 = pickle.load(handle)
    
    cutoff=0.5e6
    paths=[]
    shortPaths=[]
    #allPaths=[]
    
    unique = []
    dup = []
    
    for p in p1:  # + p2:
        if len(p) < 1: continue
    
        #remove duplicate paths
        p_tup = p[0].rid, p[0].rpos, p[0].rstrand, p[-1].qid, p[-1].qpos, p[-1].qstrand
        if p_tup not in unique:
            unique.append(p_tup)
        else:          
            dup.append(p)
            #dup.append(p1.pop(i))
            continue

        #clean up unnecessary NNNs
        #--------------------------
        if p[0].is_Nfork():
            p.pop(0)
        if p[-1].is_Nfork():
            p.pop()
            
        i=0
        toPop = []
        for fork in p:
            if fork.is_Nfork():
                if p[i-1].after_id() == p[i+1].before_id() and \
                    p[i-1].after_strand() == p[i+1].before_strand() and \
                    p[i-1].is_switch_reference() and p[i+1].is_switch_query() and \
                    p[i-1].after_pos() <= p[i+1].before_pos():
                        toPop.append(i)
            i = i+1

        c=0
        for i in toPop:
            p.pop(i-c)
            c = c +1
    
        #--------------------------

        if path_helper.path_length(p) > cutoff:
            paths.append(p)
        else:
            shortPaths.append(p)
        
    #allPaths = copy.deepcopy(paths + shortPaths)

    param = Parameters()

    refTigIds = list(refData.keys())
    queryTigIds = list(queryData.keys())

    refTigIds.sort(key=lambda x: -lengthData[x])
    queryTigIds.sort(key=lambda x: -lengthData[x])

    import contig_scaffolder as scaffolder
    import path_helper
    
    repeat = 0
    
    while repeat < 2:
        temp = []
        lemp = []
        remp = []
        scaffolds = []
        complete = dict()

        for tigId in refTigIds:
            #tigId = "tig00007642_pilon_pilon"
            if tigId in complete:
                continue
            
            scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, unitigs, lrdict, param.arks, param.LR_THRESHOLD)
            '''
            print('23094029834234')
            print(scaffold[0].rid)
            print('239048392474038')
            for i in range(0, len(paths)):
                for j in range(0, len(paths[i])):
                    if str(paths[i][j].qid) == '172':
                        k = i,j
                        temp.append(k)
            '''            
            paths = paths + leftover
            '''
            for i in range(0, len(leftover)):
                for j in range(0, len(leftover[i])):
                    if str(leftover[i][j].qid) == '172':
                        k = i,j
                        remp.append(k)
            '''
            complete[tigId] = True
            '''
            for i in range(0, len(scaffold)):
                if str(scaffold[i].qid) == '172':
                    lemp.append(i)
                        #print("+++++===+=+")
            '''
            flag = False
            while not flag:
                flag = True
                if scaffold is None or len(scaffold) <= 0:
                    break
     
                if not scaffold[0].is_Nfork():
                    tigId = scaffold[0].rid
                    if tigId not in complete:
                        scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, unitigs, lrdict, param.arks, param.LR_THRESHOLD, scaffold)
                        paths = paths + leftover
                        complete[tigId] = True
                        flag = False
    
                if not scaffold[-1].is_Nfork():
                    tigId = scaffold[-1].rid
                    if tigId not in complete:
                        scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, unitigs, lrdict, param.arks, param.LR_THRESHOLD, scaffold)
                        paths = paths + leftover
                        complete[tigId] = True
                        flag = False
    
            if scaffold is not None:
                scaffolds.append(scaffold)
            
        if repeat == 0:
            paths = scaffolds + shortPaths
        
        repeat = repeat + 1
            
    '''
     with open('scaffolds.pickle', 'wb') as handle:
        pickle.dump(scaffolds, handle)
      
    with open('scaffolds.pickle', 'rb') as handle:
        scaffolds = pickle.load(handle)

    scaffolds_ = copy.deepcopy(scaffolds)
    with open('scaffoldskeep.pickle', 'wb') as handle:
        pickle.dump(keep, handle)

    keep_ = copy.deepcopy(keep)

    scaffolds = copy.deepcopy(scaffolds_)

    keep = copy.deepcopy(keep_)
    scaffolds=keep
    '''
    
    #note: it seems tig00031379_pilon_pilon is misassembled between chr3 and 10 but is fixed? :)
    
    #scaffolds = p1
    
    scaffolds.sort(key=lambda scaffold: path_helper.path_length(scaffold))
    lengths = [path_helper.path_length(scaffold) for scaffold in scaffolds]
    
    keep = []
    
    for i in range(len(scaffolds)):
        scaffold = scaffolds[i]
        shouldKeep=True
        #print(i)
        #if (path_helper.path_length(scaffolds[i]) < 50000): continue
        #if(i<340): continue
        for j in reversed(range(i, len(scaffolds))):
            scaffold2 = scaffolds[j]
            if lengths[j] >= lengths[i] and i != j:
            
                    o1, o2 = path_helper.path_overlap(scaffold, scaffold2, lengthData, source='r')
                    o = max(o1, o2)
                    if(o > 0):
                        print("--------------")
                        print(scaffold.__repr__())
                        print(scaffold2.__repr__())
                        print(path_helper.path_overlap(scaffold, scaffold2, lengthData, printInfo=True, source='r'))
                        print(str(i) + "," + str(j))
                        
                        if o > 0.99:
                            print('deleting:' + scaffold.__repr__())
                            shouldKeep=False
                            break
                        elif o > 0:
                            print('trimming:' + scaffold.__repr__())
                            scaffolder.remove_overlap(scaffold2, scaffold, lengthData, source='r')
                            print('try again:' + scaffold.__repr__())
                            print(path_helper.path_overlap(scaffold, scaffold2, lengthData, printInfo=True, source='r'))
                            if len(scaffold) < 2:
                                shouldKeep=False

        if shouldKeep:
            print(i)
            keep.append(scaffold)
            
            
    '''
    with open('scaffoldskeep.pickle', 'wb') as handle:
        pickle.dump(keep, handle)
      
    with open('scaffoldskeep.pickle', 'rb') as handle:
        keep = pickle.load(handle)

    '''

    keep.sort(key=lambda scaffold: -path_helper.path_length(scaffold))
    
    file="/Users/allen bao/Documents/assembly_data/CF062_hybrid_LR2_strict.fa"
    sourceFile="/Users/allen bao/Documents/assembly_data/hybrid_source.fa"

    with open(file, "w+") as fasta:
        with open(sourceFile, "w+") as fastaSource:
    
            for i in range(len(keep)):
                sequence, source = output.path_to_sequence2(keep[i], seqData)
                fasta.write(">hybrid" + "_" + str(i) +"\n")
                fasta.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
                #fastaSource.write(">hybrid" + "_" + str(i) +"\n")
                #fastaSource.write(re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")
                
            fastaSource.close()
        fasta.close()

if __name__== "__main__":
  main()
  print("done")
  exit()

def checkLR(p1, lrdict):
#checking how often path-switches are found in linked reads (when canu tig changes)
    distrib = dict()
    distrib[0] = 0
    dists = dict()
    dists[0] = list()
    LRtigs = []
    pathTigs = []
    notPathTigs = []
    for p in p1:
        #if path_helper.path_length(p) < 10000000:
            #continue
        for i in range(0, len(p)-1):
            if p[i].rid not in pathTigs:
                pathTigs.append(p[i].rid)
            if p[i+1].rid not in pathTigs:
                pathTigs.append(p[i+1].rid)
            if p[i].rid != p[i+1].rid and p[i].rid != 'NNN' and p[i+1].rid != 'NNN':
                tig1 = p[i].rid
                tig2 = p[i+1].rid
                pos1 = p[i].after_pos()
                
                pos2 = p[i+1].before_pos()
                dist = abs(pos2 - pos1)
                if tig1 in lrdict:
                    a = lrdict[tig1]
                    if tig2 in a:
                        n = lrdict[tig1][tig2]
                        if n not in distrib:
                            distrib[n] = 1
                            dists[n] = list()
                            dists[n].append(dist)
                        else:
                            distrib[n] += 1
                            dists[n].append(dist)
                    else:
                        distrib[0] += 1
                        dists[0].append(dist)
                        #print('tig2: '+tig2)
                        #print(p[i])
                        #print(p[i+1])
                        #print('===============')
                else:
                    distrib[0] += 1
                    dists[0].append(dist)
                    #print('tig1: '+tig1)

                
    distribLR = dict()
    for tig in lrdict:
        if tig not in LRtigs:
            LRtigs.append(tig)
        for tig2 in lrdict[tig]:
            if lrdict[tig][tig2] not in distribLR:
                distribLR[lrdict[tig][tig2]] = 1
            else:
                distribLR[lrdict[tig][tig2]] += 1
    for tig in LRtigs:
        if tig not in pathTigs:
            notPathTigs.append(tig)