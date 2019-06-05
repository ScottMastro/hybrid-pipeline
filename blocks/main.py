import sys
sys.path.append('./align')
import contig_welder as welder
sys.path.append('./blocks')
import contig_stitcher as stitcher

import file_handler as reader
import contig_scaffolder as scaffolder
import contig_output as output
import path_helper

from parameters import Parameters
from Bio import SeqIO
import sys
import gzip
import pickle
import copy 
import re
#----Parameters-------------------------

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
    novafa = sys.argv[3]      #supernova fasta
    canufa = sys.argv[4]      #canu fasta

else:
    
    prefix1="/home/scott/Dropbox/hybrid-pipeline/blocks"
    #prefix2="/media/scott/HDD/sickkids"
    prefix2="/media/scott/Rotom/assembly_data"

    summary_file = prefix1 + "/summary2.txt"
    novafa = prefix2 + "/CF062/OSK7121_03C_supernova.pseudohap2.1.fasta"
    canufa = prefix2 + "/CF062/CF062B2D.contigs.fasta.PILON2.fasta"

    
    summary_file = prefix1 + "/summary_giab.txt"
    novafa = prefix2 + "/NA24385/NA24385_supernova.pseudohap2.1.fasta"
    canufa = prefix2 + "/NA24385/HG002_NA24385_son_57X.contigs.fasta.PILON2"



#----Parameters-------------------------

def read_fa(fa):
    if(fa[-2:] == "gz"):
        with gzip.open(fa, "rb") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fa, "fasta"))
        
    return dict(zip([r.id for r in records], [str(r.seq) for r in records]))

def main():

    print("Reading Canu fasta")
    refData = read_fa(canufa)
    print("Reading Supernova fasta")
    queryData = read_fa(novafa)    
    
    seqData = {}
    seqData.update(refData)
    seqData.update(queryData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    
    param = Parameters()
    aligndf = reader.parse_alignments(summary_file)   
    
    paths = []
    smallPaths = []
    
    for tigId in queryData.keys():
        print(tigId)
        tigId = int(tigId)
        #if(tigId < 100612): 
        #   continue
        contig = stitcher.stitch(aligndf, tigId, param)
        if contig is None: continue
        #output.plot_identity(contig, outputPath="/home/scott/Dropbox/hybrid-pipeline/blocks/idplots/")
    
        path = welder.weld(contig, seqData, lengthData, param)
        
        if len(path)>1:
            if lengthData[str(tigId)] > 10000:
                paths.append(path)
            else:
                smallPaths.append(path)




        #cleanPath = clean_path(path, lengthData, param)

        
    
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
    
    
    
    with open('giabpath10k.pickle', 'rb') as handle:
        p1 = pickle.load(handle)
    with open('giabpathsmall.pickle', 'rb') as handle:
        p2 = pickle.load(handle)
    
    cutoff=0.5e6
    paths=[]
    shortPaths=[]
    allPaths=[]
    for p in p1 + p2:
        if len(p) < 1: continue
    
    
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
        
    allPaths = copy.deepcopy(paths + shortPaths)

    param = Parameters()

    refTigIds = list(refData.keys())
    queryTigIds = list(queryData.keys())

    refTigIds.sort(key=lambda x: -lengthData[x])
    queryTigIds.sort(key=lambda x: -lengthData[x])

    repeat = 0
    
    while repeat < 2:
        
        scaffolds = []
        complete = dict()

        for tigId in refTigIds:
            if tigId in complete:
                continue
            
            scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param)
            paths = paths + leftover
            complete[tigId] = True
            
            flag = False
            while not flag:
                flag = True
                if scaffold is None or len(scaffold) <= 0:
                    break
     
                if not scaffold[0].is_Nfork():
                    tigId = scaffold[0].rid
                    if tigId not in complete:
                        scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, scaffold)
                        paths = paths + leftover
                        complete[tigId] = True
                        flag = False
    
                if not scaffold[-1].is_Nfork():
                    tigId = scaffold[-1].rid
                    if tigId not in complete:
                        scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, scaffold)
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
    
    scaffolds.sort(key=lambda scaffold: path_helper.path_length(scaffold))
    lengths = [path_helper.path_length(scaffold) for scaffold in scaffolds]
    
    keep = []
    
    for i in range(len(scaffolds)):
        scaffold = scaffolds[i]
        shouldKeep=True
        print(i)
        if(i<340): continue
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
            keep.append(scaffold)
            
            
    '''
    with open('scaffoldskeep.pickle', 'wb') as handle:
        pickle.dump(keep, handle)
      
    with open('scaffoldskeep.pickle', 'rb') as handle:
        keep = pickle.load(handle)

    '''

    keep.sort(key=lambda scaffold: -path_helper.path_length(scaffold))
    
    file="/media/scott/Rotom/assembly_data/NA24385/hybrid.fa"
    sourceFile="/media/scott/Rotom/assembly_data/NA24385/hybrid_source.fa"

    with open(file, "w+") as fasta:
        with open(sourceFile, "w+") as fastaSource:
    
            for i in range(len(keep)):
                print(i)
                sequence, source = output.path_to_sequence2(keep[i], seqData)
                fasta.write(">hybrid" + "_" + str(i) +"\n")
                fasta.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
                fastaSource.write(">hybrid" + "_" + str(i) +"\n")
                fastaSource.write(re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")
                
            fastaSource.close()
        fasta.close()
    

if __name__== "__main__":
  main()
  print("done")
