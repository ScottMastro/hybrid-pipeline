import sys
sys.path.append('./align')
import contig_welder as welder
sys.path.append('./blocks')
import contig_stitcher as stitcher

import file_handler as reader
import contig_scaffolder as scaffolder
import contig_output as output
from path_helper import clean_path
import path_helper

from parameters import Parameters
from Bio import SeqIO
import sys
import gzip
import pickle
import copy 

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
    
    cutoff=1e6
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
            paths = scaffolds + smallPaths
        
        repeat = repeat + 1
            
            
            
            
            
            
            
            
            
            
            
            
    canuTigs = dict()       
    for tigId in refData.keys():
        print(tigId)
        contig = stitcher.stitch_r(aligndf, tigId, param )
        if contig is None: continue
        canuTigs[tigId] = contig

            
    for tigId in queryTigIds:
        if tigId in complete:
            continue
        
        canuTigs['tig00001060_pilon_pilon'].mblocks[0].qid
        
        
        scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param)
        leftovers = leftovers + leftover
        complete[tigId] = True
        
        flag = False
        while not flag:
            flag = True
            if scaffold is None or len(scaffold) <= 0:
                break
 
            tigId = scaffold[0].rid
            if tigId not in complete:
                scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, scaffold)
                leftovers = leftovers + leftover
                complete[tigId] = True
                flag = False

            tigId = scaffold[-1].rid
            if tigId not in complete:
                scaffold, paths, leftover = scaffolder.scaffold(paths, tigId, lengthData, param, scaffold)
                leftovers = leftovers + leftover
                complete[tigId] = True
                flag = False

        if scaffold is not None:
            scaffolds.append(scaffold)

            
            
            
            
            
            

    
    tg = [len(get_tig_ids(x, source='r')) for x in scaffolds]
    ln = np.array([path_length(x) for x in scaffolds])
    
    def fn (path, end=True):
        x = path [-1] if end else path[0]
        if(x.after_id() == x.rid):
            pos = x.after_pos()
            if x.after_strand() == -1:
                pos = lengthData[str(x.rid)] - pos
        else:
            return -0.1
        print(str(pos) + " / " + str(lengthData[str(x.rid)]))
        return pos/lengthData[str(x.rid)]
        
    end = [fn(x, True) for x in scaffolds]

            
    
    for scaffold in leftovers:
        if scaffold[0].rid != scaffold[-1].rid:
            print(scaffold[0].rid + " - " + scaffold[-1].rid)

    
    pacbioTigs = set()
    for path in paths:
        for fork in path:
            pacbioTigs.add(fork.rid)
    
    output.output_contigs(lst, seqData, "novachr3.fasta")
    output.output_contigs(pacbioTigs, seqData, "canuchr3.fasta")
    output.path_to_sequence(paths[0], seqData, "hybridchr3.fasta")

  
if __name__== "__main__":
  main()
  print("done")

def pickler():
    with open('paths2.pickle', 'rb') as handle:
        paths = pickle.load(handle)
        
    
    
    with open('pathsx.pickle', 'wb') as handle:
        pickle.dump(paths_, handle)
    
    '''
    novaTigs = dict()
    canuTigs = dict()       
    for tigId in refData.keys():
        print(tigId)
        mblocks = stitcher.stitch_r(aligndf, tigId, param, mblockOnly=True )
        if mblocks is None or len(mblocks) < 0: continue
        canuTigs[tigId] = mblocks
    for tigId in queryData.keys():
        print(tigId)
        tigId = int(tigId)
        mblocks = stitcher.stitch(aligndf, tigId, param, mblockOnly=True )
        if mblocks is None or len(mblocks) < 0: continue
        novaTigs[tigId] = mblocks
        
    match = dict()
    
    for tigId in canuTigs:
        for mblock in canuTigs[tigId]:
            match.setdefault(tigId, []).append(mblock.qid)
    for tigId in novaTigs:
        for mblock in novaTigs[tigId]:
            match.setdefault(tigId, []).append(mblock.rid)

    x=[]
    y=[]
    for tigId in novaTigs:
        if tigId in match:
            x.append(len(match[tigId]))
            y.append(lengthData[str(tigId)])

    '''
        
