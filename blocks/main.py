import file_handler as reader
import contig_stitcher as stitcher
import contig_welder as welder
import contig_scaffolder as scaffolder
import contig_output as output

from parameters import Parameters
from Bio import SeqIO
import sys
import gzip

#----Parameters-------------------------

param = Parameters()

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
    novafa = sys.argv[3]      #supernova fasta
    canufa = sys.argv[4]      #canu fasta

else:
    prefix="/home/scott/Dropbox/hybrid-pipeline/blocks"
    summary_file = prefix + "/summary.txt"
    
    prefix = "/media/scott/Rotom/assembly_data/CF062/"
    #prefix="/media/scott/HDD/sickkids/"

    novafa = prefix + "OSK7121_03C_supernova.pseudohap2.1.fasta"
    canufa = prefix + "CF062B2D.contigs.fasta.PILON2.fasta"

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
    
    aligndf = reader.parse_alignments(summary_file)   
        
    paths = []

    lst = [0,3,4,5,6,8,9,10,11,12,13,14,15]
    
    for qid in queryData.keys(): #in lst: 
        print(qid)
        qid = int(qid)
        #if(qid < 73003): 
        #    continue
        contig = stitcher.stitch(aligndf, qid, param)
        if contig is None: continue
        path = welder.weld(contig, seqData, lengthData, param)
        paths.append(path)
    
    
    import copy 
    #paths_ = copy.deepcopy(paths)
    paths = copy.deepcopy(paths_)
    
    leftovers = []
    scaffolds = []

    complete = dict()
    refContigs = list(refData.keys())
    refContigs.sort(key=lambda x: -lengthData[x])
    
    for tigId in refContigs:
        if tigId in complete:
            continue
        
        scaffold, paths, leftover = scaffolder.scaffold2(paths, tigId, lengthData, param)
        leftovers = leftovers + leftover
        complete[tigId] = True
        
        flag = False
        while not flag:
            flag = True
            if scaffold is None or len(scaffold) <= 0:
                break
 
            tigId = scaffold[0].rid
            if tigId not in complete:
                scaffold, paths, leftover = scaffolder.scaffold2(paths, tigId, lengthData, param, scaffold)
                leftovers = leftovers + leftover
                complete[tigId] = True
                flag = False

            tigId = scaffold[-1].rid
            if tigId not in complete:
                scaffold, paths, leftover = scaffolder.scaffold2(paths, tigId, lengthData, param, scaffold)
                leftovers = leftovers + leftover
                complete[tigId] = True
                flag = False

        if scaffold is not None:
            scaffolds.append(scaffold)
    
    
    
    
    
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
