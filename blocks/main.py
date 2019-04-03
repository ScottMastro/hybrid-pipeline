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

    lst = [0,3,4,5,6,8,9,10,11,12,13,14,15,27]
    for qid in lst: #queryData.keys():
        print(qid)
        contig = stitcher.stitch(aligndf, qid, param)
        path = welder.weld(contig, seqData, lengthData, param)   
        paths.append(path)
        
    paths = scaffolder.scaffold(paths,"tig00000569_pilon_pilon", lengthData, param)
    paths = scaffolder.scaffold(paths,"tig00000002_pilon_pilon", lengthData, param)

    print(contig)
  
  
if __name__== "__main__":
  main()
  print("done")
