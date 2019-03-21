import file_handler as reader
import contig_stitcher as stitcher

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
    #prefix="/media/scott/HDD/sickkids/blocks/"
    summary_file = prefix + "/summary.txt"
    block_file = prefix + "/blocks_new.txt"
    
    prefix = "/media/scott/Rotom/assembly_data/CF062/"
    #prefix="/media/scott/HDD/sickkids/"

    novafa = prefix + "OSK7121_03C_supernova.pseudohap2.2.fasta"
    canufa = prefix + "CF062B2D.contigs.fasta.PILON2.fasta"

#----Parameters-------------------------

def read_fa(fa):
    if(fa[-2:] == "gz"):
        with gzip.open(fa, "rb") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fa, "fasta"))
        
    return dict(zip([r.id for r in records], [r.seq for r in records]))

def main():
    print("Reading Canu fasta")
    refData = read_fa(canufa)
    print("Reading Supernova fasta")
    queryData = read_fa(novafa)    
    
    seqData = {}
    seqData.update(refData)
    seqData.update(queryData)
    lengthData = {x : len(str(seqData[x])) for x in seqData.keys()}
    
    aligndf = reader.parse_alignments(summary_file)   
    #contigs = stitcher.stitch_contigs(aligndf, param)
    contig = stitcher.stitch_contig(aligndf, 6, param)
    #contig = fastener.fasten_contig(contig, seqData, lengthData, param)   

    print(contig)
  
  
if __name__== "__main__":
  main()
  print("done")
