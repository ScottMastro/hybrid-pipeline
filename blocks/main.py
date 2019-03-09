import file_handler as reader
import contig_stitcher as stitcher
from parameters import Parameters
import sys
 
#----Parameters-------------------------

param = Parameters()

if len(sys.argv) > 1 :
    summary_file=sys.argv[1]  #input
    block_file = sys.argv[2]  #output
else:
    #prefix="C:/Users/scott/Desktop/novafix/blocks"
    prefix="/home/scott/Dropbox/hybrid-pipeline/blocks"
    summary_file = prefix + "/summary.txt"
    block_file = prefix + "/blocks_new.txt"

#----Parameters-------------------------

def main():
    
    aligndf = reader.parse_alignments(summary_file)   
    #contigs = stitcher.stitch_contigs(aligndf, param)
    contig = stitcher.stitch_contig(aligndf, 0, param)

    print(contig)
  
  
if __name__== "__main__":
  main()
  print("done")
