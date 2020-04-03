import os
import utils.parameters as parameters
import utils.log as logger
import utils.fasta_handler as fasta

import polish.polish_region as polisher

def main():
    
    #--------------------------------------
    # READ IN DATA
    #--------------------------------------

    param = parameters.get_parameters_polish()
    if not os.path.exists(param.OUTPUT_DIR): os.mkdir(param.OUTPUT_DIR)
    
    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)
    logger.Logger(level=param.VERBOSE, wait=param.WAIT)

    logger.Logger().out("Reading fasta...")
    seqData = fasta.read_fasta(param.FASTA)
    tigIds = list(seqData.keys())
    lengthData = {x : len(seqData[x]) for x in tigIds}
    tigIds.sort(key=lambda x: -lengthData[x])

    #logger.Logger().out("Reading alignment data...")
    #alignDict = io.parse_alignments(param.SUMMARY)

    #--------------------------------------
    # POLISH
    #--------------------------------------

    for tigId in tigIds:
        polisher.polish_contig(tigId, seqData, lengthData, param)
        
if __name__== "__main__":
  main()
  exit()
    
        
        
        
