import os, re
import utils.log as logger
import utils.fasta_handler as fasta
import polishtools.polish as polisher

def main(param):
    
    #==================================================
    # Read in data
    #==================================================

    if not os.path.exists(param.OUTPUT_DIR): os.mkdir(param.OUTPUT_DIR)
    
    logger.Logger(clean=True, level=param.VERBOSE)
    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)

    logger.log("Reading fasta...")
    
    if param.TARGET_CONTIG is None:
        seqData = fasta.read_fasta(param.FASTA)
    else:
        target = str(param.TARGET_CONTIG)
        seqData = {target : fasta.fasta_fetch(param.FASTA, target)}
        
    tigIds = list(seqData.keys())
    lengthData = {x : len(seqData[x]) for x in tigIds}
    tigIds.sort(key=lambda x: -lengthData[x])

    #==================================================
    # Polish
    #==================================================
    
    for tigId in tigIds:
        logger.log("Polishing " + tigId + "...")

        cleanId = re.sub(r'[\\/*?:"<>|]', "_", str(tigId))
        outdir = param.OUTPUT_DIR + "/" + cleanId + "/"
        if not os.path.exists(outdir): os.mkdir(outdir)
        
        polisher.polish_contig(tigId, outdir, seqData, lengthData, param)