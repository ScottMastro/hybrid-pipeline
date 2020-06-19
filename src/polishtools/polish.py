import os
import time
import sys
sys.path.append("../..")

import utils.log as logger
import polishtools.polish_implementation as impl
import utils.file_handler as io
import utils.fasta_handler as fasta

import structures.region as rgn
import utils.external_tools as tools
    
def polish_contig(tigId, outdir, seqData, lengthData, param):
    cd = os.getcwd()
    os.chdir(outdir)
    
    startTime = time.time()
    
    seqNames = ["consensus", "hap1", "hap2"]

    region = rgn.SimpleRegion(tigId, 0, lengthData[tigId])

    #==================================================
    # Fetch reads
    #==================================================
    
    logger.log("Fetching reads from BAM files...")
    refReads = impl.fetch_ref_reads(param.REF_ALIGNED_READS, region, outdir)
    queryReads = impl.fetch_query_reads(param.QUERY_ALIGNED_READS, region, outdir) 
    hybridFa = fasta.write_fasta(outdir + "unpolished.fasta", {tigId : seqData[tigId]})

    #==================================================
    # Polishing with reads (all)
    #==================================================
    
    logger.log("Running initial polish with all reads.")
    logger.log("Polish sequence using reference reads...", indent=1)
    refPolishedFa = impl.polish_ref(hybridFa, region, outdir, param, \
                                    refAlignments=param.REF_ALIGNED_READS, outName="ref_polished")
    logger.log("Polish sequence using query reads...", indent=1)
    queryPolishedFa = impl.polish_query(refPolishedFa, region, outdir, param, \
                                        queryReads=queryReads, outName="consensus")

    consensusFa = fasta.rename_single_fasta(queryPolishedFa, seqNames[0], toUpper=True)

    #==================================================
    # Split haplotypes
    #==================================================

    logger.log("Aligning reads back to consensus.")
    logger.log("Aligning reference reads...", indent=1)
    refBam = impl.align_ref(consensusFa, refReads, outdir, param, "ref_to_consensus")
    logger.log("Aligning query reads...", indent=1)
    queryBam = impl.align_query(consensusFa, queryReads, outdir, param)

    logger.log("Getting heterozygous variants and phasing reads.")
    highConfVCF = impl.high_conf_hets(consensusFa, refBam, queryBam, outdir, param)
    refHaploBam, queryHaploBam = impl.phase_consensus(consensusFa, highConfVCF, refBam, queryBam, outdir, param)

    #==================================================
    # Polishing with reads (haplotypes)
    #==================================================

    hap1Fa, hap2Fa = consensusFa, consensusFa
    finalFa1, finalFa2 = outdir + "hap1.fasta", outdir + "hap2.fasta"
    
    if io.file_exists(finalFa1) and os.path.isfile(finalFa2):
        logger.log("Polished haplotypes found, skipping step.")
    else:
    
        niter = 1
        realign = False
        while niter > 0:
            hap1Fa = impl.haplotype_polish_ref(1, hap1Fa, refHaploBam, outdir, param, realign)
            hap2Fa = impl.haplotype_polish_ref(2, hap2Fa, refHaploBam, outdir, param, realign)
            realign=True
            
            if param.KEEP_INTERMEDIATE:
                #optional back align step
                hap1Fa = fasta.rename_single_fasta(hap1Fa, seqNames[1])
                hap2Fa = fasta.rename_single_fasta(hap2Fa, seqNames[2])
                tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refpolished")
                tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refpolished")
            
            hap1Fa = impl.haplotype_polish_query(1, hap1Fa, queryHaploBam, outdir, param, realign)
            hap2Fa = impl.haplotype_polish_query(2, hap2Fa, queryHaploBam, outdir, param, realign)
        
            if param.KEEP_INTERMEDIATE:
                #optional back align step
                hap1Fa = fasta.rename_single_fasta(hap1Fa, seqNames[1])
                hap2Fa = fasta.rename_single_fasta(hap2Fa, seqNames[2])
                tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refquerypolished")
                tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refquerypolished")
    
            niter -=1
        
        io.move_file(hap1Fa, finalFa1)
        io.move_file(hap2Fa, finalFa2)

        finalFa1 = fasta.rename_single_fasta(finalFa1, seqNames[1])
        finalFa2 = fasta.rename_single_fasta(finalFa2, seqNames[2])
        
        logger.log("Total elasped time: " + str(round(time.time() - startTime)))
    
        if not param.KEEP_INTERMEDIATE:
            impl.clean_directory(outdir)
            
        os.chdir(cd)
        return finalFa1, finalFa2