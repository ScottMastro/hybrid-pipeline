import vcf
import os
import re
import random
import math
import sys
sys.path.append("../..")

import utils.log as logger
import utils.file_handler as io
import utils.fasta_handler as fasta
import utils.longranger_vcf_correct as corrector
import utils.external_tools as tools

#==================================================
# Fetching reads
#==================================================

def fetch_ref_reads(bamFile, region, outdir, fileSuffix=""):
    refAlignments = tools.samtools_fetch(bamFile, region)
    return tools.samtools_write(refAlignments, outdir + "pacbio_reads" + fileSuffix, bamFile)

def fetch_query_reads(bamFile, region, outdir):
    io.delete_file(outdir + "10x_reads", deleteDirTree=True)
    return tools.bam2fastq_10x(bamFile, outdir, dirName="10x_reads", region=region, returnFqList=True)

#==================================================
# Aligning reads
#==================================================

def align_ref(faFile, refReads, outdir, param, outName="ref_aligned"):
    
    bam = outdir + outName + ".pbmm2.bam"
    if os.path.isfile(bam):
        logger.log("pbmm2 alignments found, skipping step.")
        return bam
    
    bam = tools.align_pacbio(faFile, refReads, outdir + outName)
    return bam

def align_query(faFile, queryReads, outdir, param, outName="query_aligned"):
    #note: outName is irrelevant
    
    queryReadsDir = "/".join(queryReads[0].split("/")[:-1])
    
    sampleName = ".".join(faFile.split(".")[:-1]).split("/")[-1] + "-longranger"

    longrangerBam = outdir + sampleName + "/phased_possorted_bam.bam"
    longrangerVCF = outdir + sampleName + "/phased_variants.vcf.gz"

    if not os.path.isfile(longrangerBam) or not os.path.isfile(longrangerVCF):
        longrangerVCF, longrangerBam = tools.align_10x(faFile, queryReadsDir, outdir, svBlacklist=True)
    else:
        print("Long Ranger files found, skipping step")
        return longrangerBam

    #correct bug in longranger 2.2.2 where hets are sometimes given GT = 1|1 
    correctedVCF = re.sub(".gz", "", longrangerVCF)
    corrector.run_correction(longrangerVCF, correctedVCF)
    io.delete_file(longrangerVCF)
    
    longrangerVCF = tools.bgzip(correctedVCF)
    #longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb"))]

    if not param.KEEP_INTERMEDIATE:
        io.delete_file(outdir + "refdata-" + outName, deleteDirTree=True)

    return longrangerBam

#==================================================
# Polishing with reads (all)
#==================================================

def polish_ref(targetFa, region, outdir, param, refReads=None, refAlignments=None, outName="ref_polished"):
        
    polishedFa = outdir + outName + ".consensus.fasta"
    if os.path.isfile(polishedFa):
        print("Arrow polished fasta found, skipping step.")
        return polishedFa
    
    deleteAlignments = False
    if refAlignments is None:
        refAlignments = align_ref(targetFa, refReads, outdir, param, outName)
        deleteAlignments = True

    polishedFa = tools.pacbio_polish(refAlignments, targetFa, outdir + outName, region=region, outputFasta=True)

    if deleteAlignments:
        io.delete_file(refAlignments)
        
    fasta.index_fasta(polishedFa)
    return polishedFa

def polish_query(targetFa, region, outdir, param, queryReads=None, queryAlignments=None, outName="query_polished"):

    polishedFa = outdir + outName + ".fasta"
    if os.path.isfile(polishedFa):
        print("Pilon polished fasta found, skipping step.")
        return polishedFa

    deleteAlignments = False
    if queryAlignments is None:
        qreads1, qreads2 = queryReads
        trimmedr1 = tools.trim_10x_barcode(qreads1)
        queryAlignments = tools.bwa_mem(targetFa, trimmedr1, qreads2, outdir + "bwa_aligned")
        deleteAlignments = True

    pilonFa = tools.pilon_polish(queryAlignments, targetFa, outdir, outputVCF=False)
    polishedFa = outdir + outName + ".fasta"

    io.move_file(pilonFa, polishedFa)
    fasta.index_fasta(polishedFa)
    
    for ext in [".sa", ".pac", ".bwt", ".ann", ".amb", ".ann"]:
        io.delete_file(polishedFa + ext)
        
    if deleteAlignments:
        io.delete_file(queryAlignments)      
        io.delete_file(trimmedr1)
        
    return polishedFa

#==================================================
# Phasing
#==================================================

def phase_region(fa, reads, prefix, region=None, writeBams=True):
    longshotVCF = tools.longshot_genotype(reads, fa, prefix, region, writeBams=writeBams)
    if not writeBams: return longshotVCF
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    return (longshotVCF, bamA, bamB, bamUnphased)

def high_conf_hets(consensusFa, refBam, queryBam, outdir, param):
    
    vcfFile = outdir + "high_confidence_hets.vcf"
    if os.path.isfile(vcfFile):
        logger.log("High confidence heterozygous variants found, skipping step.")
        return vcfFile

    longshotVCF = phase_region(consensusFa, refBam, outdir + "phased", writeBams=False)
    longrangerVCF = "/".join(queryBam.split("/")[:-1]) +  "/phased_variants.vcf.gz" 

    longshotVariants = [record for record in vcf.Reader(open(longshotVCF, "rb" if longshotVCF.endswith("gz") else "r"))]
    longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb" if longrangerVCF.endswith("gz") else "r"))]

    longshotSnps = {v.POS : v for v in longshotVariants if v.is_snp and v.samples[0].is_het}
    longrangerSnps = {v.POS : v for v in longrangerVariants if v.is_snp and v.samples[0].is_het}
    
    commonSnpPos = set()
    commonSnps = []

    for pos in longshotSnps:
        if pos in longrangerSnps and \
           longshotSnps[pos].REF == longrangerSnps[pos].REF and \
           longshotSnps[pos].ALT == longrangerSnps[pos].ALT:
               
               if len(longshotSnps[pos].FILTER) < 1 and len(longrangerSnps[pos].FILTER) < 1:
                   commonSnpPos.add(pos)
                   commonSnps.append(longrangerSnps[pos])

    longrangerOnly = [longrangerSnps[pos] for pos in longrangerSnps if pos not in commonSnpPos]
    highQualityLongranger = [v for v in longrangerOnly if len(v.FILTER) < 1 and v.QUAL > 150]

    longshotOnly = [longshotSnps[pos] for pos in longshotSnps if pos not in commonSnpPos]
    
    #TODO: find a good threshold here for quality
    highQualityLongshot = [v for v in longshotOnly if len(v.FILTER) < 1 and v.QUAL > 15]

    highConfidenceSnps = commonSnps + highQualityLongranger + highQualityLongshot
    
    longrangerNotSnp = [v for v in longrangerVariants if not v.is_snp and v.samples[0].is_het]
    
    highConfidenceIndels=[]
    for v in longrangerNotSnp:
        
        #bialleleic
        if len(v.alleles) > 2 or len(v.ALT) > 1: 
            continue
        #homopolymer
        if v.INFO["POSTHPC"][0] > 1:
            continue
        #filtered
        if len(v.FILTER) > 0:
            continue
        #complex
        if len(v.REF) != 1 and len(v.ALT[0]) != 1:
            continue
        if len(v.REF) > 8 or len(v.ALT[0]) > 8:
            continue
        #low qual
        if v.QUAL < 100:
            continue

        if not v.samples[0].is_het:
            continue
                
        highConfidenceIndels.append(v)

    highConfidenceVariants = highConfidenceSnps + highConfidenceIndels
    
    #vcfFile = outdir + "high_confidence_hets.vcf"
    writer = open(vcfFile, "w+")
    writer.write("##fileformat=VCFv4.2\n")
    header = ["#CHROM", "POS", "ID", "REF", "ALT",
               "QUAL",  "FILTER", "INFO", "FORMAT", "SAMPLE"]
    writer.write("\t".join(header) + "\n")

    for variant in sorted(highConfidenceVariants, key=lambda x: x.POS):
        line = [variant.CHROM, variant.POS, ".", variant.REF, variant.ALT[0], 
                variant.QUAL, "PASS", ".", "GT", "0/1"]       
        writer.write("\t".join([str(x) for x in line]) + "\n")

    writer.close()
    
    return vcfFile

def phase_consensus(consensusFa, highConfVCF, refBam, queryBam, outdir, param):

    #phase het variants with all reads using whatshap
    whatshapVCF = tools.whatshap_phase([refBam, queryBam], highConfVCF, consensusFa,
                                 outdir + "high_confidence",  indels=True, maxCov=17)
    
    refHaploBam = tools.whatshap_haplotag(refBam, whatshapVCF, consensusFa)    
    queryHaploBam = tools.whatshap_haplotag(queryBam, whatshapVCF, consensusFa)
    
    if not param.KEEP_INTERMEDIATE:
        io.delete_file(refBam)
    
    return refHaploBam, queryHaploBam        

#==================================================
# Polishing with reads (haplotype)
#==================================================

def split_reads(bamFile):
    reads = tools.samtools_fetch(bamFile)
    readDict = {"1" : [], "2" : [], "0" : []}
    
    for read in reads:
        try:
            hp = read.get_tag("HP")
            if hp == 1: readDict["1"].append(read)
            elif hp == 2: readDict["2"].append(read)
            else: readDict["0"].append(read)
        except:
            readDict["0"].append(read)
    return readDict

def haplotype_polish_query(hap, hapFa, haploBam, outdir, param, realign=True):
    hap = str(hap)
    reads = split_reads(haploBam)

    assignedReads = reads[hap]
    idSet = {r.query_name for r in assignedReads}

    for r in reads["1" if hap == "2" else "2"] + reads["0"]:
        if r.query_name in idSet:
            assignedReads.append(r)
    
    assignedBam = tools.samtools_write(assignedReads, outdir + hap + "_polish_input_query", haploBam)
    
    if realign:
        dirname = "10x_hap" + hap
        r1,r2 = tools.bam2fastq_10x(assignedBam, outdir, dirName=dirname, returnFqList=True)
        r1Trimmed = tools.trim_10x_barcode(r1)
        tools.bwa_index(hapFa)
        realignedBam = tools.bwa_mem(hapFa, r1Trimmed, r2, outdir + "bwa_aligned_hap" + hap)
        
        io.delete_file(outdir+dirname, deleteDirTree=True)
        io.delete_file(assignedBam)
        assignedBam = realignedBam
        for ext in [".sa", ".pac", ".bwt", ".ann", ".amb", ".ann"]:
            io.delete_file(hapFa + ext)

    hapPolishedFa = tools.pilon_polish(assignedBam, hapFa, outdir, prefix="hap" + hap + ".consensus", 
                                       changeFile=True, snpIndelOnly=True)
    fasta.index_fasta(hapPolishedFa)

    #helper.delete_file(assignedBam)
    return hapPolishedFa

def haplotype_polish_ref(hap, hapFa, haploBam, outdir, param, realign=True):

    hap = str(hap)
    reads = split_reads(haploBam)

    random.shuffle(reads["0"])
    half = math.ceil(len(reads["0"])/2)
    assignedReads = reads[hap] + reads["0"][:half]
    assignedBam = tools.samtools_write(assignedReads, outdir + hap + "_polish_input_ref", haploBam)

    if realign:
        realignedBam = tools.align_pacbio(hapFa, assignedBam, outdir + hap + "_polish_input_ref.realigned")
        io.delete_file(assignedBam)
        assignedBam = realignedBam
        
    hapPolishedFa = tools.pacbio_polish(assignedBam, hapFa, outdir + "hap" + hap, outputFasta=True)
    
    io.delete_file(assignedBam)

    return hapPolishedFa

    
def clean_directory(outdir):
    
    #delete reads
    io.delete_file(outdir + "10x_reads", deleteDirTree=True)
    io.delete_file(outdir + "pacbio_reads.bam")
