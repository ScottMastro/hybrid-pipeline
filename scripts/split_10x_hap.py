import sys
import argparse
import pysam
import os

def add_slash(string):
    if len(string) == 0: return "./"
    return string + ("" if string.endswith("/") else "/")

inputBam=sys.argv[1]
outdir=add_slash(sys.argv[2])
prefix = sys.argv[3] if len(sys.argv) > 3 else ""

def samtools_write(alignments, outName, headerBam, makeUnique=False):

    # remove redundant reads:
    if makeUnique:
        readDict = {read.qname : read for read in alignments}
        uniqueReads = [readDict[k] for k in readDict]
        alignments = uniqueReads
        
    alignments = sorted(alignments, key=lambda x: (x.is_unmapped, x.rname, x.pos))

    header = pysam.AlignmentFile(headerBam, "rb")
    bamFile = pysam.AlignmentFile(outName, "wb", template=header)
    for alignment in alignments:
         bamFile.write(alignment)

    bamFile.close()
    header.close()

    return outName

def samtools_fetch(bamFile, region=None, unaligned=True):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    alignments = [x for x in samfile.fetch(until_eof=unaligned)]
    return alignments
    
def split_reads(bamFile, outdir):
    reads = samtools_fetch(bamFile)
    readDict = {"1" : [], "2" : [], "0" : []}
    
    for read in reads:
        try:
            hp = read.get_tag("HP")
            if hp == 1: readDict["1"].append(read)
            elif hp == 2: readDict["2"].append(read)
            else: readDict["0"].append(read)
        except:
            readDict["0"].append(read)

    phasedIdSet = set()
    for hap in ["1","2"]:

        assignedReads = [x for x in readDict[hap]]
        idSet = {r.query_name for r in assignedReads}
        phasedIdSet = phasedIdSet.union(idSet)
        
        for r in readDict["1" if hap == "2" else "2"] + readDict["0"]:
            if r.query_name in idSet: assignedReads.append(r)
    
        assignedBam = samtools_write(assignedReads, outdir + prefix + "haplotype" + hap + ".bam", bamFile)

    unphasedReads = []
    for r in readDict["0"]:
        if r.query_name not in phasedIdSet: unphasedReads.append(r)
    
    assignedBam = samtools_write(unphasedReads, outdir + prefix + "unassigned.bam", bamFile)
    
split_reads(inputBam, outdir)

