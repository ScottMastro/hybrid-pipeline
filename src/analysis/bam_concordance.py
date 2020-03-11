# Original script from PacBio
# https://github.com/PacificBiosciences/hg002-ccs
# Modified by Scott Mastromatteo

#!/usr/bin/env python2
"""Measure concordance of read alignments to a reference genome"""

from __future__ import division
import argparse
import sys
import collections
import math
import concordance_helper as ch

try: pysam 
except NameError: import pysam


class AlignmentStats:
    def __init__(self, alignment):
        self.alignment = alignment
        self.matchBp = 0
        self.mismatchBp = 0
        self.nonHpInsertionBp = 0
        self.nonHpDeletionBp = 0
        self.hpInsertionBp = 0
        self.hpDeletionBp = 0
        self.ignoredBp = 0
        
        self.discordantDict = dict()
        
    def add_mismatch(self, rpos, length):
        self.discordantDict[rpos] = (ch.DIFF, length)
        self.mismatchBp += 1

    def add_deletion(self, rpos, length, homopolymer=False):
        self.discordantDict[rpos] = (ch.DEL, length)
        if homopolymer: self.hpDeletionBp += 1
        else: self.nonHpDeletionBp += 1
            
    def add_insertion(self, rpos, length, homopolymer=False):
        self.discordantDict[rpos] = (ch.INS, length)
        if homopolymer: self.hpInsertionBp += 1
        else: self.nonHpInsertionBp += 1

    @property
    def insertionBp(self):
        return self.nonHpInsertionBp + self.hpInsertionBp

    @property
    def deletionBp(self):
        return self.nonHpDeletionBp + self.hpDeletionBp

    @property
    def errorsBp(self):
        return self.mismatchBp + self.deletionBp + self.insertionBp

    @property
    def hcReadLength(self): # length of read that fell within high-confidence, non-variant positions
        return self.mismatchBp + self.matchBp + self.deletionBp

    @property
    def concordance(self): # M/(M+X+D+I) = read/reference symmetric concordance
        return 1. * self.matchBp / (self.matchBp + self.mismatchBp + self.deletionBp + self.insertionBp)

    @property
    def concordanceQv(self): # concordance in Phred-scale, capped by hcReadLength
        maxQv = 10*math.log10(self.hcReadLength+1) # QV capped by read length with pseudocount of 1
        qv = -10*math.log10(1-self.concordance) if self.errorsBp else maxQv
        return min(maxQv, qv)

    def check_discordant(self, rpos):
        if rpos in self.discordantDict:
            return self.discordantDict[rpos]
        return None
    def is_singleton(self, rpos):
        if rpos in self.isSingleton:
            return self.isSingleton[rpos]
        return None

def measure_concordance(aln, chrom, reffasta, hcRegions, hcVariants):
    """Measure the concordance of the alignment to the reference."""

    # Walk through each CIGAR block in the alignment
    stats = AlignmentStats(aln)

    rpos = aln.reference_start
    qpos = 0
    qseq = aln.query_alignment_sequence
    hcRegionsIx = 0
    for op, oplen in aln.cigartuples:
        
        global SHOW_MATCH_WARNING
        if op == ch.MATCH:
            if not SHOW_MATCH_WARNING:
                print("""Warning: BAM must be aligned with --eqx, "M" CIGAR operation not supported.""")
                print("""Treating all M as correct match ("=") from now on...""")
                SHOW_MATCH_WARNING = True
            if op == 0: op = ch.EQUAL

        if op in ch.REFOPS:
            # record reference-consuming operations
            for i in range(oplen):
                # only consider operations within the high-confidence regions
                while hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromEnd <= rpos:
                    hcRegionsIx += 1
                if hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromStart <= rpos:
                    
                    if rpos in hcVariants:
                        stats.ignoredBp += 1
                    elif op == ch.EQUAL:
                        stats.matchBp += 1
                    elif op == ch.DIFF:
                        stats.add_mismatch(rpos, oplen)
                    elif op == ch.DEL:
                        if (rpos-1) in hcVariants or (rpos-i) in hcVariants or (rpos-i-1) in hcVariants:
                            stats.ignoredBp += 1
                        else:
                            bases = reffasta.fetch(chrom, rpos-1, rpos+2).upper()
                            isHomopolymer = bases[0] == bases[1] or bases[1] == bases[2] 
                            stats.add_deletion(rpos, oplen, isHomopolymer)
                            
                rpos += 1
                qpos += 1 if op in ch.QUERYOPS else 0
        elif op == ch.INS:
            # Only consider insertions within the high-confidence regions.
            while hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromEnd <= rpos:
                hcRegionsIx += 1
            if hcRegionsIx < len(hcRegions) and hcRegions[hcRegionsIx].chromStart <= rpos:
                if rpos in hcVariants or (rpos-1) in hcVariants:
                    stats.ignoredBp += oplen
                else:
                    bases = reffasta.fetch(chrom, rpos, rpos+2).upper()
                    isHomopolymer = qseq[qpos:qpos+oplen] in (oplen*bases[0], oplen*bases[1])
                    stats.add_insertion(rpos, oplen, isHomopolymer)

            qpos += oplen

    return stats

def measure_concordance_region(alignments, startPos, endPos, svRange, genomecsv, svcsv):
    
    discordancy = dict()
    lengthDict = dict()
    totalCount = dict()
    positions = dict()
    
    start, end = startPos-svRange, endPos+svRange+1
    for rpos in range(start, end):
        discordancy[rpos] = []
        lengthDict[rpos] = []
        totalCount[rpos] = 0

    for aln, alnStats in alignments:
        s = max(aln.reference_start,start)
        e = min(aln.reference_end,end)
        
        for rpos in range(s, e+1): totalCount[rpos] +=1

        for (rpos,(discordType, length)) in alnStats.discordantDict.items():
            if rpos >= start and rpos <= end:
                discordancy[rpos].append(discordType)
                lengthDict[rpos].append(length)
                positions[rpos] = True

    start, end = startPos, endPos+1

    for rpos in positions.keys():
        discordantList = discordancy[rpos]
        depth = totalCount[rpos]
        lengths = lengthDict[rpos]
        
        if len(discordantList) > 3 and depth > 12:
            
            diffCount = sum([discord == ch.DIFF for discord in discordantList])
            delCount = sum([discord == ch.DEL for discord in discordantList])
            insCount = sum([discord == ch.INS for discord in discordantList])
            singletonCount = sum([l == 1 for l in lengths])
    
            

            genomecsv.write("%s,%d,%d,%d,%d,%d,%d,%d\n" % (aln.reference_id, rpos, 
                totalCount[rpos], len(discordantList), diffCount, delCount, insCount,
                singletonCount))

        
        
        
    for (rpos, discordantList), (_, singletonList) in \
        zip(discordancy.items(), lengthDict.items()):
        
            
        if len(discordantList) > 3 and totalCount[rpos] > 12:
            diffCount = sum([discord == ch.DIFF for discord in discordantList])
            delCount = sum([discord == ch.DEL for discord in discordantList])
            insCount = sum([discord == ch.INS for discord in discordantList])
            singletonCount = sum(singletonList)

            genomecsv.write("%s,%d,%d,%d,%d,%d,%d,%d\n" % (aln.reference_id, rpos, 
                totalCount[rpos], len(discordantList), diffCount, delCount, insCount,
                singletonCount))
        




class Bed:
    def __init__(self, line):
        cols = line.rstrip("\n").split()
        self.chrom = cols[0]
        self.chromStart = int(cols[1])
        self.chromEnd = int(cols[2])

class Vcf:
    def __init__(self, chrom, chromStart, chromEnd, ref, alt):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.ref = ref
        self.alt = alt

    @staticmethod
    def LineToVcfs(line):
        cols = line.rstrip("\n").split()

        chrom = cols[0]
        chromStart = int(cols[1]) - 1
        ref = cols[3]
        chromEnd = chromStart + len(ref)
        result = []
        for alt in cols[4].split(","):
            result.append(Vcf(chrom, chromStart, chromEnd, ref, alt))

        return result

class DevNullFile:
    def __init__(self):
        pass

    def close(self):
        pass

    def write(self, s):
        pass

def parseParam(args, dummy=True):
    """Usage"""
    argparser = argparse.ArgumentParser("""Measure concordance of read alignments to a reference genome""",
        epilog="""
    CONCORDANCE
        Measure the concordance of alignments to the reference genome.
        Consider only non-variant positions in high-confidence regions.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    
    if dummy:
        DEFAULT_DIR = "/media/scott/HDD/sickkids/NA24385/"
        DEFAULT_FA = DEFAULT_DIR + "hg38.fa"
        DEFAULT_BAM = DEFAULT_DIR + "HG002.15kb.ccs.bam"
        DEFAULT_READ_CSV = DEFAULT_DIR + "read_concordance_ccs.csv"
        DEFAULT_GENOME_CSV = DEFAULT_DIR + "genome_concordance_ccs.csv"
        DEFAULT_SV_CSV = DEFAULT_DIR + "genome_concordance_sv_ccs.csv"

        #DEFAULT_DIR = "/media/scott/Rotom/assembly_data/hg002/"
        #DEFAULT_FA = DEFAULT_DIR + "hybrid.fa"
        #DEFAULT_BAM = DEFAULT_DIR + "hybrid_pacbio_ccs_15k.mm2.bam"
        #DEFAULT_READ_CSV = DEFAULT_DIR + "read_concordance_ccs.csv"
        #DEFAULT_GENOME_CSV = DEFAULT_DIR + "genome_concordance_ccs.csv"




        chrom= None #"hybrid_10"
        
        args = args + [DEFAULT_FA, DEFAULT_BAM, DEFAULT_READ_CSV, DEFAULT_GENOME_CSV, DEFAULT_SV_CSV] + \
            (["--chrom", chrom] if chrom is not None else [])
        
    argparser.add_argument("reffasta", metavar="ref.fasta")
    argparser.add_argument("inbam", metavar="in.bam")
    argparser.add_argument("readcsv", metavar="read.csv", help="Per aligment concordance statistics, CSV")
    argparser.add_argument("genomecsv", metavar="genome.csv", help="Concordance across reference genome, CSV")
    argparser.add_argument("svcsv", metavar="sv.csv", help="Structural differences across reference genome, CSV")

    argparser.add_argument("--hcregions", metavar="hcregions.bed.gz", help="High-confidence regions, BED.gz (.tbi required)")
    argparser.add_argument("--hcvariants", metavar="hcvariants.vcf.gz", help="High-confidence variants, VCF.gz (.tbi required)")
    argparser.add_argument("--chrom", metavar="S", help="Limit analysis to the specified chromosome")
    argparser.add_argument("--outbam", metavar="out.bam", help="Input BAM with concordance annotations for relevant records")
    
    return argparser.parse_args(args)

#def main():
args = parseParam(sys.argv[1:], dummy=True)

# Open files
inbam = pysam.AlignmentFile(args.inbam)
reffasta = pysam.FastaFile(args.reffasta)

readcsv = open(args.readcsv, "w")

genomecsv = open(args.genomecsv, "w")
svcsv = open(args.svcsv, "w")
genomecsv.write("#chrom,pos,alignmentCount,diffCount,substitutionCount,deletionCount,insertionCount,singletonCount\n")
svcsv.write("#chrom,pos,alignmentCount,diffCount,substitutionCount,deletionCount,insertionCount,singletonCount\n")

hcregionsbed = pysam.TabixFile(args.hcregions) if args.hcregions else None
hcvariantsvcf = pysam.TabixFile(args.hcvariants) if args.hcvariants else None
outbam = pysam.AlignmentFile(args.outbam, "wb", template=inbam) if args.outbam else None

activeChrom,activeChromLength,activeHcRegions, = None, None, None
activeHcVariants = collections.defaultdict(list)
activeAlns, rpos, regionSize, svRange = [], 1, 50000, 100
count = 0

# Process alignments in sort order.
for aln in (inbam if not args.chrom else inbam.fetch(args.chrom)):
    
    count += 1
    if count % 1000 == 0: print(count)
    
    # For a new chromosome, update the list of HC regions and HC variants.
    if "chr" + aln.reference_name != activeChrom:
        activeChrom = "chr" + aln.reference_name
        activeChromLength = reffasta.get_reference_length(activeChrom)
        print(activeChrom)

        if hcregionsbed:
            # List of non-overlapping, high-confidence regions, sorted by chromStart.
            activeHcRegions = collections.deque(Bed(x) for x in hcregionsbed.fetch(aln.reference_name))
        else: # treat the whole chromosome as a high-confidence region
            activeHcRegions = collections.deque([Bed("%s\t0\t%d" % (activeChrom, activeChromLength))])

        activeHcVariants = collections.defaultdict(list)
        if hcvariantsvcf:
            # Map from chromStart to known variants.
            for line in hcvariantsvcf.fetch(aln.reference_name):
                for vcf in Vcf.LineToVcfs(line):
                    activeHcVariants[vcf.chromStart].append(vcf)

    # Pop high confidence regions that end before the current alignment.
    while len(activeHcRegions) and activeHcRegions[0].chromEnd <= aln.reference_start:
        activeHcRegions.popleft()

    # Only measure concordance for non-secondary alignments that overlap a high confidence region.
    if (not aln.is_secondary and
        len(activeHcRegions) and
        max(activeHcRegions[0].chromStart, aln.reference_start) < min(activeHcRegions[0].chromEnd, aln.reference_end)):

        subreadPasses = str(aln.get_tag("np")) if aln.has_tag("np") else "NA"
        predictedConcordance = "%0.6f" % (aln.get_tag("rq")) if aln.has_tag("rq") else "NA"

        alnStats = measure_concordance(aln, activeChrom, reffasta, activeHcRegions, activeHcVariants)
        activeAlns.append((aln, alnStats))

        if False and alnStats.hcReadLength > 0:
            readcsv.write("%s,%d,%s,%s,%s,%d,%d,%0.6f,%0.2f,%d,%d,%d,%d,%d\n" % (aln.query_name, aln.query_length, subreadPasses, predictedConcordance,
                "Supplementary" if aln.is_supplementary else "Primary", aln.mapping_quality, alnStats.hcReadLength,
                alnStats.concordance, alnStats.concordanceQv, alnStats.mismatchBp, alnStats.nonHpInsertionBp, alnStats.nonHpDeletionBp,
                alnStats.hpInsertionBp, alnStats.hpDeletionBp))


        if aln.reference_start - rpos > regionSize + svRange:
            measure_concordance_region(activeAlns, rpos, aln.reference_start-1, svRange, genomecsv, svcsv)          
            rpos = aln.reference_start-1
            
            while len(activeAlns) > 0:
                aln, _ = activeAlns[0]
                if aln.reference_end <= rpos - svRange: activeAlns.pop(0)
                else: break
            

            #statscsv.write("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%0.10f\n" % (aln.query_name, aln.query_length, alnStats.hcReadLength, alnStats.errors,
            #    alnStats.matches, alnStats.mismatches, alnStats.deletions, alnStats.hpDeletions, alnStats.insertions, alnStats.hpInsertions, alnStats.concordance))
            #aln.set_tag("rc", alnStats.concordance, "f")
            #aln.set_tag("tl", alnStats.hcReadLength, "i")

    

    measure_concordance_region(activeAlns, rpos, activeChromLength, svRange, genomecsv, svcsv)          


    if outbam:
        outbam.write(aln)


# Close input files
inbam.close()
reffasta.close()

# Close output file
if outbam:
    outbam.close()
readcsv.close()
genomecsv.close()


#if __name__ == "__main__":
#    main()
