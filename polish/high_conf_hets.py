import sys
import vcf
import os 
import gzip

longshotVCF=sys.argv[1]
longrangerVCF=sys.argv[2]
outVCF=sys.argv[3]


def is_empty(fname):
    ''' Test if gzip file fname is empty
        Return True if the uncompressed data in fname has zero length
        or if fname itself has zero length
        Raises OSError if fname has non-zero length and is not a gzip file
    '''
    if fname.endswith(".gz"):
        with gzip.open(fname, 'rb') as f:
            data = f.read(1)
            return len(data) == 0
    return os.path.getsize(fname) == 0


if is_empty(longshotVCF):
    longshotVariants = []
else:
    longshotVariants = [record for record in vcf.Reader(open(longshotVCF, "rb" if longshotVCF.endswith(".gz") else "r"))]

if is_empty(longrangerVCF):
    longrangerVariants = []    
else:
    longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb" if longrangerVCF.endswith(".gz") else "r"))]

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

writer = open(outVCF, "w+")
writer.write("##fileformat=VCFv4.2\n")
header = ["#CHROM", "POS", "ID", "REF", "ALT",
       "QUAL",  "FILTER", "INFO", "FORMAT", "SAMPLE"]
writer.write("\t".join(header) + "\n")

for variant in sorted(highConfidenceVariants, key=lambda x: x.POS):
    line = [variant.CHROM, variant.POS, ".", variant.REF, variant.ALT[0], 
        "100", "PASS", ".", "GT", "0/1"]       
    writer.write("\t".join([str(x) for x in line]) + "\n")

writer.close()
