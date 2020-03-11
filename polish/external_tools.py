import pysam
import subprocess
import re
from Bio import SeqIO

WHATSHAP = "whatshap"
SAMTOOLS = "samtools"
BGZIP = "bgzip"
TABIX = "tabix"
LONGSHOT = "longshot"
PBSV = "pbsv"
GENOMIC_CONSENSUS = "gcpp"
CANU = "/home/scott/bin/canu-1.9/Linux-amd64/bin/canu"
PBMM2 = "pbmm2"
MM2 = "/home/scott/bin/minimap2-2.17_x64-linux/minimap2"
NANOPLOT = "NanoPlot"

PICARD=["java","-Xmx8G","-jar","/home/scott/bin/picard.jar"]

#picard
#========================================================

def unalign_bam(bamFile, prefix):    
    
    revert = PICARD + ["RevertSam"]
    revert.append("I=" + bamFile)
    unaligned = prefix + ".unaligned.bam"
    revert.append("O=" + unaligned)
    
    subprocess.call(revert)
    return unaligned 

#samtools
#========================================================

def samtools_index(bamFile):    
    """Index a bam file."""
    subprocess.call([SAMTOOLS, "index", bamFile])
    return bamFile + ".bai"

def samtools_faidx(faFile):    
    """Index a bam file."""
    subprocess.call([SAMTOOLS, "faidx", faFile])
    return faFile + ".fai"

def samtools_subset(bamFile, region, prefix):
    '''View and index a subset of a bam file'''
    
    outName = prefix + ".subset.bam"

    writer = open(outName, 'w+')    
    subprocess.call([SAMTOOLS, "view", "-b", bamFile, str(region)], stdout=writer)
    writer.close()
        
    samtools_index(outName)
    return outName

def samtools_fetch(bamFile, region=None):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    if region is None:
        alignments = [x for x in samfile.fetch()]
    else:
        alignments = [x for x in samfile.fetch(region.chrom, region.start, region.end)]
    return alignments

def samtools_write(alignments, prefix, headerBam, makeUnique=False):
    outName = prefix + ".bam"

    # remove redundant reads:
    if makeUnique:
        readDict = {read.qname : read for read in alignments}
        uniqueReads = [readDict[k] for k in readDict]
        alignments = uniqueReads

    alignments = sorted(alignments, key=lambda x: (x.rname, x.pos))

    header = pysam.AlignmentFile(headerBam, "rb")
    bamFile = pysam.AlignmentFile(outName, "wb", template=header)
    for alignment in alignments:
         bamFile.write(alignment)

    bamFile.close()
    header.close()
    samtools_index(outName)
    return outName

#conversion
#========================================================

def bam2fastq(bamFile):

    outName = ".".join(bamFile.split(".")[:-1]) + ".fasta"
    writer = open(outName, 'w+')    
    subprocess.call([SAMTOOLS, "bam2fq", bamFile], stdout=writer)
    writer.close()
    return outName

def reads2fasta(alignments, prefix):

    fastaLines = dict()
    fastaFile = prefix + ".fasta"
    
    for alignment in alignments:
        fastaLines[alignment.qname] = alignment.seq
    
    writer = open(fastaFile, 'w+')
    for readId in fastaLines:
        writer.write(">" + readId + "\n" + fastaLines[readId] + "\n")
    writer.close()

    return fastaFile

def fasta2dict(faFile, toUpper=False):
    fastaDict = dict()
    seqs = SeqIO.parse(open(faFile),'fasta')
    for f in seqs:
        fastaDict[f.id] = str(f.seq)
        if toUpper: fastaDict[f.id] = fastaDict[f.id].upper()
    return fastaDict

def dict2fasta(fastaDict, prefix, toUpper=False, index=True):
    
    outFile = prefix + ".fasta"
    writer = open(outFile, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")

    writer.close()
    if index: samtools_faidx(outFile)
    return outFile

#bgzip
#========================================================

def bgzip(vcfFile):
    '''Call bgzip to compress a VCF file.'''
    subprocess.call([BGZIP, "-f", vcfFile])
    subprocess.call([TABIX, vcfFile + ".gz"])
    return(vcfFile + ".gz")

#longshot
#========================================================

LONGSHOT_BAM_PREFIX = ".longshot"

def longshot_genotype(bamFile, refFa, prefix, region=None, coverageAware=True, writeBams=False):
    
    outName = prefix + ".longshot.vcf"
    
    longshot = [LONGSHOT]
    if coverageAware:
        longshot.append("-A")
    
    if writeBams:
        longshot.extend(["--hap_bam_prefix", prefix + LONGSHOT_BAM_PREFIX])

    #force overwrite of output file
    longshot.append("-F")

    if region is not None:
        longshot.extend(["--region", region.str_base1()])

    longshot.extend(["--bam", bamFile, "--ref", refFa, "--out", outName])
    subprocess.call(longshot)
    
    if writeBams:
        bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
        hpA= bamFilePrefix + ".hap1.bam"
        hpB = bamFilePrefix + ".hap2.bam"
        unphased = bamFilePrefix + ".unassigned.bam"
        samtools_index(hpA)
        samtools_index(hpB)
        samtools_index(unphased)
        
    return outName

def get_longshot_phased_reads(prefix):
    '''
    --hap_bam_prefix
    Write haplotype-separated reads to 3 bam files using this prefix:
   <prefix>.hap1.bam, <prefix>.hap2.bam, <prefix>.unassigned.bam
   '''
   
    bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
    h1= bamFilePrefix + ".hap1.bam"
    h2 = bamFilePrefix + ".hap2.bam"
    unphased = bamFilePrefix + ".unassigned.bam"
        
    return [h1, h2, unphased]

#whatshap
#========================================================

def whatshap_phase(inputList, vcfFile, refFa, prefix, genotype=True, indels=True):
    
    outName = prefix + ".whatshap.phase.vcf"

    whatshapPhase = [WHATSHAP, "phase", "-o", outName, "--ignore-read-groups"]
    
    if genotype:
        whatshapPhase.append("--distrust-genotypes")
    if indels:
        whatshapPhase.append("--indels")

    whatshapPhase.extend(["--reference", refFa, vcfFile])
    whatshapPhase.extend(inputList)

    subprocess.call(whatshapPhase)
    return outName

def whatshap_haplotag(bamFile, vcfFile, refFa, prefix):
    
    outName = prefix + ".whatshap.haplotag.vcf"

    whatshapTag = [WHATSHAP, "haplotag", "-o", outName, "--ignore-read-groups"]
    whatshapTag.extend(["--reference", refFa, vcfFile, bamFile])

    subprocess.call(whatshapTag)
    samtools_index(outName)
    return outName


def whatshap_genotype(bamFile, vcfFile, refFa, prefix, indels=True):
    
    outName = prefix + ".whatshap.genotype.vcf"

    whatshapGenotype = [WHATSHAP, "genotype", "-o", outName, "--ignore-read-groups"]
    
    whatshapGenotype.extend(["--reference", refFa])
    if indels:
        whatshapGenotype.append("--indels")

    whatshapGenotype.extend([vcfFile, bamFile])
    subprocess.call(whatshapGenotype)
    return outName

'''
def whatshap_find_snv(outname, reference, vcf, bamfile):
    
    whatshapSNV = [WHATSHAP, "find_snv_candidates", reference, bamfile, "-o", outname]
    subprocess.call(whatshapSNV)
    return outname
'''

#pbsv
#========================================================

def find_sv(bamFile, refFa, prefix, threads=6, indelOnly=False):

    outName = prefix + ".pbsv.vcf"
    signatureFile = prefix + "pbsv.svsig.gz"

    discover = ["pbsv", "discover", "--sample", "hybrid", bamFile, signatureFile]
    call = ["pbsv", "call", "-j", str(threads)]

    if indelOnly:
        call.extend(["--types", "DEL,INS"])
        
    call.extend([refFa, signatureFile, outName])

    subprocess.call(discover)
    subprocess.call(call)
    return outName

#genomicconsensus (arrow)
#========================================================

def pacbio_polish(bamFile, refFa, prefix, region=None, outputFasta=False, outputGff=False, useMapFilter=True, chunkSize=None):
    #samtools view -H $BAM | sed "s/VN:1.3/VN:1.3\tpb:3.0.4/" | samtools reheader - $BAM
        
    '''What is MapQV and why is it important?

    MapQV is a single scalar Phred-scaled QV per aligned read that reflects the mapper's degree of certainty that the read aligned to this part of the reference and not some other. Unambigously mapped reads will have a high MapQV (typically 255), while a read that was equally likely to have come from two parts of the reference would have a MapQV of 3.
    MapQV is pretty important when you want highly accurate variant calls. Quiver and Plurality both filter out aligned reads with a MapQV below 20 (by default), so as not to call a variant using data of uncertain genomic origin.
    This can be problematic if using quiver/arrow to get a consensus sequence. If the genome of interest contains long (relative to the library insert size) highly-similar repeats, the effective coverage (after MapQV filtering) may be reduced in the repeat regions---this is termed these MapQV dropouts. If the coverage is sufficiently reduced in these regions, quiver/arrow will not call consensus in these regions---see `What do quiver/arrow do for genomic regions with no effective coverage?`_.
    If you want to use ambiguously mapped reads in computing a consensus for a denovo assembly, the MapQV filter can be turned off entirely. In this case, the consensus for each instance of a genomic repeat will be calculated using reads that may actually be from other instances of the repeat, so the exact trustworthiness of the consensus in that region may be suspect. The next section describes how to disable the MapQV filter.
    How can the MapQV filter be turned off and when should it be?
    The MapQV filter can be disabled using the flag --mapQvThreshold=0 (shorthand: -m=0). If running a quiver/arrow job via SMRT Portal, this can be done by unchecking the "Use only unambiguously mapped reads" option. Consider this in de novo assembly projects, but it is not recommended for variant calling applications.
    '''    
    
    outVcf = prefix + ".consensus.vcf"
    outFa =  prefix + ".consensus.fasta"
    outGff = prefix + ".consensus.gff"

    outputFiles = outVcf

    if outputFasta:
        outputFiles = outputFiles + "," + outFa
    if outputGff:
        outputFiles = outputFiles + "," + outGff
    arrow = [GENOMIC_CONSENSUS, "-r", refFa, "-o", outputFiles ]
    arrow.extend(["--algorithm", "arrow"])
    if region is not None:
        arrow.extend(["--windows", str(region)])
        
    if not useMapFilter:
        arrow.append("-m=0")

    if chunkSize is not None:
        arrow.extend(["-C", str(chunkSize)])

    if outputGff:
        arrow.append("--annotate-gff")

    arrow.append(bamFile)
    
    print(arrow)
    subprocess.call(arrow)
    
    if outputFasta: 
        samtools_faidx(outFa)    
        return outFa
    
    return outVcf

#mm2
#========================================================

def align_pacbio(refFa, readsFile, prefix, minConcordance=None, medianFilter=False):

    outName = prefix + ".pbmm2.bam"
    
    pbmm2Align = [PBMM2, "align", "--sort"]
    
    if minConcordance is not None:
       pbmm2Align.extend(["--min-concordance-perc", str(minConcordance)])
       
    if medianFilter:
       pbmm2Align.append("--median-filter")

    pbmm2Align.append("--unmapped")

    pbmm2Align.extend([refFa, readsFile, outName])
    subprocess.call(pbmm2Align)
    samtools_index(outName)
    
    return outName

def align_paf(refFa, readsFile, prefix, asm="asm5", r=500, secondary=False):
    '''
     -r
     Bandwidth used in chaining and DP-based alignment [500].
     This option approximately controls the maximum gap size.
     '''
    outName = prefix + ".paf"
    
    mm2Align = [MM2, "-cx", asm, "--cs", "-r", str(r)]
    
    if not secondary:
        mm2Align.append("--secondary=no")
    
    mm2Align.extend(["-o", outName, refFa, readsFile])
    subprocess.call(mm2Align)
    
    return outName

def align_paf_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=10000, asm="asm10", secondary=False)
def align_paf_very_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=20000, asm="asm20", secondary=False)


#canu
#========================================================

def canu_correct(fastaFile, prefix, size, trim=True, rawErrorRate=None):
    name = prefix.split("/")[-1]
    
    canuCorrect = [CANU, "-correct", "-p", name, "-d", prefix]
    canuCorrect.append("genomeSize=" + str(size/1000) + "k")
    if rawErrorRate is not None:
        canuCorrect.append("-rawErrorRate=" + str(rawErrorRate))

    canuCorrect.append("-pacbio-raw")
    canuCorrect.append(fastaFile)

    subprocess.call(canuCorrect)
    
    if trim:
        canuTrim = canuCorrect
        canuTrim[1] = "-trim"
        subprocess.call(canuTrim)

def canu_assemble(fastaFile, prefix, size):
    name = prefix.split("/")[-1]

    canuAssemble = [CANU, "-assemble", "-p", name, "-d", prefix]
    canuAssemble.append("genomeSize=" + str(size/1000) + "k")
    canuAssemble.append("-pacbio-corrected")
    canuAssemble.append(fastaFile)

    subprocess.call(canuAssemble)
    







