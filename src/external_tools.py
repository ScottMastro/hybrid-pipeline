import pysam
import pyfaidx
import subprocess
import re
import os
from Bio import SeqIO
import glob
import shutil
import gzip

WHATSHAP = "whatshap"
SAMTOOLS = "samtools"
BCFTOOLS = "bcftools"
BEDTOOLS = "bedtools"
BGZIP = "bgzip"
TABIX = "tabix"
LONGSHOT = "longshot"
PBSV = "pbsv"
GENOMIC_CONSENSUS = "gcpp"
CANU = "/home/scott/bin/canu-1.9/Linux-amd64/bin/canu"
PBMM2 = "pbmm2"
BWA = "bwa"

MM2 = "/home/scott/bin/minimap2-2.17_x64-linux/minimap2"
NANOPLOT = "NanoPlot"
VG="/home/scott/bin/vg"
GRAPHVIZ_DOT="dot"
HAP_PY="/home/scott/bin/hap.py-install/bin/hap.py"
PRE_PY="/home/scott/bin/hap.py-install/bin/pre.py"

RTG = "/home/scott/bin/rtg-tools-3.10.1/rtg"
LONGRANGER="/home/scott/bin/longranger-2.2.2/longranger"
LONGRANGER_CORES=16
LONGRANGER_MEM=16
BAM2FQ_10X="/home/scott/bin/bamtofastq-1.2.0"

#required for longranger wgs
#version of GATK must be 3.3-3.8, or 4 except 3.6
GATK_10x_JAR="/home/scott/bin/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar"

PILON=["java","-Xmx16G","-jar","/home/scott/bin/pilon-1.23.jar"]

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

def create_seq_dict(fastaFa):
    
    outName = ".".join(fastaFa.split(".")[:-1]) + ".dict"    
    createDict = PICARD + ["CreateSequenceDictionary", "R=" + fastaFa, "O=" + outName]
    
    subprocess.call(createDict)
    return outName 


#samtools
#========================================================

def samtools_index(bamFile):    
    """Index a bam file."""
    subprocess.call([SAMTOOLS, "index", bamFile])
    return bamFile + ".bai"

def sam2bam(samFile):
    view = [SAMTOOLS, "view", "-S", "-b", samFile]
    base = ".".join(samFile.split(".")[:-1]) 
    unsorted = base + ".unsorted.bam"
    
    writer = open(unsorted, 'w+')    
    subprocess.call(view, stdout=writer)
    writer.close()

    outName = base + ".sorted.bam"
    sort = [SAMTOOLS, "sort", unsorted, "-o", outName]
    subprocess.call(sort)

    os.remove(unsorted)

    samtools_index(outName)
    return outName

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

def samtools_fetch(bamFile, region=None, unaligned=True):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    if region is None:
        alignments = [x for x in samfile.fetch(until_eof=unaligned)]
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
        
    alignments = sorted(alignments, key=lambda x: (x.is_unmapped, x.rname, x.pos))

    header = pysam.AlignmentFile(headerBam, "rb")
    bamFile = pysam.AlignmentFile(outName, "wb", template=header)
    for alignment in alignments:
         bamFile.write(alignment)

    bamFile.close()
    header.close()
    #samtools_index(outName)
    return outName

def faidx_fetch(faFile, region):
    fa = pyfaidx.Faidx(faFile)
    return str(fa.fetch(region.chrom, region.start, region.end))

#conversion
#========================================================
def bam2fastq_10x(bamFile, prefix, dirName=None, region=None, returnFqList=False):
    
    bam2fq = [BAM2FQ_10X]
            
    if region is None:
        outDir = prefix + "10x_reads"
    else:
        outDir = prefix + re.sub(":", "_", str(region) + "_10x")
        bam2fq.append("--locus=" + str(region))
    
    if dirName is not None:
        outDir = prefix + dirName

    bam2fq = bam2fq + [bamFile, outDir]

    subprocess.call(bam2fq)
    
    files = glob.glob(outDir + "/*/*")
    extraDir = glob.glob(outDir + "/*")[0]
    for file in files:
        shutil.move(file, extraDir + "/..")
    try:
        shutil.rmtree(extraDir)
    except: 
        pass
    
    if returnFqList:
        return sorted(glob.glob(outDir + "/*R[1,2]*.f*q.gz"))

    return outDir + "/"

def bam2fastq(bamFile):

    outName = ".".join(bamFile.split(".")[:-1]) + ".fq"    
    writer = open(outName, 'w+')    
    subprocess.call([SAMTOOLS, "bam2fq", bamFile], stdout=writer)
    writer.close()
    return outName

def reads2fasta(alignments, prefix):

    fastaLines = dict()
    fastaFile = prefix + (".fasta" if not prefix.endswith(".fasta") else "")

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
    
    outFile = prefix + (".fasta" if not prefix.endswith(".fasta") else "")
    writer = open(outFile, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")

    writer.close()
    if index: samtools_faidx(outFile)
    return outFile

def get_fasta_seq(faFile, region=None, toUpper=False):   
    faDict = fasta2dict(faFile, toUpper=toUpper)
    
    if region is None:
        fid = list(faDict.keys())[0]
        return faDict[fid]

    fid = region.chrom
    return faDict[fid][region.start:region.end]

def rename_fasta(faFile, name, toUpper=False):
    seq = get_fasta_seq(faFile, toUpper=toUpper)
    faDict = {name : seq}
    temp = faFile + "__temporaryfa__.fasta"
    
    dict2fasta(faDict, temp, index=False)
    shutil.move(temp, faFile)
    
    samtools_faidx(faFile)

    return faFile


#bcftools
#========================================================

def sort_vcf(vcfFile, removeOriginal=False, index=False):
    outName = re.sub(".vcf", ".sorted.vcf", vcfFile)
    sort = [BCFTOOLS, "sort"]
    if index:
        sort.append("-Ob")
        outName.append(".gz")
        
    sort = sort + ["-o", outName, vcfFile]

    print(sort)
    subprocess.call(sort)
    if removeOriginal:
        os.remove(vcfFile)
    return outName

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

def whatshap_phase(inputList, vcfFile, refFa, prefix, genotype=True, indels=True, maxCov=15):
    
    outName = prefix + ".whatshap.vcf"

    whatshapPhase = [WHATSHAP, "phase", "-o", outName, "--ignore-read-groups"]
    
    if genotype:
        whatshapPhase.append("--distrust-genotypes")
    if indels:
        whatshapPhase.append("--indels")

    whatshapPhase.extend(["--max-coverage", str(maxCov)])

    whatshapPhase.extend(["--reference", refFa, vcfFile])
    whatshapPhase.extend(inputList)

    subprocess.call(whatshapPhase)
    
    return bgzip(outName)

def whatshap_haplotag(bamFile, vcfFile, refFa):
    
    outName = re.sub(".bam", ".haplotag.bam", bamFile)

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

#pilon 
#========================================================

def pilon_polish(bamFile, refFa, outdir, prefix=None, snpIndelOnly=False,
                 outputVCF=False, changeFile=False, tracks=False):
    
    pilon = [x for x in PILON]
    
    pilon.extend(["--genome", refFa])
    pilon.extend(["--frags", bamFile])

    if outputVCF:
        pilon.append("--vcf")
        
    if changeFile:
        pilon.append("--changes")
    if tracks:
        pilon.append("--tracks")

    if snpIndelOnly:
        pilon.extend(["--fix", "snps,indels"])

    pilon.extend(["--outdir", outdir])

    subprocess.call(pilon)
    
    outName = outdir + "pilon.fasta"
    if prefix is not None:
        newOutName = outdir + prefix + ".fasta"
        shutil.move(outName, newOutName)
        outName = newOutName
        
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

    arrow.append(bamFile)
    
    print(arrow)
    subprocess.call(arrow)
    
    if outputFasta: 
        samtools_faidx(outFa)    
        return outFa
    
    return outVcf

#longranger
#========================================================

def mkref_10x(refFa):
    
    faName = refFa.split("/")[-1]
    faBase = ".".join(faName.split(".")[:-1]) 

    print(os.getcwd())
    indexDirName = "./refdata-" + faBase
    #remove existing
    try:
        shutil.rmtree(indexDirName)
    except:
        pass
    
    subprocess.call([LONGRANGER, "mkref", refFa])
    create_seq_dict(indexDirName + "/fasta/genome.fa")

    return indexDirName
    
def align_10x(refFa, fqDir, prefix, useFreebayes=False, correctBug=True):
    
    sampleName = ".".join(refFa.split(".")[:-1]).split("/")[-1] + "-longranger"
    cwd=os.getcwd()
    os.chdir(prefix)

    reference = mkref_10x(refFa)
    
    longrangerWGS = [LONGRANGER, "wgs", "--id=" + sampleName,
                     "--fastqs=" + fqDir, "--reference=" + reference,
                     "--jobmode=local", "--disable-ui", "--sex=female",
                     "--nopreflight", "--noloupe"]
    
    if useFreebayes:
        longrangerWGS.append("--vcmode=freebayes")
    else:
        longrangerWGS.append("--vcmode=gatk:" + GATK_10x_JAR)

    longrangerWGS.extend(["--localcores=" + str(LONGRANGER_CORES),
                          "--localmem=" + str(LONGRANGER_MEM)])

    subprocess.call(longrangerWGS)
    
    os.chdir(cwd)

    #clean up output directory
    files = glob.glob(prefix + sampleName + "/*")
    nfiles = len(files)
    files = [file for file in files if not file.endswith("outs")]
    if len(files) +1 == nfiles:
        for file in files:
            try:
                os.remove(file)
            except:
                try:
                    shutil.rmtree(file)
                except:
                    pass
    
    outs = glob.glob(prefix + sampleName + "/*/*")
    for out in outs:
        try:
            shutil.move(out, prefix + sampleName + "/")
        except:
            pass

    try:
        shutil.rmtree(prefix + sampleName + "/outs")
    except:
        pass
    

    return [prefix + sampleName + "/phased_variants.vcf.gz",
            prefix + sampleName + "/phased_possorted_bam.bam"]
    
#the 16-base 10x barcode plus 7 additional bases
def trim_10x_barcode(fastq, trim=23):
    reader = gzip.open(fastq, "rt" if fastq.endswith("gz") else "r")
    trimmedFq = "/".join(fastq.split("/")[:-1]) +  "/trimmed" + str(trim) + "_" + fastq.split("/")[-1] 
    writer = gzip.open(trimmedFq, 'wt')

    line = reader.readline()

    while(line):        
        if not line.startswith("@") and not line.startswith("+"):
            line = line[trim:]            
        writer.write(line)
        line = reader.readline()
        
    reader.close()
    writer.close()
    return trimmedFq

#bwa
#========================================================

def bwa_index(refFa):
    index = [BWA, "index", refFa]
    subprocess.call(index)
    

def bwa_mem(refFa, read1, read2, prefix):
    
    bwa = [BWA, "mem", refFa, read1, read2]
    bwa_index(refFa)
    
    samFile = prefix + ".sam"
    writer = open(samFile, 'w+')
    print(bwa)
    subprocess.call(bwa, stdout=writer)
    writer.close()

    outName = sam2bam(samFile)
    os.remove(samFile)
    return outName
    
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

def align_bam(refFa, readsFile, prefix, asm="asm5", r=500, secondary=False):
    
    outName = prefix + ".bam"    
    mm2Align = [MM2, "-acx", asm, "--eqx", "-r", str(r)]
    
    if not secondary:
        mm2Align.append("--secondary=no")
    
    mm2Align.extend(["-o", outName, refFa, readsFile])
    subprocess.call(mm2Align)
    
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
   
#hap.py
#========================================================
def hap_py(refFa, vcf1, vcf2, prefix, outPrefix="happy", engine="xcmp"):
    outDir=prefix + "happy/"
    try:
        os.mkdir(outDir)
    except: 
        pass
    
    happy = [HAP_PY, "-r", refFa, "-L", "--preprocess-truth", "--engine", engine]
    
    if engine == "vcfeval":
        happy.append("--engine-vcfeval-path")
        happy.append(RTG)

    happy = happy + ["-o", outDir + outPrefix]    
    happy.extend([vcf1, vcf2])
    
    print(happy)
    subprocess.call(happy)
    
    return outDir + outPrefix + ".vcf.gz"

def pre_py(refFa, vcfFile, prefix=None, decompose=True, noGz=False):
    
    if prefix is None:
        outname = re.sub(".gz", "", vcfFile)
        outName = re.sub(".vcf", ".normalized.vcf.gz", outname)
    else:
        outName = prefix + ".vcf.gz"

    prepy = [PRE_PY, "-r", refFa]
    if decompose:
        prepy.append("--decompose")
    prepy.extend([vcfFile, outName])
    subprocess.call(prepy)
    
    if noGz:
        newName = re.sub(".gz", "", outName)
        unzip = [BGZIP, "-d", outName]
        subprocess.call(unzip)
        return newName

    return outName

#RTG
#========================================================
def make_sdf(refFa, sdfName):    
    rtgformat = [RTG, "format", "-o", sdfName, refFa]
    subprocess.call(rtgformat)
    
    return sdfName

def vcfeval(refFa, baselineVCF, callsVCF, prefix, dirname="vcfeval_out", 
            squashPloidy=False, useQUAL=True, useFilter=False, outputType="split"):
    #rtgtools vcfeval -b consensus-longranger/phased_variants.vcf.gz -c AB_temp_merged.vcf -o rtg
    vcfeval = [RTG, "vcfeval", "-b", baselineVCF, "-c", callsVCF]

    if squashPloidy:
        vcfeval.append("--squash-ploidy")   
    vcfeval.append("--output-mode="+outputType)   

    vcfeval.append("--vcf-score-field=QUAL")   
    
    if not useFilter:
        vcfeval.append("--all-records")   


    sdf = prefix + "sdf"
    try:
        sdf = make_sdf(refFa, sdf)
    except:
        print("SDF already exists:" + prefix + "temp_sdf")
        
    outdir = prefix + dirname

    try:
        shutil.rmtree(outdir,ignore_errors=True)
    except:
        pass
    
    vcfeval = vcfeval + ["-t", sdf, "-o", outdir]
    subprocess.call(vcfeval)
    
    return outdir + "/"

#vg
#========================================================
def construct_graph(faFile, vcfFiles, prefix):

    if not isinstance(vcfFiles, list):
        vcfFiles = [vcfFiles]
        
    construct = [VG, "construct", "-a", "-r", faFile, "-v"]
    construct = construct + vcfFiles
    
    outName = prefix + ".vg"
    
    writer = open(outName, 'w+')    
    subprocess.call(construct, stdout=writer)
    writer.close()

    return outName

def construct_graph_msga(faFiles, prefix, normalize=False,
                         bedFile=None, renameSeqs=None, baseSeq=None):

    rename = [x for x in renameSeqs]
    if not isinstance(faFiles, list):
        combinedFa = faFiles
    else:
        faDict = dict()
        for fa in faFiles:
            d = fasta2dict(fa)
            if rename == None:
                faDict.update(d)
            else:
                fid = list(d.keys())[0]
                faDict[rename.pop(0)] = d[fid]
            
        combinedFa = dict2fasta(faDict, prefix + "_temp_msga_fasta_")

    msga = [VG, "msga", "-f", combinedFa]
    if baseSeq is not None:
        msga.extend(["-b", baseSeq])
    if bedFile is not None:
        msga.extend(["-R", bedFile])
    if normalize:
        msga.append("-N")

    outName = prefix + (".vg" if not prefix.endswith(".vg") else "")
    
    writer = open(outName, 'w+')    
    subprocess.call(msga, stdout=writer)
    
    if isinstance(faFiles, list):
        os.remove(combinedFa)
        os.remove(combinedFa + ".fai")

    writer.close()
    return outName
    
def index_graph(vgFile, xg=True, gcsa=True):

    prefix = re.sub(".vg", "", vgFile)
    index = [VG, "index"]
    
    if xg:
        index.extend(["-x", prefix + ".xg"])
    if gcsa:
        index.extend(["-g", prefix + ".gcsa"])

    index.append(vgFile)
    subprocess.call(index)

    return prefix


def index_graph_haplotypes(vgFile, phasedVCF, writeHaps=None):

    prefix = re.sub(".vg", "", vgFile)
    index = [VG, "index", "-T", "-v", phasedVCF, "-G",  prefix + ".gbwt"]
    if writeHaps is not None:
        index.extend(["-H", writeHaps])

    index.append(vgFile)
    subprocess.call(index)

    return prefix

def map_graph(inFile, vgFile, name=None):

    typeFlag = None
    
    prefix = re.sub(".vg", "", vgFile)
    
    if inFile.endswith(".fa") or inFile.endswith(".fasta"):
        typeFlag = "-F"
    if inFile.endswith(".bam") or inFile.endswith(".sam"):
        typeFlag = "-b"
    
    map = [VG, "map", "-x", prefix + ".xg",
             "-g", prefix + ".gcsa"]
    
    if name is not None:
        map = map + ["-V", str(name)]

    map = map + [typeFlag, inFile]

    outName = ".".join(inFile.split(".")[:-1]) + ".gam"

    writer = open(outName, 'w+')
    subprocess.call(map, stdout=writer)
    writer.close()

    return outName

def view_graph(vgFile, node, width=10):

    graphPrefix = re.sub(".vg", "", vgFile)

    find = [VG, "find", "-n", str(node), 
            "-x", graphPrefix + ".xg", "-c", str(width)]
    view = [VG, "view", "-dp", "-"]  

    pdfName = "_".join([graphPrefix, str(node), "w" + str(width)]) + ".pdf"
    dot = [GRAPHVIZ_DOT, "-Tpdf", "-o", pdfName]

    p1 = subprocess.Popen(find, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(view, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()  
    writer = open(pdfName, 'w+')    
    subprocess.Popen(dot, stdin=p2.stdout, stdout=writer)
    p2.stdout.close() 
    writer.close()

    return pdfName

def call_graph(vgFile, gamFile):

    graphPrefix = re.sub(".vg", "", vgFile)
    gamPrefix = re.sub(".gam", "", gamFile)

    packFile = gamPrefix + ".pack"
    pack = [VG, "pack", "-g",  gamFile,
            "-x", graphPrefix + ".xg", "-o", packFile]

    subprocess.call(pack)

    call = [VG, "call", graphPrefix + ".xg", "-k",  packFile]

    outName = gamPrefix + ".vcf"

    writer = open(outName, 'w+')
    subprocess.call(call, stdout=writer)
    writer.close()
    
    return outName

def augment_graph(vgFile, gamFile):
    
    augment = [VG, "augment", "-B",  vgFile, gamFile]

    tempOutName = vgFile + ".temp.vg"

    writer = open(tempOutName, 'w+')
    subprocess.call(augment, stdout=writer)
    writer.close()
    
    os.replace(tempOutName, vgFile)
    
    return vgFile

def add_path_graph(vgFile, faFile, name, indexed=False):
    
    if not indexed:
        index_graph(vgFile)
        
    gam = map_graph(faFile, vgFile, name=name)
    return augment_graph(vgFile, gam)

def graph_to_vcf(vgFile, refPaths, altPaths, prefix, index=False):
    if not isinstance(refPaths, list):
        refPaths = [refPaths]
    if not isinstance(altPaths, list):
        altPaths = [altPaths]
        
    deconstruct = [VG, "deconstruct", "-e"]
    deconstruct = deconstruct + ["-p", ",".join(refPaths)]
    deconstruct = deconstruct + ["-A", ",".join(altPaths)]
    deconstruct.append(vgFile)
        
    outName = prefix + (".vcf" if not prefix.endswith(".vcf") else "")

    writer = open(outName, 'w+')
    subprocess.call(deconstruct, stdout=writer)
    writer.close()
    
    outName = sort_vcf(outName, removeOriginal=False, index=index)        

    return outName