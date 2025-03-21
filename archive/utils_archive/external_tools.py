import pysam, pyfaidx
import subprocess
import re, glob, gzip
import os, shutil, copy
import sys
sys.path.append("../..")

import utils.fasta_handler as fasta

import utils.environment as env
#import utils.environment_hpf as env

def parse(command):
    return command.split()

#==================================================
# Bam helpers
#==================================================

def samtools_index(bamFile):    
    cmd = parse(env.SAMTOOLS)
    subprocess.call(cmd + ["index", bamFile])
    return bamFile + ".bai"

def sam_to_bam(samFile):
    cmd = parse(env.SAMTOOLS)
    view = cmd + ["view", "-S", "-b", samFile]
    base = ".".join(samFile.split(".")[:-1]) 
    unsorted = base + ".unsorted.bam"
    
    writer = open(unsorted, 'w+')    
    subprocess.call(view, stdout=writer)
    writer.close()

    outName = base + ".sorted.bam"
    sort = cmd + ["sort", unsorted, "-o", outName]
    subprocess.call(sort)

    os.remove(unsorted)

    samtools_index(outName)
    return outName

def samtools_fetch(bamFile, region=None, unaligned=True, asDict=False):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    if region is None:
        alignments = [x for x in samfile.fetch(until_eof=unaligned)]
    else:
        alignments = [x for x in samfile.fetch(region.chrom, region.start, region.end)]
    
    if asDict:
        alignDict = dict()
        for a in alignments:
            if a.qname not in alignDict: alignDict[a.qname] = []
            alignDict[a.qname].append(a)
        return alignDict

    return alignments

def samtools_write(alignments, prefix, headerBam, makeUnique=False, index=True):
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
    
    if index:
        samtools_index(outName)
    return outName

def bgzip(vcfFile):
    cmd1 = parse(env.BGZIP)
    cmd2 = parse(env.TABIX)

    '''Call bgzip to compress a VCF file.'''
    subprocess.call(cmd1 + ["-f", vcfFile])
    subprocess.call(cmd2 + [vcfFile + ".gz"])
    return(vcfFile + ".gz")

#==================================================
# Long Ranger
#==================================================

def bam2fastq_10x(bamFile, prefix, dirName=None, region=None, returnFqList=False):
    
    cmd = parse(env.BAM2FQ_10X)
    bam2fq = cmd
    
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

def create_seq_dict(fastaFa):
    
    cmd = parse(env.PICARD)

    outName = ".".join(fastaFa.split(".")[:-1]) + ".dict"    
    createDict = cmd + ["CreateSequenceDictionary", "R=" + fastaFa, "O=" + outName]
    
    subprocess.call(createDict)
    return outName 

def mkref_10x(refFa, svBlacklist=False):
    
    cmd = parse(env.LONGRANGER)

    faName = refFa.split("/")[-1]
    faBase = ".".join(faName.split(".")[:-1]) 

    print(os.getcwd())
    indexDirName = "./refdata-" + faBase
    
    #remove existing
    try:
        shutil.rmtree(indexDirName)
    except:
        pass
            
    subprocess.call(cmd + ["mkref", refFa])
    create_seq_dict(indexDirName + "/fasta/genome.fa")

    if svBlacklist:
        faDict = fasta.read_fasta(refFa)
        blacklist = indexDirName + "/regions/sv_blacklist.bed"
        writer = open(blacklist, 'w')

        for tig in faDict:
            line="\t".join([tig, "0", str(len(faDict[tig])-1), "contig"])
            writer.write(line)
        writer.close()

    return indexDirName
    
def align_10x(refFa, fqDir, prefix, useFreebayes=False, svBlacklist=False):
    
    cmd = parse(env.LONGRANGER)
    gatkJar = str(env.GATK_10x_JAR)
    cores = str(env.CORES)
    mem = str(env.MEM)

    sampleName = ".".join(refFa.split(".")[:-1]).split("/")[-1] + "-longranger"
    cwd=os.getcwd()
    os.chdir(prefix)

    reference = mkref_10x(refFa, svBlacklist)
    
    longrangerWGS = cmd + ["wgs", "--id=" + sampleName,
                     "--fastqs=" + fqDir, "--reference=" + reference,
                     "--jobmode=local", "--disable-ui", "--sex=female",
                     "--nopreflight", "--noloupe"]
    
    if useFreebayes:
        longrangerWGS.append("--vcmode=freebayes")
    else:
        longrangerWGS.append("--vcmode=gatk:" + gatkJar)

    longrangerWGS.extend(["--localcores=" + cores, "--localmem=" + mem])
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

#==================================================
# Longshot
#==================================================

LONGSHOT_BAM_SUFFIX = ".longshot.bam"

def longshot_genotype(bamFile, refFa, prefix, region=None, coverageAware=True, writeBams=False):
    
    outName = prefix + ".longshot.vcf"
    cmd = parse(env.LONGSHOT)
    longshot = cmd
    
    if coverageAware:
        longshot.append("-A")
    
    if writeBams:
        #old version of longshot
        #longshot.extend(["--hap_bam_prefix", prefix + LONGSHOT_BAM_PREFIX])
        longshot.extend(["--out_bam", prefix + LONGSHOT_BAM_SUFFIX])

    #force overwrite of output file
    longshot.append("-F")

    if region is not None:
        longshot.extend(["--region", region.str_base1()])

    longshot.extend(["--bam", bamFile, "--ref", refFa, "--out", outName])
    subprocess.call(longshot)
    
    if writeBams:
        '''
        bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
        hpA = bamFilePrefix + ".hap1.bam"
        hpB = bamFilePrefix + ".hap2.bam"
        unphased = bamFilePrefix + ".unassigned.bam"
        samtools_index(hpA)
        samtools_index(hpB)
        samtools_index(unphased)
        '''
        samtools_index(prefix + LONGSHOT_BAM_SUFFIX)

    return outName

def get_longshot_phased_reads(prefix, keepOriginal=False):
   
    #old version of longshot
    #--hap_bam_prefix
    #Write haplotype-separated reads to 3 bam files using this prefix:
    #<prefix>.hap1.bam, <prefix>.hap2.bam, <prefix>.unassigned.bam

    '''
    bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
    h1= bamFilePrefix + ".hap1.bam"
    h2 = bamFilePrefix + ".hap2.bam"
    unphased = bamFilePrefix + ".unassigned.bam"
    '''

    h1 = prefix + ".longshot.hap1"
    h2 = prefix + ".longshot.hap2"
    unphased = prefix + ".longshot.unassigned"

    if os.path.isfile(h1 +".bam" ) and os.path.isfile(h2 +".bam" ) and \
        os.path.isfile(unphased +".bam" ):
            return [x+".bam" for x in (h1, h2, unphased)]

    bamFile = prefix + LONGSHOT_BAM_SUFFIX 
    hap = {x:[] for x in ["0","1","2"]}

    
    for read in samtools_fetch(bamFile):
        try:
            hp = str(read.get_tag("HP"))
        except:
            hp = "0"
        
        hap[hp].append(read)
        

    samtools_write(hap["1"], h1, bamFile)
    samtools_write(hap["2"], h2, bamFile)
    samtools_write(hap["0"], unphased, bamFile)
    
    if not keepOriginal:
        try:
            os.remove(bamFile)
            os.remove(bamFile + ".bai")
        except:
            pass
        
    return [x+".bam" for x in (h1, h2, unphased)]

#==================================================
# Whatshap
#==================================================

def whatshap_phase(inputList, vcfFile, refFa, prefix, genotype=True, indels=True, maxCov=15):

    cmd = parse(env.WHATSHAP)
    outName = prefix + ".whatshap.vcf"

    whatshapPhase = cmd + ["phase", "-o", outName, "--ignore-read-groups"]
    
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
    
    cmd = parse(env.WHATSHAP)
    outName = re.sub(".bam", ".haplotag.bam", bamFile)

    whatshapTag = cmd + ["haplotag", "-o", outName, "--ignore-read-groups"]
    whatshapTag.extend(["--reference", refFa, vcfFile, bamFile])

    subprocess.call(whatshapTag)
    samtools_index(outName)
    return outName

#==================================================
# Polishing
#==================================================

def pilon_polish(bamFile, refFa, outdir, prefix=None, snpIndelOnly=False,
                 outputVCF=False, changeFile=False, tracks=False):
    
    cmd = parse(env.PILON)
    pilon = cmd + ["--genome", refFa, "--frags", bamFile]
    
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
    cmd = parse(env.GENOMIC_CONSENSUS)
    
    outVcf = prefix + ".consensus.vcf"
    outFa =  prefix + ".consensus.fasta"
    outGff = prefix + ".consensus.gff"

    outputFiles = outVcf

    if outputFasta:
        outputFiles = outputFiles + "," + outFa
    if outputGff:
        outputFiles = outputFiles + "," + outGff
    arrow = cmd + ["-r", refFa, "-o", outputFiles ]
    arrow.extend(["--algorithm", "arrow"])
    if region is not None:
        arrow.extend(["--windows", str(region)])
        
    if not useMapFilter:
        arrow.append("-m=0")

    if chunkSize is not None:
        arrow.extend(["-C", str(chunkSize)])

    arrow.append(bamFile)

    subprocess.call(arrow)
    
    if outputFasta: 
        pyfaidx.Faidx(outFa)    
        return outFa
    
    return outVcf

#==================================================
# BWA
#==================================================

def bwa_index(refFa):
    cmd = parse(env.BWA)
    subprocess.call(cmd + ["index", refFa])

def bwa_mem(refFa, read1, read2, prefix):
    
    cmd = parse(env.BWA)

    bwa = cmd + ["mem", refFa, read1, read2]
    bwa_index(refFa)
    
    samFile = prefix + ".sam"
    writer = open(samFile, 'w+')
    print(bwa)
    subprocess.call(bwa, stdout=writer)
    writer.close()

    outName = sam_to_bam(samFile)
    os.remove(samFile)
    return outName
    
#==================================================
# minimap2
#==================================================

def align_pacbio(refFa, readsFile, prefix, minConcordance=None, medianFilter=False):

    cmd = parse(env.PBMM2)
    outName = prefix + ".pbmm2.bam"
    
    pbmm2Align = cmd + ["align", "--sort"]
    
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
    
    cmd = parse(env.MM2)

    outName = prefix + ".bam"    
    mm2Align = cmd + ["-acx", asm, "--eqx", "-r", str(r)]
    
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
    cmd = parse(env.MM2)
    outName = prefix + ".paf"
    
    mm2Align = cmd + ["-cx", asm, "--cs", "-r", str(r)]
    
    if not secondary:
        mm2Align.append("--secondary=no")
    
    mm2Align.extend(["-o", outName, refFa, readsFile])
    subprocess.call(mm2Align)
    
    return outName

def align_paf_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=10000, asm="asm10", secondary=False)

def align_paf_very_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=20000, asm="asm20", secondary=False)

def paf_liftover(pafFile, regions, minMapQuality=5, maxDivergance=1, minAlnLen=2500):
    '''
    -q INT    min mapping quality [5]
    -l INT    min alignment length [50000]
    -d FLOAT  max sequence divergence (>=1 to disable) [1]
    '''
    cmd = parse(env.PAFTOOLS)
    bedFile = os.path.dirname(os.path.realpath(pafFile)) + "/_tempbed_.bed"
    
    if not isinstance(regions, list): regions = [regions]

    f = open(bedFile, "w")
    for region in regions:
        f.write("\t".join([str(x) for x in (region.chrom, region.start, region.end)]) + "\n")
    f.close()
    
    paftoolsLiftover = cmd + ["liftover"]  + \
                    ["-d", str(maxDivergance)] + \
                    ["-l", str(minAlnLen)] + \
                    ["-q", str(minMapQuality)]

    paftoolsLiftover.extend([pafFile])
    paftoolsLiftover.extend([bedFile])
    
    outFile = os.path.dirname(os.path.realpath(pafFile)) + "/_result_.bed"

    writer = open(outFile, 'w+')
    subprocess.call(paftoolsLiftover, stdout=writer)
    writer.close()

    liftover = []    
    for line in open(outFile, "r"):
        r = copy.deepcopy(regions[0])
        s = line.split("\t")
        r.chrom = s[0] ; r.start = int(s[1]); r.end = int(s[2])
        liftover.append(r)


    #todo: outfile missing matches?
    os.remove(bedFile)
    os.remove(outFile)

    return liftover

def paf_call(pafFile, referenceFa, prefix, vcfSample="sample", minAlnLen=2500):
    '''
      -f FILE   reference sequences (enabling VCF output) [null]
      -s NAME   sample name in VCF header [sample]
    '''
    cmd = parse(env.PAFTOOLS)

    paftoolsCall = cmd + ["call"]  + \
                    ["-f", str(referenceFa)] + \
                    ["-L", str(minAlnLen)] + \
                    ["-s", str(vcfSample)]

    paftoolsCall.extend([pafFile])
    
    outFile = prefix + ".vcf"

    print(" ".join(paftoolsCall) )
    writer = open(outFile, 'w+')
    subprocess.call(paftoolsCall, stdout=writer)
    writer.close()

    return outFile

#==================================================
# Graph
#==================================================
    
def seqwish_graph(faFile, pafFile, prefix, minMatchLen=None, repeatMax=None):
    # $SEQWISH -k 16 -s $FASTA -p $PAF -g $prefix.gfa

    outName = prefix + ".gfa"
    seqwish = parse(env.SEQWISH)
    
    if minMatchLen is not None:
        seqwish = seqwish + ["-k", str(minMatchLen)]
    if repeatMax is not None:
        seqwish = seqwish + ["-r", str(repeatMax)]

    seqwish + ["-r", str(repeatMax)]

    seqwish = seqwish + ["-s", faFile, "-p", pafFile, "-g", outName]
    subprocess.call(seqwish)

    return outName

def vg_normalize_gfa(gfaFile, prefix, generateGfa=False):

    #vg view -Fv $gfaFile | vg mod -X 256 - | vg mod -n - | vg ids -c -| vg mod -X 256 - > $prefix.vg
    
    outName = prefix + "_normalized.gfa"
    cmd = parse(env.VG)
    
    view = cmd + ["view", "-Fv", prefix]
    mod_size = cmd + ["mod", "-X", "256", "-"]
    mod_normalize = cmd + ["mod", "-n", "-"]
    ids = cmd + ["ids", "-c", "-"]

    p1 = subprocess.Popen(view, stdout=subprocess.PIPE)
    
    p2 = subprocess.Popen(mod_size, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    
    p3 = subprocess.Popen(mod_normalize, stdin=p2.stdout, stdout=subprocess.PIPE)
    p2.stdout.close()  
    
    p4 = subprocess.Popen(ids, stdin=p3.stdout, stdout=subprocess.PIPE)
    p3.stdout.close()  

    writer = open(outName, 'w+')    
    subprocess.Popen(mod_size, stdin=p4.stdout, stdout=writer)
    p4.stdout.close()
    writer.close()
    
    if generateGfa:
        #vg view -Vg $prefix.vg > $prefix_normalized.gfa

        gfaOut = prefix + "_normalized.gfa"
        gfaview = cmd + ["view", "-Vg", outName]

        writer2 = open(gfaOut, 'w+')    
        subprocess.Popen(gfaview, stdin=p4.stdout, stdout=writer2)
        writer.close()

    return outName

def vg_index(vgFile, xg=True, gcsa=False):

    cmd = parse(env.VG)

    prefix = re.sub(".vg", "", vgFile)
    index = cmd + ["index"]
    
    if xg:
        index.extend(["-x", prefix + ".xg"])
    if gcsa:
        index.extend(["-g", prefix + ".gcsa"])

    index.append(vgFile)
    subprocess.call(index)

    return prefix

def vg_paths_to_gam(xgFile, prefix):
    #vg paths -X -x prefix.xg > prefix.gam

    cmd = parse(env.VG)
       
    paths = cmd + ["paths", "-X", "-x", xgFile]
    
    outName = prefix + ".gam"
    
    writer = open(outName, 'w+')    
    subprocess.call(paths, stdout=writer)
    writer.close()

    return outName

#==================================================
# Legacy graph code
#==================================================

def construct_graph(faFile, vcfFiles, prefix):

    cmd = parse(env.VG)

    if not isinstance(vcfFiles, list):
        vcfFiles = [vcfFiles]
        
    construct = cmd + ["construct", "-a", "-r", faFile, "-v"]
    construct = construct + vcfFiles
    
    outName = prefix + ".vg"
    
    writer = open(outName, 'w+')    
    subprocess.call(construct, stdout=writer)
    writer.close()

    return outName

def construct_graph_msga(faFiles, prefix, normalize=False,
                         bedFile=None, renameSeqs=None, baseSeq=None):
    cmd = parse(env.VG)

    rename = [x for x in renameSeqs]
    if not isinstance(faFiles, list):
        combinedFa = faFiles
    else:
        faDict = dict()
        for fa in faFiles:
            d = fasta.read_fasta(fa)
            if rename == None:
                faDict.update(d)
            else:
                fid = list(d.keys())[0]
                faDict[rename.pop(0)] = d[fid]
            
        combinedFa = fasta.write_fasta(prefix + "_temp_msga_", faDict, index=True)

    msga = cmd + ["msga", "-f", combinedFa]
    msga.extend(["-t", str(env.CORES)])

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


def index_graph_haplotypes(vgFile, phasedVCF, writeHaps=None):

    cmd = parse(env.VG)

    prefix = re.sub(".vg", "", vgFile)
    index = cmd + ["index", "-T", "-v", phasedVCF, "-G",  prefix + ".gbwt"]
    if writeHaps is not None:
        index.extend(["-H", writeHaps])

    index.append(vgFile)
    subprocess.call(index)

    return prefix

def map_graph(inFile, vgFile, name=None):

    cmd = parse(env.VG)
    typeFlag = None
    
    prefix = re.sub(".vg", "", vgFile)
    
    if inFile.endswith(".fa") or inFile.endswith(".fasta"):
        typeFlag = "-F"
    if inFile.endswith(".bam") or inFile.endswith(".sam"):
        typeFlag = "-b"
    
    map = cmd + ["map", "-x", prefix + ".xg",
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

    cmd1 = parse(env.VG)
    cmd2 = parse(env.GRAPHVIZ_DOT)

    graphPrefix = re.sub(".vg", "", vgFile)

    find = cmd1 + ["find", "-n", str(node), 
            "-x", graphPrefix + ".xg", "-c", str(width)]
    view = cmd1 + ["view", "-dp", "-"]  

    pdfName = "_".join([graphPrefix, str(node), "w" + str(width)]) + ".pdf"
    dot = cmd2 + ["-Tpdf", "-o", pdfName]

    p1 = subprocess.Popen(find, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(view, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()  
    writer = open(pdfName, 'w+')    
    subprocess.Popen(dot, stdin=p2.stdout, stdout=writer)
    p2.stdout.close() 
    writer.close()

    return pdfName

def call_graph(vgFile, gamFile):

    cmd = parse(env.VG)

    graphPrefix = re.sub(".vg", "", vgFile)
    gamPrefix = re.sub(".gam", "", gamFile)

    packFile = gamPrefix + ".pack"
    pack = cmd + ["pack", "-g",  gamFile,
            "-x", graphPrefix + ".xg", "-o", packFile]

    subprocess.call(pack)

    call = cmd + ["call", graphPrefix + ".xg", "-k",  packFile]

    outName = gamPrefix + ".vcf"

    writer = open(outName, 'w+')
    subprocess.call(call, stdout=writer)
    writer.close()
    
    return outName

def augment_graph(vgFile, gamFile):
    
    cmd = parse(env.VG)

    augment = cmd + ["augment", "-B",  vgFile, gamFile]

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
