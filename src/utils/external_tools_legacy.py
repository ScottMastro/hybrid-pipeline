import subprocess

from . import environment as env
#from . import environment_hpf as env

'''
WHATSHAP = "whatshap"
HAP_PY="/home/scott/bin/hap.py-install/bin/hap.py"
PRE_PY="/home/scott/bin/hap.py-install/bin/pre.py"
RTG = "/home/scott/bin/rtg-tools-3.10.1/rtg"
NANOPLOT = "NanoPlot"
GRAPHVIZ_DOT="dot"
'''

def parse(command):
    return command.split()

def unalign_bam(bamFile, prefix):    
    
    cmd = parse(env.PICARD)
    revert = cmd + ["RevertSam"]
    revert.append("I=" + bamFile)
    unaligned = prefix + ".unaligned.bam"
    revert.append("O=" + unaligned)
    
    subprocess.call(revert)
    return unaligned 

def samtools_subset(bamFile, region, prefix):
    '''View and index a subset of a bam file'''
    
    cmd = parse(env.SAMTOOLS)
    outName = prefix + ".subset.bam"

    writer = open(outName, 'w+')    
    subprocess.call(cmd + ["view", "-b", bamFile, str(region)], stdout=writer)
    writer.close()
        
    samtools_index(outName)
    return outName

def bam2fastq(bamFile):

    cmd = parse(env.SAMTOOLS)

    outName = ".".join(bamFile.split(".")[:-1]) + ".fq"    
    writer = open(outName, 'w+')    
    subprocess.call(cmd + ["bam2fq", bamFile], stdout=writer)
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

def sort_vcf(vcfFile, removeOriginal=False, index=False):
    cmd = parse(env.BCFTOOLS)
    
    outName = re.sub(".vcf", ".sorted.vcf", vcfFile)
    sort = cmd + ["sort"]
    if index:
        sort.append("-Ob")
        outName.append(".gz")
        
    sort = sort + ["-o", outName, vcfFile]

    print(sort)
    subprocess.call(sort)
    if removeOriginal:
        os.remove(vcfFile)
    return outName

def whatshap_genotype(bamFile, vcfFile, refFa, prefix, indels=True):
    
    cmd = parse(env.WHATSHAP)
    outName = prefix + ".whatshap.genotype.vcf"

    whatshapGenotype = cmd + ["genotype", "-o", outName, "--ignore-read-groups"]
    
    whatshapGenotype.extend(["--reference", refFa])
    if indels:
        whatshapGenotype.append("--indels")

    whatshapGenotype.extend([vcfFile, bamFile])
    subprocess.call(whatshapGenotype)
    return outName

def find_sv(bamFile, refFa, prefix, threads=6, indelOnly=False):

    cmd = parse(env.PBSV)

    outName = prefix + ".pbsv.vcf"
    signatureFile = prefix + "pbsv.svsig.gz"

    discover = cmd + ["discover", "--sample", "hybrid", bamFile, signatureFile]
    call = cmd + ["call", "-j", str(threads)]

    if indelOnly:
        call.extend(["--types", "DEL,INS"])
        
    call.extend([refFa, signatureFile, outName])

    subprocess.call(discover)
    subprocess.call(call)
    return outName

def run_nucdiff(refFa, queryFa, outDir, filePrefix=None):
    
    cmd = parse(env.NUCDIFF)
    
    if filePrefix is None:
        filePrefix="nucdiff"
        
    nucdiff = cmd + [refFa, queryFa, outDir, filePrefix]    
    subprocess.call(nucdiff)
            
    return outDir

def get_nucdiff_stats(outDir, filePrefix=None, toInt=False):
    
    if filePrefix is None:
        filePrefix="nucdiff"
    
    statsFile = outDir+"/results/" + filePrefix + "_stat.out"
    reader = open(statsFile, "r")
    
    def to_int(x):
        if not toInt: return x
        try: 
            return int(x)
        except:
            return x
    
    def clean(x):
        x = re.sub("-", "_", x)
        x = "_".join(x.split())
        x = x.lower()
        return x
    
    d = dict()
    for line in reader:
        split = line.strip().split("\t")
        if len(split) != 2: continue
        key, value = split        
        d[clean(key)] = to_int(value)
        
    reader.close()
            
    return d

def canu_correct(fastaFile, prefix, size, trim=True, rawErrorRate=None):
    cmd = parse(env.CANU)
    name = prefix.split("/")[-1]
    
    canuCorrect = cmd + ["-correct", "-p", name, "-d", prefix]
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
    
    cmd = parse(env.CANU)
    name = prefix.split("/")[-1]

    canuAssemble = cmd + ["-assemble", "-p", name, "-d", prefix]
    canuAssemble.append("genomeSize=" + str(size/1000) + "k")
    canuAssemble.append("-pacbio-corrected")
    canuAssemble.append(fastaFile)

    subprocess.call(canuAssemble)
   
def hap_py(refFa, vcf1, vcf2, prefix, outPrefix="happy", engine="xcmp"):
    cmd = parse(env.HAP_PY)

    outDir=prefix + "happy/"
    try:
        os.mkdir(outDir)
    except: 
        pass
    
    happy = cmd + ["-r", refFa, "-L", "--preprocess-truth", "--engine", engine]
    
    if engine == "vcfeval":
        happy.append("--engine-vcfeval-path")
        happy.append(env.RTG)

    happy = happy + ["-o", outDir + outPrefix]    
    happy.extend([vcf1, vcf2])
    
    print(happy)
    subprocess.call(happy)
    
    return outDir + outPrefix + ".vcf.gz"

def pre_py(refFa, vcfFile, prefix=None, decompose=True, noGz=False):
    cmd1 = parse(env.PRE_PY)
    cmd2 = parse(env.BGZIP)

    if prefix is None:
        outname = re.sub(".gz", "", vcfFile)
        outName = re.sub(".vcf", ".normalized.vcf.gz", outname)
    else:
        outName = prefix + ".vcf.gz"

    prepy = cmd1 + ["-r", refFa]
    if decompose:
        prepy.append("--decompose")
    prepy.extend([vcfFile, outName])
    subprocess.call(prepy)
    
    if noGz:
        newName = re.sub(".gz", "", outName)
        unzip = cmd2 + ["-d", outName]
        subprocess.call(unzip)
        return newName

    return outName

def make_sdf(refFa, sdfName):  
    cmd = parse(env.RTG)
    rtgformat = cmd + ["format", "-o", sdfName, refFa]
    subprocess.call(rtgformat)
    
    return sdfName

def vcfeval(refFa, baselineVCF, callsVCF, prefix, dirname="vcfeval_out", 
            squashPloidy=False, useQUAL=True, useFilter=False, outputType="split"):
    cmd = parse(env.RTG)

    #rtgtools vcfeval -b consensus-longranger/phased_variants.vcf.gz -c AB_temp_merged.vcf -o rtg
    vcfeval = cmd + ["vcfeval", "-b", baselineVCF, "-c", callsVCF]

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

