import sys
sys.path.append("..")
import file_handler as io
import subprocess
import pysam
import pyfaidx 
import numpy as np
import re
from matplotlib_venn import venn2, venn3
import mmh3
import matplotlib.pyplot as plt

MM2 = "/home/scott/bin/minimap2-2.17_x64-linux/minimap2"
REFERENCE_GENOME_INDEX="/media/scott/Rotom/reference/hg38/hg38.mmi" #"/media/scott/HDD/sickkids/hg38.mmi"
REFERENCE_GENOME="/media/scott/Rotom/reference/hg38/hg38.fa" #"/media/scott/HDD/sickkids/hg38.fa"
REF_ALIGNMENTS="/media/scott/Rotom/hybrid2/CF062_19/pacbio_to_canu_purged.reheader.bam" #"/media/scott/HDD/sickkids/hg38.fa"


rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

def samtools_fetch(bamFile, chrom=None, start=None, end=None):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    if chrom is None:
        alignments = [x for x in samfile.fetch()]
    else:
        alignments = [x for x in samfile.fetch(chrom, start, end)]
    return alignments


def faidx_fetch(faFile, chrom, start, end):
    fa = pyfaidx.Faidx(faFile)
    fa.fetch(chrom, start, end)
    return str(fa.fetch(chrom, start, end))

def mm2_align(refFa, readsFile, prefix, asm="asm5", secondary=False):
    
    outName = prefix + ".bam"    
    mm2Align = [MM2, "-acx", asm, "--eqx"]
    
    if not secondary:
        mm2Align.append("--secondary=no")
    
    mm2Align.extend(["-o", outName, refFa, readsFile])
    subprocess.call(mm2Align)
    
    return outName

Q="query"
R="ref"
H="hybrid"
HG="hg38"
N="_nogap"
QN=Q+N
HN=H+N

LEN="len"
ALEN="alen"
MATCH="match"
INDEL="indel"
MISMATCH="mismatch"
GAP="gap"
CLIP="clip"
NCOUNT="Ns"
GAPCOUNT="ngaps"
COMPRESSINDEL="cindel"
DV="divergence"
NM="nm"
REGION="region"
STRANDCOUNT="strandcount"
STRAND="strand"

def get_cigar_dict(dictionary, key):
    return { name : dictionary[name][key] for name in dictionary }

def cigar_summary(alignmentList, block=False):
    '''
    M   BAM_CMATCH 	0
    I 	BAM_CINS 	1
    D 	BAM_CDEL 	2
    N 	BAM_CREF_SKIP 	3
    S 	BAM_CSOFT_CLIP 	4
    H 	BAM_CHARD_CLIP 	5
    P 	BAM_CPAD 	6
    = 	BAM_CEQUAL 	7
    X 	BAM_CDIFF 	8
    B 	BAM_CBACK 	9
    NM 	NM tag 	10
    
    total, indel, clip, match
    '''
    ops = "MIDNSHP=XB"
    
    if not isinstance(alignmentList, list):
        alignmentList = [alignmentList]
        
    listSummary = dict()
    
    for op in ops: 
        listSummary[op] = 0
        
    for alignment in alignmentList:
        if alignment is None:
            continue
    
        x = alignment.get_cigar_stats()
        
        for i,op in enumerate(ops):
            listSummary[op] += x[0][i]
    
    def calc(s):
        s["total"] = sum([s[key] for key in s])
        s[INDEL] = s["I"] + s["D"] 
        s[CLIP]  = s["S"] + s["H"] 
        s[MATCH] = s["M"] + s["="] 
        s[MISMATCH] = s["X"]

        return s
    
    return calc(listSummary)


def summarize_alignments(alns, seq, name):
   
    summary = dict()
    
    summary[LEN] = len(seq)
    summary[ALEN] = sum([aln.alen for aln in alns])
    summary[GAP] = sum([aln.get_tag("nn") for aln in alns])
    summary[NCOUNT] = seq.upper().count("N")
    
    summary[STRANDCOUNT] = (sum([aln.alen for aln in alns if not aln.is_reverse]),
                       sum([aln.alen for aln in alns if aln.is_reverse]))

    summary[STRAND] = 1 if summary[STRANDCOUNT][0] > summary[STRANDCOUNT][1] else -1 
    
    summary[GAPCOUNT] = 0
    if summary[NCOUNT] > 0:
        gapCompress = re.sub("[nN]+", "N", seq)
        summary[GAPCOUNT] = gapCompress.count("N")
    

    #summary[DV] = np.mean([aln.get_tag("de") * (aln.alen)/summary[ALEN] for aln in alns])

    summary[NM] = sum([aln.get_tag("NM") for aln in alns])
    
    summary[REGION] = []
    for aln in alns:
        summary[REGION].append(aln.reference_name + ":" + 
              str(aln.reference_start) + "-" + str(aln.reference_end))

    cig = cigar_summary(alns)
    
    summary[MATCH] = cig[MATCH]
    summary[INDEL] = cig[INDEL]
    summary[CLIP] = cig[CLIP]
    summary[MISMATCH] = cig[MISMATCH]
    
    summary[COMPRESSINDEL] = aln.cigarstring.count("I") + aln.cigarstring.count("D")

    return summary
    
    print(name, ":")
    print("\talignment % :", round(100*alen/len(seq), 2))
    print("\tmatch % :", round(100*cig['=']/len(seq), 2))
    print("\tindel % :", round(100*cig['indel']/alen, 2))
    print("\tmismatch % :", round(100*cig['X']/alen, 2))
    print("\tNM :", nm)
    print("\tNNNs :", nn)
    print("\tdivergence :", round(100*de, 2))


def jaccard_similarity(a, b):
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))
    union = len(a.union(b))

    return intersection / union

def jaccard_containment(a, b):
    a = set(a)
    b = set(b)

    intersection = len(a.intersection(b))

    return intersection / len(a)

def hash_kmer(kmer):
    rc_kmer = reverse_complement(kmer)
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer

    hsh = mmh3.hash64(canonical_kmer, 42)[0]
    if hsh < 0: hsh += 2**64
    return hsh

def build_kmers(seq, k=21, excludeN=True):
    kmers = []
    n = len(seq) - k + 1

    for i in range(n):
        kmer = seq[i:i + k].upper()
        if excludeN and "N" in kmer:
            continue
        kmers.append(hash_kmer(kmer))

    return kmers

def analyze_blockpaths(blockPaths, seqData, param):

    for blockPath in blockPaths:
        print(blockPath)
        print(blockPath[0], blockPath[-1])
        qid, rid = blockPath[0].qid, blockPath[0].rid
        if qid != blockPath[-1].qid or rid != blockPath[-1].rid:
            print("Blockpath IDs don't match...")
        
        qSeq = seqData[qid][blockPath[0].qpos:blockPath[-1].qpos]
        rSeq = seqData[rid][blockPath[0].rpos:blockPath[-1].rpos]
        hSeq, source = io.path_to_sequence(blockPath, seqData)
        hSeq, source = "".join(hSeq), "".join(source)
        
        print(len(qSeq), len(rSeq), len(hSeq))
        
        seqDict = { Q: qSeq,
                    R: rSeq,
                    H: hSeq,
                    QN: re.sub("[nN]+", "", qSeq),
                    HN: re.sub("[nN]+", "", hSeq)}
        
        prefix = param.OUTPUT_DIR + "/"
                
        fa = io.dict2fasta(prefix + "seq.fa", seqDict)
        bam = mm2_align(REFERENCE_GENOME_INDEX, fa, prefix + "seq")
        alignments = samtools_fetch(bam)
        summaries = dict()
        
        for x in seqDict:
            seq = seqDict[x]
            alns = [a for a in alignments if a.qname == x]
            summaries[x] = summarize_alignments(alns, seq, x)
            
        sourcePercent = (source.count("q")/len(source), source.count("r")/len(source))
        baseCount = (summaries[QN][LEN], summaries[HN][LEN])
        
        nCount = (summaries[Q][NCOUNT], summaries[H][NCOUNT])
        gapCount = (summaries[Q][GAPCOUNT], summaries[H][GAPCOUNT])
        
        editDistance = (summaries[Q][NM], summaries[R][NM], summaries[H][NM] )

        indels = (summaries[QN][COMPRESSINDEL], summaries[R][COMPRESSINDEL], summaries[HN][COMPRESSINDEL] )
        mismatch = (summaries[QN][MISMATCH], summaries[R][MISMATCH], summaries[HN][MISMATCH] )

        print(ALEN, summaries[Q][ALEN]/summaries[Q][LEN],  summaries[R][ALEN]/summaries[R][LEN], summaries[H][ALEN]/summaries[H][LEN] )

        print(NM, summaries[Q][NM], summaries[R][NM], summaries[H][NM] )
        print(MATCH, summaries[Q][MATCH], summaries[R][MATCH], summaries[H][MATCH] )

        print(summaries[H])

        chrom=None
        start,end=None,None        
        for x in [Q,R,H]:
            regions = summaries[x][REGION]
            if len(regions) < 1:
                print(x, "has no alignment to reference genome!")
            for region in regions:
                c = region.split(":")[0]
                s,e =  region.split(":")[1].split("-")
                s,e = min(int(s), int(e)), max(int(s), int(e))
                
                if chrom is None: 
                    chrom,start,end = c,s,e
                    continue                    
                
                if chrom != c:
                    print(region, chrom, "ambiguous chromosome alignment!")
                    
                start,end=min(start,s), max(end,e)
                
        hgSeq = faidx_fetch(REFERENCE_GENOME, chrom, start, end)
        
        print("hybrid vs reference length", summaries[H][LEN], len(hgSeq))

        
        k=21
        kmers = { Q: build_kmers(qSeq, k),
                  R: build_kmers(rSeq, k),
                  H: build_kmers(hSeq, k),
                  HG: build_kmers(hgSeq, k)
                }       

        
        #print(k, jaccard_similarity(kmers[Q], kmers[R]))
        #print(k, jaccard_similarity(kmers[Q], kmers[H]))
        #print(k, jaccard_similarity(kmers[R], kmers[H]))
        
        print(k, jaccard_similarity(kmers[Q], kmers[HG]))
        print(k, jaccard_similarity(kmers[R], kmers[HG]))
        print(k, jaccard_similarity(kmers[H], kmers[HG]))

        q = set(kmers[Q])
        r = set(kmers[R])
        hg = set(kmers[HG])
        
        diffq,diffr = (q-hg),(r-hg)
        diffqr = diffq.intersection(diffr)

        intersection = len(a.intersection(b))
        
        _all = set(q).intersection(r).intersection(hg)
        _hgq =  set(q).intersection(hg) - r
        _hgr =  set(r).intersection(hg) - q
        _qr =  set(r).intersection(q) - hg
        _q =  q - r - hg
        _r =  r - q - hg
        _hg =   hg - r - q

        _unique =  _r.union(_q).union(_hg)


        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        pie = {'All':_all, 
               'Q+R':_qr, 'HG38 + Q' : _hgq, 'HG38 + R': _hgr,
               'Unique':_unique 
              }
        labels = [l for l in pie]
        sizes = [len(pie[l]) for l in pie]

        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        
        plt.show()
        
        
        venn3([set(kmers[Q]), set(kmers[R]), set(kmers[H])])
        
        #samyia22@yiustrange.com
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
'''     
        
import helper
import external_tools as tools
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def get_tag(reads, tag):
    return np.array([read.get_tag(tag) for read in reads])

def get_tag_from_dict(dictionary, tag, average=True):
    tagDict = dict()
    for name in dictionary:
        tagDict[name] = np.mean(get_tag(dictionary[name], tag))
    return tagDict
    
def join_dicts(dictionary1, dictionary2):
    list1, list2 = [], []
    list1Unmapped, list2Unmapped = [], []

    for name in dictionary1:
        if name in dictionary2:
            list1.append(dictionary1[name])
            list2.append(dictionary2[name])
        else:
            list1Unmapped.append(dictionary1[name])
    for name in dictionary2:
        if name not in dictionary1:
            list2Unmapped.append(dictionary2[name])
    
    return (np.array(list1), np.array(list2), np.array(list1Unmapped), np.array(list2Unmapped))
     
def to_array(dictionary, withNames=False):
    arrayList = []
    nameList = []
    for name in dictionary:
        arrayList.append(dictionary[name])
        nameList.append(name)

    return np.array(arrayList) if not withNames else (np.array(arrayList), np.array(nameList))
               
def to_dict(reads):
    dictionary = dict()
    for read in reads:
        name = read.qname
        if name in dictionary:
            dictionary[name].append(read)
        else:
            dictionary[name] = [read]
    return dictionary

def plot_identity(x, y):
    low, high = min(min(x), min(y)), max(max(x), max(y))
    x=np.linspace(low,high,10) 
    plt.plot(x,x,'r-')
    
def plot_scatter(x, y, alpha=0.5, colour="c"):
    plt.scatter(x, y, alpha=alpha, color=colour)
sns.set(style="whitegrid")

def plot_hist(x, bins=40, alpha=0.5, colour=None):
    if colour is not None:
        plt.hist(x, bins=bins, alpha=0.5, color=colour)
    else:
        plt.hist(x, bins=bins, alpha=0.5)

def plot_mc(originalReads, newReads, simple=True):
    originalReadDict = to_dict(originalReads)
    newReadDict = to_dict(newReads)

    originalMCDict = get_tag_from_dict(originalReadDict, "mc")
    newMCDict = get_tag_from_dict(newReadDict, "mc")

    mc1, mc2, mc1_, mc2_ = join_dicts(originalMCDict, newMCDict)

    if not simple:
        plot_identity(mc1, mc2)
        plot_scatter(mc1, mc2)
        plt.show()
    
    plot_hist(mc1, bins=30, colour="grey")
    plot_hist(mc2, bins=30, colour="forestgreen")
    plt.show()
    print("mc")
    print(np.mean(mc1), np.mean(mc2))


def plot_cigar(originalReads, newReads, simple=True):
    originalReadDict = to_dict(originalReads)
    newReadDict = to_dict(newReads)
    
    originalCigarDict = dict()
    newCigarDict = dict()

    for name in originalReadDict:
        originalCigarDict[name] = helper.cigar_summary(originalReadDict[name], block=False)
    for name in newReadDict:
        newCigarDict[name] = helper.cigar_summary(newReadDict[name], block=False)

    def get_cigar_dict(dictionary, key):
        return { name : dictionary[name][key] for name in dictionary }
    
    indel1, indel2, indel1Unmapped, indel2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "indel"), 
                                                get_cigar_dict(newCigarDict, "indel"))
    plot_identity(indel1, indel2)
    plot_scatter(indel1, indel2)
    
    total1, total2, total1Unmapped, total2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "total"), 
                                            get_cigar_dict(newCigarDict, "total"))
    
    clip1, clip2, clip1Unmapped, clip2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "clip"), 
                                            get_cigar_dict(newCigarDict, "clip"))

    norm1, norm2 = (total1 - clip1), (total2 - clip2)

        
    I1, I2, I1Unmapped, I2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "I"), 
                                                get_cigar_dict(newCigarDict, "I"))
    
    D1, D2, D1Unmapped, D2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "D"), 
                                                get_cigar_dict(newCigarDict, "D"))
    
    M1, M2, M1Unmapped, M2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "match"), 
                                                get_cigar_dict(newCigarDict, "match"))

    X1, X2, X1Unmapped, X2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "X"), 
                                                get_cigar_dict(newCigarDict, "X"))

    if not simple:

        normI1, normI2 = I1/norm1, I2/norm2
        plot_identity(normI1, normI2)
        plot_scatter(normI1, normI2)
        plt.show()
        plot_hist(normI1, colour="grey")
        plot_hist(normI2, colour="forestgreen")
        plt.show()
        print("insertion")
        print(round(np.mean(normI1),4), round(np.mean(normI2),4))
    
        normD1, normD2 = D1/norm1, D2/norm2
        plot_identity(normD1, normD2)
        plot_scatter(normD1, normD2)
        plt.show()
        plot_hist(normD1, colour="grey")
        plot_hist(normD2, colour="forestgreen")
        plt.show()
        print("deletion")
        print(round(np.mean(normD1),4), round(np.mean(normD2),4))
    
        print("")
        normX1, normX2 = X1/norm1, X2/norm2
        plot_identity(normX1, normX2)
        plot_scatter(normX1, normX2)
        plt.show()
    
    discord1, discord2 = (I1+D1+X1)/norm1, (I2+D2+X2)/norm2
    plot_identity(discord1, discord2)
    plot_scatter(discord1, discord2)
    plt.show()
    plot_hist(discord1, colour="grey")
    plot_hist(discord2, colour="forestgreen")
    plt.show()
    print("discordance")
    print(round(np.mean(discord1),4), round(np.mean(discord2),4))

    concord1, concord2 = M1/norm1, M2/norm2
    if not simple:
        plot_identity(concord1, concord2)
        plot_scatter(concord1, concord2)
        plt.show()
    plot_hist(concord1, colour="grey")
    plot_hist(concord2, colour="forestgreen")
    plt.show()
    print("concordance")
    print(round(np.mean(concord1),4), round(np.mean(concord2),4))

def calculate_difference(originalReads, newReads):
    plot_mc(originalReads, newReads)
    plot_cigar(originalReads, newReads)
    

def assess_difference(segment, seqData, param):
    
    headerBam = segment.bamA if segment.bamA is not None else segment.bamB
    headerBam = segment.bamU if headerBam is None else headerBam

    if segment.regionA is not None:
        originalFa = helper.isolate_region_fasta(segment.regionA, seqData, param)
        newFaDict = {"polished" : segment.seqA}
        prefix = helper.file_prefix(segment.regionA, param)
        newFa = tools.dict2fasta(newFaDict, prefix + "_polished_TEMP_")
        
        if segment.bamA is not None:
            originalBamA = None if segment.bamA is None else tools.align_pacbio(originalFa, segment.bamA, prefix + "_original_TEMP_A")
            newBamA = None if segment.bamA is None else tools.align_pacbio(newFa, segment.bamA, prefix + "_new_TEMP_A")
            originalReads = tools.samtools_fetch(originalBamA)
            newReads = tools.samtools_fetch(newBamA)
            calculate_difference(originalReads, newReads)
            
        if segment.bamU is not None:
            originalBamU = None if segment.bamU is None else tools.align_pacbio(originalFa, segment.bamU, prefix + "_original_TEMP_Au")
            newBamU = None if segment.bamU is None else tools.align_pacbio(newFa, segment.bamU, prefix + "_new_TEMP_Au")
            originalReads = tools.samtools_fetch(originalBamU)
            newReads = tools.samtools_fetch(newBamU)
            calculate_difference(originalReads, newReads)

    if segment.regionB is not None:
        originalFa = helper.isolate_region_fasta(segment.regionB, seqData, param)
        newFaDict = {"polished" : segment.seqB}
        prefix = helper.file_prefix(segment.regionB, param)
        newFa = tools.dict2fasta(newFaDict, prefix + "_polished_TEMP_B")
        
        if segment.bamB is not None:
            originalBamB = None if segment.bamB is None else tools.align_pacbio(originalFa, segment.bamB, prefix + "_original_TEMP_B")
            newBamB = None if segment.bamB is None else tools.align_pacbio(newFa, segment.bamB, prefix + "_new_TEMP_B")
            originalReads = tools.samtools_fetch(originalBamB)
            newReads = tools.samtools_fetch(newBamB)
            calculate_difference(originalReads, newReads)

        if segment.bamU is not None:
            originalBamU = None if segment.bamU is None else tools.align_pacbio(originalFa, segment.bamU, prefix + "_original_TEMP_Bu")
            newBamU = None if segment.bamU is None else tools.align_pacbio(newFa, segment.bamU, prefix + "_new_TEMP_Bu")
            originalReads = tools.samtools_fetch(originalBamU)
            newReads = tools.samtools_fetch(newBamU)       
            calculate_difference(originalReads, newReads)


def plot_alignments(faDict, bamDict, prefix, faDictOrder=None, image="svg"):
    
    def get_reads(fasta, bam):
        if bam is None:
            return []
        file = tools.align_pacbio(fasta, bam, prefix + "_TEMP_realigned_plot")
        reads = tools.samtools_fetch(file)
        helper.delete_file(file)
        return reads

    mc, bamName, faName, readName = [], [], [], []
    
    for faFileName in faDict:
        faFile = faDict[faFileName]
    
        for bamFileName in bamDict:
            bamFile = bamDict[bamFileName]
            reads = get_reads(faFile, bamFile)
            readDict = to_dict(reads)
            readMC,names = to_array(get_tag_from_dict(readDict, "mc"), withNames=True)
            readMC,names = list(readMC),list(names)

            mc = mc + readMC
            bamName = bamName + [bamFileName]*len(readMC)
            faName = faName + [faFileName]*len(readMC)
            readName = readName +  names
            
    plotDf = pd.DataFrame({"mc" : mc, "alignment" : bamName, "reference" : faName, "read" : readName})
    
    palette = "muted"
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style("whitegrid")
    
    data = plotDf.groupby(['alignment', 'reference']).mean().reset_index()
    #sd = plotDf.groupby(['alignment', 'reference']).std().reset_index()
    #data["mc_sd"] = sd["mc"]
    mean = np.mean(plotDf["mc"])
    
    plt.figure(figsize=(12, 8))
    ax = sns.pointplot(x="reference", y="mc",data=plotDf, palette=palette,
                         hue="alignment", size=10, linewidth=1, edgecolor='black',
                         order=faDictOrder, dodge=True, join=False, capsize=.1,
                         errwidth=2)
    ax.axhline(y = mean, color='red', linewidth=2, alpha=.4)
    ax.set(xlabel="", ylabel="Average Mapped Concordance Percentage")
    plt.legend(title='Reads')
    figure = ax.get_figure()
    figure.savefig(prefix + "_average_mapped_concordance." + image)
    
    for faFileName in faDict:

        data = plotDf[plotDf["reference"] == faFileName]
        plt.figure(figsize=(12, 8))

        sns.swarmplot("alignment", y="mc", data=data,
                 color="black", edgecolor="black", size=5, alpha=0.3)

        ax = sns.violinplot(x="alignment", y="mc", 
                             data=data, palette=palette, inner="quartile")
        ax.set_title("Reads Aligned Against " + faFileName)
        ax.set(xlabel="Reads", ylabel="Mapped Concordance Percentage")

        figure = ax.get_figure()
        figure.savefig(prefix + "_" + faFileName + "_violin." + image)

    
        plt.figure(figsize=(12, 8))
        ax = sns.boxplot(x="alignment", y="mc", data=data, palette=palette)
        ax.set_title("Reads Aligned Against " + faFileName)
        ax.set(xlabel="Reads", ylabel="Mapped Concordance Percentage")
        figure = ax.get_figure()
        figure.savefig(prefix + "_" + faFileName + "_boxplot." + image)



''' 
    