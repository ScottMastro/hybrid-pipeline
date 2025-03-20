import mmh3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as PathEffects
import matplotlib.colors as mc
import colorsys
import external_tools as tools

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

def compress_homopolymer(seq):
    if len(seq) < 1: return ""
    compressed = [seq[0]]
    prev = seq[0]
    
    for s in seq[1:]:
        if s != prev: 
            compressed.append(s)
            prev = s
            
    return "".join(compressed)

def canonical_kmer(kmer):
    rc_kmer = reverse_complement(kmer)
    return kmer if kmer < rc_kmer else rc_kmer

def hash_kmer(kmer):
    hsh = mmh3.hash64(canonical_kmer(kmer), 42)[0]
    if hsh < 0: hsh += 2**64
    return hsh

class KmerSet:    
    def __init__(self, seq, k=21, excludeN=True, compressed=False):
        self.k = k
        self.kmers = set()
        self.compressed = compressed
        self.seqLen = None

        if seq is None: return
        
        seq = seq if not compressed else compress_homopolymer(seq)
        self.seqLen = len(seq)
        
        n = len(seq) - k + 1
        for i in range(n):
            kmer = seq[i:i + k]
            if excludeN and "N" in kmer: continue
            self.kmers.add(hash_kmer(kmer))
                    
    def jaccard_similarity(self, other):
        a, b = self.kmers, other.kmers

        intersection = len(a.intersection(b))
        union = len(a.union(b))
        return intersection / union

    def jaccard_containment(self, other, compressed=False):
        a, b = self.kmers, other.kmers
    
        intersection = len(a.intersection(b))
        return intersection / len(a)

    def __sub__(self, other): #set diff
        new = KmerSet(None, self.k)
        new.kmers = self.kmers - other.kmers
        return new
   
    def __add__(self, other):  #union
        new = KmerSet(None, self.k)
        new.kmers = self.kmers.union(other.kmers)
        return new

    def __and__(self, other):  #intersection
        new = KmerSet(None, self.k)
        new.kmers = self.kmers.intersection(other.kmers)
        return new

    def __len__(self):
        return len(self.kmers)
    def __repr__(self):
        string = ""
        if self.seqLen is not None:
            string = string + ("Compressed s" if self.compressed else "S") + \
                "equence of length " + str(self.seqLen) + " with " 
                
        string = string + str(len(self.kmers)) + " kmers, k = " + str(self.k)
        return string
    
    def __str__(self):
        return self.__repr__()

def kmers_from_bam(bamFile, k=21):
    
    alignments = tools.samtools_fetch(bamFile)
    kmersA = dict()
    kmersB = dict()
    kmersU = dict()

    def add_kmers(kmers, kdict):
        for kmer in kmers:
            if kmer in kdict: kdict[kmer] += 1
            else: kdict[kmer] = 1

    for alignment in alignments:
        qseq = alignment.query_sequence
        n = len(qseq) - k + 1
        kmers = [qseq[i:i + k] for i in range(n)]
        kmers = [canonical_kmer(kmer) for kmer in kmers]

        
        try:
            hp = alignment.get_tag("HP")
            if hp == 1:
                add_kmers(kmers, kmersA)
            if hp == 2:
                add_kmers(kmers, kmersB)            
        except:
            add_kmers(kmers, kmersU)

    return kmersA, kmersB, kmersU


def lighten_color(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def draw_pie(kmerSetTuple, showPc=True, total=1, plotPath=None):
            
    fig, ax = plt.subplots(figsize=(12, 9), subplot_kw=dict(aspect="equal"))

    labels = [tup[0] for tup in kmerSetTuple]
    values = [tup[1] for tup in kmerSetTuple]
    colours = [tup[2] for tup in kmerSetTuple]

    def func(pct, allvals):
        realPct = pct*total
        return "{:.1f}%".format(realPct)
    
    if showPc:
        wedges, texts, autotexts = ax.pie(values, colors=colours,
                                          wedgeprops=dict(width=0.5),
                                          startangle=-40, pctdistance=0.75,
                                          autopct=lambda pct: func(pct, values),
                                          textprops={"color":"k"})
    else:
        wedges, texts = ax.pie(values, colors=colours, pctdistance=0.75,
                               wedgeprops=dict(width=0.5), startangle=-40)

    
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")
    
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.1*np.sign(x), 1.2*y),
                    horizontalalignment=horizontalalignment, size=20, **kw)
    if showPc:
        plt.setp(autotexts, size=28) #, weight="bold")
        for i,text in enumerate(autotexts):
            text.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='#424242')])
            plt.setp(text, color=lighten_color(colours[i]))

    if plotPath is not None:
        plt.savefig(plotPath)
        plt.close(fig)
    else:
        plt.show()

def draw_pie2(name1, kmer1, name2, kmer2, plotPath, fileType=".svg"):
    
    data = []
    shared = kmer1 & kmer2
    only1 = kmer1 - kmer2
    only2 = kmer2 - kmer1
    data.append(("Shared", len(shared), "#a4a698"))
    data.append((name1 + " only", len(only1), "#a7dbd7"))
    data.append((name2 + " only", len(only2), "#f38631"))

    draw_pie(data, plotPath=plotPath + fileType)

    info=dict()
    info["kmer_overlap_" + name1 + "_" + name2] = len(kmer1 & kmer2)
    info["kmer_overlap_" + name1] = len(only1)
    info["kmer_overlap_" + name2] = len(only2)
    return info


def draw_pie3_counts(name1, name2, name3,
                         shared, shared12, shared13, shared23, 
                         only1, only2, only3, plotPath, fileType=".svg"):
    
    data1, data2 = [], []
    data1.append(("Shared by all", shared, "#a4a698"))
    notShared = shared12 + shared13 + shared23 + only1 + only2 + only3
    data1.append(("Other", notShared, "#e0e4cd"))
    
    red="#fc3912"
    blue="#12b2fc"
    yellow= "#fff355"
    
    blueYellow="#98b797"
    redBlue= "#a298b1"
    redYellow="#cda675"

    data2.append((name1 + " only", only1, red))
    data2.append((name1 + " and " + name2, shared12, redBlue))
    data2.append((name2 + " only", only2, blue))
    data2.append((name2 + " and " + name3, shared23, blueYellow))
    data2.append((name3 + " only", only3, yellow))
    data2.append((name1 + " and " + name3, shared13, redYellow))

    draw_pie(data1, plotPath=plotPath + "_1" + fileType)
    portion = notShared/(notShared + shared)
    draw_pie(data2, total=portion, plotPath=plotPath + "_2" + fileType)
    
    info=dict()
    info["kmer_overlap_" + name1 + "_" + name2 + "_" + name3] = shared
    info["kmer_overlap_" + name1 + "_" + name2] = shared12
    info["kmer_overlap_" + name1 + "_" + name3] = shared13
    info["kmer_overlap_" + name2 + "_" + name3] = shared23
    info["kmer_overlap_" + name1] = only1
    info["kmer_overlap_" + name2] = only2
    info["kmer_overlap_" + name3] = only3
    
    return info

def draw_pie3(name1, kmer1, name2, kmer2, name3, kmer3, plotPath, fileType=".svg"):
    
    shared = kmer1 & kmer2 & kmer3

    u1 = kmer1 - shared
    u2 = kmer2 - shared
    u3 = kmer3 - shared

    only1 = u1 - u2 - u3
    only2 = u2 - u1 - u3
    only3 = u3 - u1 - u2

    shared12 = u1 & u2
    shared13 = u1 & u3
    shared23 = u2 & u3
    
    return draw_pie3_counts(name1, name2, name3,
                         len(shared), len(shared12), len(shared13), len(shared23), 
                         len(only1), len(only2), len(only3), plotPath, fileType)