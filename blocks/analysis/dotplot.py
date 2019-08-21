import flexidot3 as dot
import glob, os, re
import contig_output as output
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from pyfaidx import Fasta
from reverse_complement import reverse_complement
from itertools import islice
from collections import defaultdict
import numpy as np

blastdir = "/home/scott/blast/"
blastdb = "hg38"
hg38fasta = "hg38.fa.gz"

outdir = "/home/scott/Dropbox/hybrid-pipeline/blocks/analysis/dot/"

hg38 = None
QUERY, REF, HYBRID, HREF = "query", "reference", "hybrid", "hg38"

def build_temp_fasta(ids, seqs, path):
    f = open(path, "w")

    for fid, seq in zip(ids, seqs):
        f.write(">" + fid + "\n")
        f.write(re.sub("(.{64})", "\\1\n", "".join(seq), 0, re.DOTALL) + "\n")
    f.close()  
    
def build_temp_gff(seqNames, annotations, ranges, path):
    f = open(path, "w")
    t = "\t"
    i = 1
    for name, annotation, (start, end) in zip(seqNames, annotations, ranges):
        f.write(str(name) +t+ "manual_annotations" +t+ str(annotation) +t+ \
                str(start) +t+ str(end) +t+ "." +t+ "+" +t+ "." +t+ "ID=" + str(i) + "\n")
        i = i + 1
    f.close()   

def build_temp_gff_config(colorDict, alphaDict, path):
    f = open(path, "w")
    t = "\t"
    
    f.write("#annotation_type" +t+ "color" +t+ "alpha" +t+ "zoom" + "\n")
    for annotation in colorDict.keys():
        f.write(str(annotation) +t+ str(colorDict[annotation]) +t+  \
                str(alphaDict[annotation]) +t+ "0" + "\n")
        
    f.close()   

def remove_temp_file(path): os.remove(path)
    
def search_reference(sequence, hg38=hg38, blastdir=blastdir):

    #bgzip hg38 fasta file
    if hg38 is None:
        hg38 = Fasta(blastdir + hg38fasta, as_raw=True, sequence_always_upper=True)
        
    #makeblastdb -in hg38.fa.gz -dbtype nucl -title hg38.p12 -out hg38
    tempblastout = blastdir + "__blast_out__.xml"
    blastquery = blastdir + "__blast__.fa"

    ids = ["blast_start", "blast_end"]
    seqs = [sequence[:200], sequence[-200:]]
    build_temp_fasta(ids, seqs, path=blastquery)
    
    blast = NcbiblastnCommandline(query=blastquery, db=blastdir + blastdb, evalue=0.001,
                                      outfmt=5, out=tempblastout)
    stdout, stderr = blast()
    
    results = open(tempblastout)
    remove_temp_file(blastquery)
    os.remove(tempblastout)

    start_record, end_record = NCBIXML.parse(results) 
    
    if len(start_record.alignments) < 1 or len(end_record.alignments) < 1:
        print("BLAST alignment failed. Missing alignment.")
        return (None, None, None, None)
    
    startChr = start_record.alignments[0].hit_def
    endChr = end_record.alignments[0].hit_def
    
    if startChr != endChr:
        print("BLAST alignment failed. Different chromsomes.")
        return (None, None, None, None)
    
    start = start_record.alignments[0].hsps[0].sbjct_start
    end = end_record.alignments[0].hsps[0].sbjct_end
    
    if end < start: start, end = end, start
    return (hg38[startChr][start:end], startChr, start, end)
    

    
def generate_gff(seqs, nameDict, hybridSourceList, gff, gffconfig):
   
    names, ranges, sources = [], [], []
    pos=0
    for x in hybridSourceList:
        size = len(x)
        if size > 0:
            sources.append(QUERY if x[0]=="q" else REF)
            names.append(HYBRID)
            ranges.append((pos+1, pos+size+1))
            pos = pos + size + 1
        
    names.append(nameDict[QUERY])
    sources.append(QUERY)
    ranges.append( (1, len(seqs[0])) )
    names.append(nameDict[REF])
    sources.append(REF)
    ranges.append( (1, len(seqs[1])) )
    names.append(nameDict[HREF])
    sources.append(HREF)
    ranges.append( (1, len(seqs[3])) )

    colorDict = {QUERY:"skyblue", REF:"coral", HREF:"orchid"}
    a = 0.7
    alphaDict = {QUERY:a, REF:a, HREF:a}

    build_temp_gff(names, sources, ranges, path=gff)
    build_temp_gff_config(colorDict, alphaDict, path=gffconfig)


def get_kmers(sequence, k):
    it = iter(sequence)
    result = tuple(islice(it, k))
    if len(result) == k:
        yield "".join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield "".join(result)

def generate_matrix(ids, seqs, path, k=7):
    kmerCounts = dict()
    
    for idx, seq in enumerate(seqs):
        kmers = defaultdict(int)
        for kmer in get_kmers(seq, k): 
            if 'N' in kmer: continue
            kmers[kmer] += 1
        kmerCounts[idx] = kmers
    
    kmerSet = set()
    for idx in range(len(seqs)):
        for kmer in kmerCounts[idx].keys(): kmerSet.add(kmer)
            
    kmerSpectra = dict()
    for idx in range(len(seqs)):
        kmerSpectra[idx] = [0 if kmer not in kmerCounts[idx] else kmerCounts[idx][kmer] for kmer in kmerSet]
        kmerSpectra[idx] = np.array(kmerSpectra[idx])
        kmerSpectra[idx] = kmerSpectra[idx] / np.sum(kmerSpectra[idx])
        
    f = open(path, "w")
    t = "\t"

    f.write(t + t.join([ids[idx] for idx in range(len(seqs))]) + '\n')
    
    def corr(x, y):
        return str(np.round(np.corrcoef(x,y)[1,0], 2))

    for i in range(len(seqs)):
        f.write(ids[i] + t + t.join( \
                [ corr(kmerSpectra[i], kmerSpectra[j]) for j in range(len(seqs))])
                + '\n')

    f.close()   
    
def draw_dotplot(nameDict, seqs, hybridSourceList, outdir=outdir, prefix=None, wordSize=None):

    filePrefix = nameDict[QUERY] + "_" + nameDict[REF] + "_" + nameDict[HREF]
    if prefix is not None: filePrefix = prefix
    
    print("Writing to FASTA...")
    tempfasta = outdir + "__temp__.fa"
    tempgff = outdir + "__temp__.gff"
    tempgffconfig = outdir + "__temp__.config"
    tempmatrix = outdir + "__temp__.matrix.txt"

    ids = [nameDict[QUERY], nameDict[REF], nameDict[HYBRID], nameDict[HREF]]

    build_temp_fasta(ids, seqs, path=tempfasta)
    generate_gff(seqs, nameDict, hybridSourceList, tempgff, tempgffconfig)
    generate_matrix(ids, seqs, tempmatrix)
    
    if wordSize is None:
            wordSize = int(max(7, min(250, max([len(x) for x in seqs])/1000)))

    print("Building dotplot...")
    dot.main(tempfasta, wordSize, prefix=outdir, lcs_shading_ori=2, user_matrix_print=True,  \
             gff=tempgff, gff_color_config_file=tempgffconfig, input_user_matrix_file=tempmatrix)
    file1=outdir + "-Polydotplot_LCS_Shading_Legend_*.png"
    file2=outdir + "-Polydotplot_lcs_data_file_6shades_ref0_ori2.txt"
    file3=outdir + "-Polydotplot_wordsize" + str(wordSize) + "_matrix_6shades_ref0_ori2.png"
    
    for x in glob.glob(file1): os.remove(x)
    for x in glob.glob(file2): os.remove(x)

    #os.rename(file2, outdir + "data" + "-" + refString + "-ws" + str(wordSize) + ".txt")
    os.rename(file3, outdir + filePrefix + "_ws" + str(wordSize) + ".png")
    remove_temp_file(tempfasta)
    remove_temp_file(tempgff)
    remove_temp_file(tempgffconfig)
    remove_temp_file(tempmatrix)

    
def plot_gap(startFork, endFork, seqData, lengthData, \
                     out=outdir, prefix=None, wordSize=None, buffer=1000):
    nameDict = dict()
    
    def get_seq(q=True):
        s = endFork.get_pos_norm(lengthData, q=q)
        e = startFork.get_pos_norm(lengthData, q=q)
        if s > e: s,e=e,s
        tigId = startFork.get_id(q)
        s = max(0, s - buffer)
        e = min(lengthData[tigId], e + buffer)
        nameDict[QUERY if q else REF] = str(tigId).replace("_pilon", "")+ ":" + str(s) + "-" + str(e)
        strand = startFork.get_strand(q=q)
        seq = seqData[tigId][s:e]
        return reverse_complement(seq) if strand == -1 else seq
    
    print("Extracting query and reference sequence...")
    seqs = [get_seq(q=True), get_seq(q=False)]
    
    print("Extracting hybrid sequence...")
    hybridSourceList = ["q"*buffer, "r"*(len(seqs[1])-buffer*2), "q"*buffer]
    seqs.append(seqs[0][:buffer] + seqs[1][buffer:len(seqs[1])-buffer] + seqs[0][-buffer:])
    nameDict[HYBRID] = "hybrid"
    
    print("BLASTing against hg38...")
    seq, chrom, start, end = search_reference(seqs[0])
    if seq is None: return
    nameDict[HREF] = chrom + ":" + str(start) + "-" + str(end)
    seqs.append(seq)

    draw_dotplot(nameDict, seqs, hybridSourceList, outdir=out, prefix=prefix, wordSize=wordSize)

def plot_path(path, seqData, lengthData, \
                      out=outdir, prefix=None, wordSize=None, buffer=0):

    nameDict = dict()
    if len(path) < 2: return
    
    def get_seq(q=True):
        pos = np.array([f.get_pos_norm(lengthData, q=q) for f in path])
        s, e = np.min(pos), np.max(pos)
        tigId = path[0].get_id(q)
        s = max(0, s - buffer)
        e = min(lengthData[tigId], e + buffer)
        nameDict[QUERY if q else REF] = str(tigId).replace("_pilon", "")+ ":" + str(s) + "-" + str(e)
        return seqData[tigId][s:e]
    
    print("Extracting query and reference sequence...")
    seqs = [get_seq(q=True), get_seq(q=False)]
    hybridSeqList, hybridSourceList = output.path_to_sequence2(path, seqData)
    nameDict[HYBRID] = "hybrid"

    seqs.append("".join(hybridSeqList))
    
    print("BLASTing against hg38...")
    seq, chrom, start, end = search_reference(seqs[2])
    if seq is None: return
    nameDict[HREF] = chrom + ":" + str(start) + "-" + str(end)
    seqs.append(seq)
    
    draw_dotplot(nameDict, seqs, hybridSourceList, outdir=out, prefix=prefix, wordSize=wordSize)    
    
    
def gaps_dotplot(blockPaths, seqData, lengthData, outdir=outdir):
    
    prevPath = None
    for pid,path in enumerate(blockPaths):
        
        prevFork = None
        for fid,fork in enumerate(path):
            prefix = fork.get_id(q=True) + "_" + fork.get_id(q=False).replace("_pilon", "") + "_" + \
                    "path" + str(pid) + "_" + "fork" + str(fid-1) + "-" + str(fid) 
            if prevFork is not None and prevFork.is_switch_reference():
                plot_gap(prevFork, fork, seqData, lengthData, prefix=prefix)
                    
            prevFork = fork
            
        prefix = fork.get_id(q=True) + "_" + fork.get_id(q=False).replace("_pilon", "") + "_" + \
            "path" + str(pid-1) + "-" + str(pid) 
        if prevPath is not None:
            plot_gap(prevPath[-1], path[0], seqData, lengthData, prefix=prefix, buffer=2500)
    
    prevPath = path
            
def paths_dotplot(blockPaths, seqData, lengthData, outdir=outdir, gaps=True, paths=True):
    
    for pid,path in enumerate(blockPaths):
        prefix = path[0].get_id(q=True) + "_" + path[0].get_id(q=False).replace("_pilon", "") + \
            "_path" + str(pid)
        plot_path(path, seqData, lengthData, prefix=prefix, buffer=0)

            
def megapath_dotplot(megaPath, seqData, lengthData, outdir=outdir):

        prefix = megaPath[0].get_id(q=True) + "_" + megaPath[0].get_id(q=False).replace("_pilon", "") + \
            "_megapath"
        plot_path(megaPath, seqData, lengthData, prefix=prefix, buffer=0)
