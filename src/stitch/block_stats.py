import re
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("..")

import utils.log as logger
from structures.region import SimpleRegion
import dependencies.external_tools as tools
import utils.fasta_handler as fasta

Q,R,H,HG = "query", "ref", "hybrid", "hg38"

def N50(lengths):
    lengths.sort(reverse=True)
    lenSum = sum(lengths)
    
    target = lenSum/2
    cumSum = 0
    index = None
    for i,l in enumerate(lengths):
        cumSum += l
        if cumSum > target:
            index = i
            break
        
    if index is None:
        return None
    return lengths[index]

def parse_interval_bed(fileName, param, refKey=True):
    qreader = open(param.OUTPUT_DIR + "/q_" + fileName, "r")
    rreader = open(param.OUTPUT_DIR + "/r_" + fileName, "r")

    def parse_block(qline, rline):
        qparts = qline.split("\t")
        rparts = rline.split("\t")
        q = SimpleRegion(qparts[0], int(qparts[1]), int(qparts[2]))
        r = SimpleRegion(rparts[0], int(rparts[1]), int(rparts[2]))

        name = qparts[3]
        d = dict()
        for string in name.split("_")[1:]:
            key, value = string.split("=")
            d[key] = float(value)

        return (q,r,d)
    
    regionBlocks = [parse_block(qline, rline) for qline, rline in zip(qreader, rreader)]
    qreader.close() ; rreader.close()

    blockDict = dict()
    if refKey:
        for q,r,d in regionBlocks:
            if r.chrom not in blockDict:
                blockDict[r.chrom] = []    
            blockDict[r.chrom].append((q,r,d))
    else:
        for q,r in regionBlocks:
            if q.chrom not in blockDict:
                blockDict[q.chrom] = []    
            blockDict[q.chrom].append((q,r,d))

    return blockDict

def sequence_assessment(seq, prefix):
    
    def gap_count(seq):
        gapCompress = re.sub("[N]+", "N", seq)
        return gapCompress.count("N")
    def gc_count(seq):
        return (seq.count("G") + seq.count("C")) / (len(seq) - seq.count("N"))
    
    info = dict()
        
    info["length_" + prefix] = len(seq)
    info["ambiguous_bases_" + prefix] = seq.count("N")
    info["gaps_" + prefix] = gap_count(seq)
    info["gc_percent_" + prefix] = gc_count(seq)

    return info

def analyze_blocks(fileName, seqData, lengthData, param):
    
    infos = []
    blockDict = parse_interval_bed(fileName, param)
    
    for tigId in blockDict:
        print(tigId)
        for qRegion, rRegion, propDict in blockDict[tigId]:
                            
            info = dict()
            info.update(propDict)
            #print(info['cov']/len(qRegion), info['cov']/len(qRegion))
            qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
            rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    
            info.update(sequence_assessment(qSeq, Q))
            info.update(sequence_assessment(rSeq, R))
            infos.append(info)

            
    def vector(key):
        return np.array([info[key] for info in infos])
    qlens = vector("length_" + Q)
    rlens = vector("length_" + R)
    pcid = vector("id")
    coverage = vector("cov")


    plt.hist(pcid, density=True, bins=100)
    plt.xlabel("percent identity")
    
    plt.hist(np.log10(coverage), density=True, bins=50)
    plt.xlabel("coverage")

    plt.scatter(np.log10(qlens), np.log10(coverage), alpha=0.1)
    plt.xlabel("coverage")
    
    plt.hist(coverage/qlens, bins=100)
    plt.xlabel("coverage")
    plt.hist(coverage/rlens, bins=100)
    plt.xlabel("coverage")

    plt.hist(np.log10(qlens), density=True, bins=30)
    plt.xlabel("log10 qlen")
    plt.hist(np.log10(rlens), density=True, bins=30)
    plt.xlabel("log10 rlen")

    print(N50(qlens))
    print(N50(rlens))

            
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]

    
    info.update(sequence_assessment(qSeq, Q))
    info.update(sequence_assessment(rSeq, R))
    
    print(info)
    
    '''
    qFa = fasta.dict2fasta(param.OUTPUT_DIR + "/ref", {Q:qSeq})
    rFa = fasta.dict2fasta(param.OUTPUT_DIR + "/qur" + rRegion.clean_str(), {R:rSeq})
    diffDir = tools.run_nucdiff(rFa, qFa, param.OUTPUT_DIR + "/nucdiff", filePrefix=None)
    stats = tools.get_nucdiff_stats(diffDir, filePrefix=None, toInt=True)            
    print([(x, stats[x]) for x in stats if stats[x] > 0])
    input()
    '''

            

        

        
    if len(blockPaths) > 1:
        print("Analyzing inter-block sequence...")
        for i,blockPath in enumerate(blockPaths[1:]):
            prevBlockPath = blockPaths[i]
            analyze_inter_block(prevBlockPath, blockPath, seqData, lengthData, param, skipIfDone=True)







class AlignmentInfo:

    def _parse_cigar(self, aln):
        ops = "MIDNSHP=XB"
        cigars = aln.get_cigar_stats()
        return { op : cigars[0][i] for i,op in enumerate(ops) }       

    def _get_tag(self, tag, alns):
        try: 
            return [aln.get_tag(tag) for aln in alns]
        except:
            return None

    def __init__(self, alns):
        
        if not isinstance(alns, list):
            alns = [alns]
            
        alns = [aln for aln in alns if aln.reference_name is not None]

        self.alignments = len(alns)
        self.qlen = [aln.qlen for aln in alns]
        self.alen = [aln.alen for aln in alns]
        self.qsize = alns[0].infer_query_length() if len(alns) > 0 else 0
        self.strand = [-1 if aln.is_reverse else 1 for aln in alns]
        self.regions = [rgn.SimpleRegion(aln.reference_name, aln.reference_start, aln.reference_end) for aln in alns]                    

        self.cigars = [self._parse_cigar(aln) for aln in alns]
        self.compressedIndels = [aln.cigarstring.count("I") + aln.cigarstring.count("D") for aln in alns]
        self.concordance = self._get_tag("mc", alns)
        self.editDist = self._get_tag("NM", alns)

    def aligned_length(self):
        return sum(self.alen)
    def query_aligned_length(self):
        return sum(self.qlen)

    def indel_count(self):
        return sum([cig["I"] + cig["D"] for cig in self.cigars])
    def match_count(self):
        return sum([cig["M"] + cig["="] for cig in self.cigars])
    def mismatch_count(self):
        return sum([cig["X"] for cig in self.cigars])
    def clip_count(self):
        return sum([cig["S"] + cig["h"] for cig in self.cigars])

    def weighted_concordance(self):
        return sum([c*l for c,l in zip(self.concordance, self.qlen)])/self.query_aligned_length()

    def get_strand_count(self):
        plus, minus = 0, 0
        for s, alen in zip(self.strand, self.alen):
            plus += 0 if s < 0 else alen
            minus += 0 if s > 0 else alen
            
        return (plus, minus)
   
    def contiguous_region(self, misalignTolerance=0.95, sizeTolerance=1.25):        
        
        chromDict = {region.chrom: [] for region in self.regions}
        for region in self.regions: chromDict[region.chrom].append(region)
        totalSize = sum([len(region) for region in self.regions])
                
        #ignore small alignments to other chroms
        choice = None
        for chrom in chromDict: 
            percentage = sum([len(region) for region in chromDict[chrom]])/totalSize
            if  percentage > misalignTolerance:
                choice = chrom ; break
        
        if choice is None: return False
        chromRegion = chromDict[choice]
        chrom = choice

        start = min([region.start for region in chromRegion])
        end = max([region.end for region in chromRegion])
        rsize = end - start
        
        #ignore far-away small alignments
        while(True):

            if rsize > self.qsize*sizeTolerance:
                excludeStart = [region for region in chromRegion if region.start != start]
                excludeStartSize = max([region.end for region in excludeStart]) - min([region.start for region in excludeStart])
                excludeEnd = [region for region in chromRegion if region.end != end]
                excludeEndSize = max([region.end for region in excludeEnd]) - min([region.start for region in excludeEnd])
    
                #try removing start alignment
                if excludeStartSize < excludeEndSize:
                    s = sum([len(region) for region in excludeStart])
                    if s/totalSize > misalignTolerance:
                        chromRegion = excludeStart
                        start = min([region.start for region in chromRegion])
                        end = max([region.end for region in chromRegion])
                        rsize = end - start
                        if rsize >= self.qsize/sizeTolerance:
                            continue
                    
                #try removing end alignment
                s = sum([len(region) for region in excludeEnd])
                if s/totalSize > misalignTolerance:
                    chromRegion = excludeEnd
                    start = min([region.start for region in chromRegion])
                    end = max([region.end for region in chromRegion])
                    rsize = end - start
                    if rsize >= self.qsize/sizeTolerance:
                        continue

                return False
            
            else:
                break
            
        #print(rsize, self.qsize*tolerance, self.qsize/tolerance)
        if rsize <= self.qsize*sizeTolerance and rsize >= self.qsize/sizeTolerance:
            return rgn.SimpleRegion(chrom, start, end)
        return False
    
    #def __len__(self):
    #def __repr__(self):
    #def __str__(self):



def make_read_dict(alns):
    readDict = {aln.qname : [] for aln in alns}
    for aln in alns:
        readDict[aln.qname].append(aln)
    return readDict
  
def pacbio_align_assessment(seqDict, rRegion, outdir, param):
    
    alignments = tools.samtools_fetch(param.REF_ALIGNED_READS, rRegion)
    
    def coverage(alignDict, length):
        return sum([sum([aln.alen for aln in alignDict[id]]) for id in alignDict]) / length
    
    readsBam = tools.samtools_write(alignments, outdir + "pacbio_alignments", param.REF_ALIGNED_READS)
    bams,alns,readDict = dict(), dict(), dict()

    for x in seqDict:
        xfa = io.dict2fasta(outdir + x + "_seq.fa", {x:seqDict[x]})
        bams[x] = tools.align_pacbio(xfa, readsBam, outdir + x + "_pacbio")
        alns[x] = tools.samtools_fetch(bams[x])
        readDict[x] = make_read_dict(alns[x])
        
    QUALITY_THRESHOLD = 80 #%
    SIZE_THRESHOLD = 30 #%
    EDGE_SIZE = 2000 #bp
    EDGE_THRESHOLD = 90 #%

    readQualityDict,sizePercentDict,isEdgeRead = dict(),dict(),dict()

    for x in readDict:
        for readId in readDict[x]:
            if readId not in readQualityDict:
                readQualityDict[readId] = []
                sizePercentDict[readId] = []
                isEdgeRead[readId] = []
                
            totalAlen = sum([aln.alen for aln in readDict[x][readId]])
            weightedQuality = sum([(aln.get_tag("mc")*aln.alen)/totalAlen for aln in readDict[x][readId]])
            percentAligned = sum([100*aln.qlen/aln.infer_read_length() for aln in readDict[x][readId]])

            start, end = 0, len(seqDict[x])
            startEdge, endEdge = EDGE_SIZE, max(0, end - EDGE_SIZE)


            startEdgeOverlap = sum([100*max(0, (min(startEdge, aln.reference_end) - max(start, aln.reference_start)))/EDGE_SIZE \
                               for aln in readDict[x][readId]])
            endEdgeOverlap = sum([100*max(0, (min(end, aln.reference_end) - max(endEdge, aln.reference_start)))/EDGE_SIZE \
                               for aln in readDict[x][readId]])
                        
            isEdgeOverlap = startEdgeOverlap > EDGE_THRESHOLD or endEdgeOverlap > EDGE_THRESHOLD

            readQualityDict[readId].append(weightedQuality)
            sizePercentDict[readId].append(percentAligned)
            isEdgeRead[readId].append(isEdgeOverlap)

    qualityFilterList = [readId for readId in readQualityDict if max(readQualityDict[readId]) < QUALITY_THRESHOLD ]
    sizeFilterList = [readId for readId in sizePercentDict if max(sizePercentDict[readId]) < SIZE_THRESHOLD ]
    edgeReadList = [readId for readId in isEdgeRead if sum(isEdgeRead[readId]) == len(isEdgeRead[readId]) ]

    filterSet = set(qualityFilterList).union(set(sizeFilterList) - set(edgeReadList))
    
    info = dict()
    for x in readDict:
        filteredReads = {readId : readDict[x][readId] for readId in readDict[x] if readId not in filterSet}
        filteredAlns = []
        for readId in filteredReads: filteredAlns.extend(filteredReads[readId])

        n50 = N50(filteredAlns)
        info["pacbio_align_n50_" + x] = 0 if n50 is None else n50
        info["pacbio_align_count_" + x] = len(filteredAlns)
        info["pacbio_read_count_" + x] = len(filteredReads)
        info["pacbio_read_count_" + x] = len(filteredReads)
        
        summaries = [AlignmentInfo(filteredReads[readId]) for readId in filteredReads]

        info["pacbio_align_bp_" + x] = sum([summary.query_aligned_length() for summary in summaries])
        info["pacbio_concordance_" + x] = 0 if info["pacbio_align_bp_" + x] == 0 else \
            sum([summary.weighted_concordance()*summary.query_aligned_length() for summary in summaries])/ \
                                            info["pacbio_align_bp_" + x]
            
        info["pacbio_match_" + x] = sum([summary.match_count() for summary in summaries])
        info["pacbio_indel_" + x] = sum([summary.indel_count() for summary in summaries])
        info["pacbio_mismatch_" + x] = sum([summary.mismatch_count() for summary in summaries])

    return info


def kmer_assessment(seqDict, outdir, param):
    
    k=21
    
    def assess(compress):
        kmers = dict()
        for x in seqDict:
            kmers[x] = kt.KmerSet(seqDict[x], k, compressed=compress)
        
        seqList = [x for x in seqDict]
        for i,x in enumerate(seqList):
            for y in seqList[i+1:]:
                pre = "compressed_" if compress else ""
                info[pre + "kmer_similarity_" + x + "_" + y] = kmers[x].jaccard_similarity(kmers[y])
    
        if len(seqList) == 3:
            kmerBreakdown = kt.draw_pie3(seqList[0], kmers[seqList[0]],
                                         seqList[1], kmers[seqList[1]],
                                         seqList[2], kmers[seqList[2]],
                                         outdir + pre + "kmer_breakdown")
        elif len(seqList) == 2:
            kmerBreakdown = kt.draw_pie2(seqList[0], kmers[seqList[0]],
                             seqList[1], kmers[seqList[1]],
                             outdir + pre + "kmer_breakdown")
            
        return kmerBreakdown
    
    
    info = dict()
    info.update(assess(False))
    info.update(assess(True))

    return info

def write_info(info, outdir, outputFileName, transpose=True):
    
    order = [key for key in info] ; order.sort()
    values = [str(info[x]) for x in order]
    
    out = open(outdir + outputFileName, 'w')
    out.write("\t".join(order) + "\n")
    out.write("\t".join(values))
    out.close()
    
    if transpose:
        out = open(outdir + "transpose_" + outputFileName, 'w')    
        out.write("\n".join([o + ":\t" + v for o,v in zip(order,values)]))    
        out.close()
        
def get_outdir(prefix, param, fork1, fork2):
    
    outdir = param.OUTPUT_DIR + "/" + prefix + "_" + fork1.qid + "_" + \
                              str(fork1.qpos) + "_" + \
                              str(fork2.qpos) + "/"
                              
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    return outdir
    
def get_regions(fork1, fork2, lengthData):
    
    qid, rid = fork1.qid, fork1.rid
    if qid != fork2.qid: 
        qRegion = None
    else:
        qp = [fork1.get_pos_norm(lengthData, q=True), 
              fork2.get_pos_norm(lengthData, q=True)]
        qRegion = rgn.SimpleRegion(qid, min(qp[0], qp[1]), max(qp[0], qp[1]))
        
    if rid != fork2.rid:
        rRegion = None
    else:
        rp = [fork1.get_pos_norm(lengthData, q=False), 
              fork2.get_pos_norm(lengthData, q=False)]
        rRegion = rgn.SimpleRegion(rid, min(rp[0], rp[1]), max(rp[0], rp[1]))

    return (qRegion, rRegion)

def analyze_block(blockPath, seqData, lengthData, param, skipIfDone=True):
    
    outputFileName = "info.txt"
    info = dict()

    outdir = get_outdir("block", param, blockPath[0], blockPath[-1])
    
    print("Dumping info to " + outdir)
    if skipIfDone and os.path.isfile(outdir + outputFileName):
        print("Analysis has already been done, skipping.")
        return
    
    qRegion, rRegion = get_regions(blockPath[0], blockPath[-1], lengthData) 
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    hSeq, source = io.path_to_sequence(blockPath, seqData)
    hSeq, source = "".join(hSeq), "".join(source)
    info["hybrid_r_composition"] = source.lower().count("r")/len(source)

    print("LENGTHS:", len(qSeq), len(rSeq), len(hSeq))
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    seqSummary = sequence_assessment(seqDict, param)
    info.update(seqSummary)
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    refGenomeSummary, hgRegion = ref_genome_assessment(seqDict, outdir, param, requireContinuous=[Q,R])

    info.update(refGenomeSummary)
    
    if hgRegion is None:
        return
        print(" WAIT! ")
        input()   
    else: print(hgRegion)

    hgSeq = tools.faidx_fetch(param.HG38, hgRegion) if hgRegion is not None else ""
    print("hybrid vs reference length", len(hSeq), len(hgSeq))

    seqDict = { HG: hgSeq }
    hgSeqSummary = sequence_assessment(seqDict, param)
    info.update(hgSeqSummary)

    seqDict = { Q: qSeq, R: rSeq, HG: hgSeq }
    kmerSummary = kmer_assessment(seqDict, outdir, param)
    info.update(kmerSummary)

    seqDict = { Q: qSeq, R: rSeq, H: hSeq, HG: hgSeq }
    pbSummary = pacbio_align_assessment(seqDict, rRegion, outdir, param)
    info.update(pbSummary)
    
    write_info(info, outdir, outputFileName)

    
def analyze_inter_block(prevBlockPath, blockPath, seqData, lengthData, param, skipIfDone=True):

    outputFileName = "info.txt"
    info = dict()

    outdir = get_outdir("inter", param, prevBlockPath[-1], blockPath[0])

    print("Dumping info to " + outdir)
    if skipIfDone and os.path.isfile(outdir + outputFileName):
        print("Analysis has already been done, skipping.")
        return

    qid, rid = prevBlockPath[-1].qid, prevBlockPath[-1].rid
    if qid != blockPath[0].qid or rid != blockPath[0].rid:
        print("Blockpath IDs don't match...")
    
    qRegion, rRegion = get_regions(prevBlockPath[-1], blockPath[0], lengthData) 
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    hSeq, source = rSeq, "r"*len(rSeq)
    info["hybrid_r_composition"] = source.lower().count("r")/len(source)

    print("LENGTHS:", len(qSeq), len(rSeq), len(hSeq))
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    seqSummary = sequence_assessment(seqDict, param)
    info.update(seqSummary)
    
    seqDict = { Q: qSeq, R: rSeq }
    info.update(refGenomeSummary)
    
    if hgRegion is None:
        return
        print(" WAIT! ")
        input()   
    else: print(hgRegion)

    hgSeq = tools.faidx_fetch(param.HG38, hgRegion) if hgRegion is not None else ""
    print("hybrid vs reference length", len(hSeq), len(hgSeq))
    
    seqDict = { HG: hgSeq }
    hgSeqSummary = sequence_assessment(seqDict, param)
    info.update(hgSeqSummary)

    seqDict = { R: rSeq, HG: hgSeq }
    kmerSummary = kmer_assessment(seqDict, outdir, param)
    info.update(kmerSummary)

    seqDict = { Q: qSeq, R: rSeq, HG: hgSeq }
    pbSummary = pacbio_align_assessment(seqDict, rRegion, outdir, param)
    info.update(pbSummary)

    write_info(info, outdir, outputFileName)
