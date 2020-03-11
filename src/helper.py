import external_tools as tools
import structures.region as rgn
import variants as var
import os.path
import re
import os
import shutil

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

def copy_file(file, newFile):
    shutil.copyfile(file, newFile)
    return newFile

def move_file(file, newFile):
    shutil.move(file, newFile)
    return newFile

def rename_file(file, newFile):
    os.rename(file, newFile)
    return newFile

def make_dir(file, directory):
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

def delete_file(file, deleteDirTree=False):
    if os.path.isfile(file):
        os.remove(file)
    if os.path.isfile(file + ".bai"):
        os.remove(file + ".bai")
    if os.path.isfile(file + ".fai"):
        os.remove(file + ".fai")
    if os.path.isfile(file + ".tbi"):
        os.remove(file + ".tbi")

    if os.path.isdir(file):
        try:
            os.rmdir(file)
        except:
            if deleteDirTree: shutil.rmtree(file)


def file_prefix(region, param):
    chromName = re.sub('[^a-zA-Z0-9 \n\.]', '-', region.chrom)
    
    return param.OUTPUT_DIR + "/" + "_".join([chromName, str(region.start), str(region.end)])

def contig_fasta(tigName, seqData, param):
    fastaPrefix = param.OUTPUT_DIR + "/" + tigName + "_full"
    fastaName = fastaPrefix + ".fasta"
    if not os.path.isfile(fastaName):
        faDict = { tigName : seqData[tigName] }
        return tools.dict2fasta(faDict, fastaPrefix)
    return fastaName

def contig_region(tigName, lengthData, param):
    return rgn.SimpleRegion(tigName, 0, lengthData[tigName]-1, lengthData)

def get_fasta_len(faFile, fid=None):   
    faDict = tools.fasta2dict(faFile)
    
    if fid is None:
        fid = list(faDict.keys())[0]
    
    return len(faDict[fid])

def get_fasta_seq(faFile, region=None):   
    faDict = tools.fasta2dict(faFile)
    
    if region is None:
        fid = list(faDict.keys())[0]
        return faDict[fid]

    fid = region.chrom
    return faDict[fid][region.start:region.end]

def get_fasta_id(faFile, region=None):   
    faDict = tools.fasta2dict(faFile)
    fid = list(faDict.keys())[0]
    return fid


def fetch_read_ids(bamFile):
    alignments = tools.samtools_fetch(bamFile)

    # remove redundant reads:
    readDict = {read.qname : read for read in alignments}
    uniqueIds = [k for k in readDict]
    return uniqueIds

def combine_reads(bamList, prefix):
    bam = []
    for bamFile in bamList:
        reads = tools.samtools_fetch(bamFile)        
        bam.extend(reads)
        
    return tools.samtools_write(bam, prefix + "_combined", bamList[0])
    

def append_flanking_region(fullFaFile, subFaFile, subRegion, prefix, newFaId=None):
    
    fullFaDict = tools.fasta2dict(fullFaFile)
    subFaDict = tools.fasta2dict(subFaFile)

    #assumes exactly one sequence in fasta files:
    fullId = list(fullFaDict.keys())[0]
    subId = list(subFaDict.keys())[0]

    newSeq = fullFaDict[fullId][:subRegion.start] + subFaDict[subId] + fullFaDict[fullId][subRegion.end+1:]
    
    if newFaId is None:
        newFaId = fullId
        
    newFaDict = {newFaId : newSeq}
    return tools.dict2fasta(newFaDict, prefix, toUpper=True)

def get_average_divergance(bamFile, tag="mc", normalize=False, region=None):    
    reads = tools.samtools_fetch(bamFile, region=region)
    mc = [read.get_tag(tag) for read in reads]
    rq = [1 for read in reads]
    if normalize:
        rq = [100*read.get_tag("rq") for read in reads]
    norm = [m/q for m,q in zip(mc,rq)]
    return sum(norm)/len(norm)

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
    listBlockSummary = dict()
    
    for op in ops: 
        listSummary[op] = 0
        listBlockSummary[op] = 0
        
    for alignment in alignmentList:
    
        if alignment is None:
            continue
    
        x = alignment.get_cigar_stats()
        
        for i,op in enumerate(ops):
            listSummary[op] += x[0][i]
            listBlockSummary[op] += x[1][i]
    
    def calc(s):
        s["total"] = sum([s[key] for key in s])
        s["indel"] = s["I"] + s["D"] 
        s["clip"]  = s["S"] + s["H"] 
        s["match"] = s["M"] + s["="] 
        return s
    
    return calc(listBlockSummary) if block else calc(listSummary)


def phase_sv_reads(regionFa1, readsBam1, regionFa2, readsBam2, prefix, param, readsBamUnphased=None):  
    reads1, reads2, reads3 = tools.samtools_fetch(readsBam1), tools.samtools_fetch(readsBam2), []

    if readsBamUnphased is not None:
        reads3 = tools.samtools_fetch(readsBamUnphased)

    allReads = reads1 + reads2 + reads3
    readIds = set([read.qname for read in allReads])
    readIdMap = {read.qname : read for read in allReads}

    combinedBam = tools.samtools_write(allReads, prefix + "_TEMP_combined", param.REF_ALIGNED_READS, makeUnique=True)
    uBam = tools.unalign_bam(combinedBam, prefix)
    alignedR1Bam = tools.align_pacbio(regionFa1, uBam, prefix + "_TEMP_region1")
    alignedR2Bam = tools.align_pacbio(regionFa2, uBam, prefix + "_TEMP_region2")
    readDict1, readDict2 = dict(), dict()
    
    for read in tools.samtools_fetch(alignedR1Bam):
        if read.qname in readDict1: readDict1[read.qname].append(read)
        else: readDict1[read.qname] = [read]
    for read in tools.samtools_fetch(alignedR2Bam):
        if read.qname in readDict2: readDict2[read.qname].append(read)
        else: readDict2[read.qname] = [read]

    readAssignment1 = set()
    readAssignment2 = set()
    unassignedReads = set()

    for readId in readIds:
        aln1 = readDict1[readId] if readId in readDict1 else None
        aln2 = readDict2[readId] if readId in readDict2 else None
        summary1 = cigar_summary(aln1)
        summary2 = cigar_summary(aln2)
        
        concordanceRatio1 = summary1["match"]/(summary1["X"] + summary1["clip"] + summary1["indel"] + 1)
        concordanceRatio2 = summary2["match"]/(summary2["X"] + summary2["clip"] + summary2["indel"] + 1)

        concordDiff = (concordanceRatio1 + 1e-8) / (concordanceRatio2 + 1e-8)

        #todo: these numbers might need to be adjusted...
        if concordDiff > 1.429:
            readAssignment1.add(readId)
        elif concordDiff < 0.7:
            readAssignment2.add(readId)
        else:
            unassignedReads.add(readId)
            
        '''
        print("=============================")
        print(cigar_summary(aln1))
        print(cigar_summary(aln2))
        print(readId)
        print(round(concordDiff,4), round(concordanceRatio1,4),round(concordanceRatio2, 4))
        '''
        
    #here we phase using all reads aligned to region1, might be better
    #to sometimes exclude the reads assigned to region2?
    tools.longshot_genotype(alignedR1Bam, regionFa1, prefix, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)

    readsA = tools.samtools_fetch(bamA)
    readsB = tools.samtools_fetch(bamB)
    readsUnphased = tools.samtools_fetch(bamUnphased)
    
    readIdsA = set([read.qname for read in readsA])
    readIdsB = set([read.qname for read in readsB])

    evidenceA = len(readAssignment1.intersection(readIdsA)) - len(readAssignment1.intersection(readIdsB))
    evidenceB = len(readAssignment2.intersection(readIdsB)) - len(readAssignment2.intersection(readIdsA))

    if (evidenceA + evidenceB) < 0:
        readIdsA, readIdsB = readIdsB, readIdsA
        
    readIdsA = (readIdsA - readAssignment2).union(readAssignment1)
    readIdsB = (readIdsB - readAssignment1).union(readAssignment2)
    readIdsUnphased = readIds - readIdsA - readIdsB
    
    readsA = [readIdMap[readId] for readId in readIdsA]
    readsB = [readIdMap[readId] for readId in readIdsB]
    readsUnphased = [readIdMap[readId] for readId in readIdsUnphased]
    bamA = tools.samtools_write(readsA, prefix + ".A", param.REF_ALIGNED_READS, makeUnique=False)
    bamB = tools.samtools_write(readsB, prefix + ".B", param.REF_ALIGNED_READS, makeUnique=False)
    bamUnphased = tools.samtools_write(readsUnphased, prefix + ".unphased", param.REF_ALIGNED_READS, makeUnique=False)
    bamA = tools.unalign_bam(bamA, prefix + ".A")
    bamB = tools.unalign_bam(bamB, prefix + ".B")
    bamUnphased = tools.unalign_bam(bamUnphased, prefix + ".unphased")
    
    return (bamA, bamB, bamUnphased)
    
def isolate_region_fasta(region, seqData, param):

    prefix = file_prefix(region, param)
    
    seqDict = dict()
    fid = "_".join([region.chrom, str(region.start), str(region.end)])
    seqDict[fid] = seqData[region.chrom][region.start:region.end]

    isolationFa = tools.dict2fasta(seqDict, prefix)
    tools.samtools_faidx(isolationFa)

    return isolationFa


def isolate_region(seqData, bamFile, region, prefix, unaligned=False):
    
    seqDict = dict()
    fid = "_".join([region.chrom, str(region.start), str(region.end)])
    seqDict[fid] = seqData[region.chrom][region.start:region.end]

    isolationFa = tools.dict2fasta(seqDict, prefix)
    tools.samtools_faidx(isolationFa)

    subsetBam = tools.samtools_subset(bamFile, region, prefix + "_TEMP_")
    
    if unaligned:
        unalignedBam = tools.unalign_bam(subsetBam, prefix)
        delete_file(subsetBam)
        return [isolationFa, unalignedBam]
    
    realignedBam = tools.align_pacbio(isolationFa, subsetBam, prefix)
    delete_file(subsetBam)
    return [isolationFa, realignedBam]

def isolate_phased_region(region, seqData, bamA, bamB, bamUnphased, param, unaligned=False):

    prefix = file_prefix(region, param)
    
    seqDict = dict()
    fid = "_".join([region.chrom, str(region.start), str(region.end)])
    seqDict[fid] = seqData[region.chrom][region.start:region.end]

    isolationFa = tools.dict2fasta(seqDict, prefix)
    tools.samtools_faidx(isolationFa)

    subsetBamA = tools.samtools_subset(bamA, region, prefix + "_A_TEMP")
    subsetBamB = tools.samtools_subset(bamB, region, prefix + "_B_TEMP")
    subsetBamU = tools.samtools_subset(bamUnphased, region, prefix + "_U_TEMP")

    if unaligned:
        unalignedBamA = tools.unalign_bam(subsetBamA, prefix + "_A")
        unalignedBamB = tools.unalign_bam(subsetBamB, prefix + "_B")
        unalignedBamU = tools.unalign_bam(subsetBamU, prefix + "_U")
        
        delete_file(subsetBamA) ; delete_file(subsetBamB) ; delete_file(subsetBamU)
        
        return [isolationFa, unalignedBamA, unalignedBamB, unalignedBamU]
    
    realignedBamA = tools.align_pacbio(isolationFa, subsetBamA, prefix + "_A")
    realignedBamB = tools.align_pacbio(isolationFa, subsetBamB, prefix + "_B")
    realignedBamU = tools.align_pacbio(isolationFa, subsetBamU, prefix + "_U")
    
    delete_file(subsetBamA) ; delete_file(subsetBamB) ; delete_file(subsetBamU)

    return [isolationFa, realignedBamA, realignedBamB, realignedBamU]


LONGSHOT="LS"
ARROW="AR"
QUERYTIG="QT"

def extract_confident_calls(candidates):
    calls = []
    filteredCandidates = []
    for candidate in candidates:
        if candidate.get_confidence() == var.Confidence.CALL or \
           candidate.get_confidence() == var.Confidence.CONSENSUS:
               calls.append(candidate.to_variant_call())
        else:
            filteredCandidates.append(candidate)
            
    return (calls, filteredCandidates)

'''

variantSetLS = var.vcf_to_variantcallset(longshotVCF, region)

variantSetPolishA = var.vcf_to_variantcallset(vcfA, region)
variantSetPolishB = var.vcf_to_variantcallset(vcfB, region) 
variantSetPolish = var.combine_as_genotypes(variantSetPolishA, variantSetPolishB, phased=True)

#variantSetPacbio = helper.combine_variantsets_pacbio(variantSetLS, variantSetPolish)

variantSetQueryA = paf.parse_cs_string(alnA["cs"], alnA["rid"], alnA["rstart"], rSeq, noComplex=True)
variantSetQueryB = paf.parse_cs_string(alnB["cs"], alnB["rid"], alnB["rstart"], rSeq, noComplex=True)          
variantSetQuery = var.combine_as_genotypes(variantSetQueryA, variantSetQueryB, phased=True)

pos=19030

'''
def combine_variantsets_2haplotigs(variantSetLS, variantSetPolish, variantSetQuery):

    candidates = []
    
    #todo: normalize phase!
        
    for pos in variantSetLS.get_positions():
        
        snp = variantSetLS.get(pos) #.pop(pos)
        polishedVariant = variantSetPolish.get(pos) #.pop(pos)
        queryVariant = variantSetQuery.get(pos) #.pop(pos)
        
        candidate = var.VariantCandidate()
        candidate.add_call(LONGSHOT, snp)
        candidate.add_call(ARROW, polishedVariant)
        candidate.add_call(QUERYTIG, queryVariant)
        candidates.append(candidate)

        pacbioResult = var.compare(snp, polishedVariant)
        queryResult = var.compare(snp, queryVariant)

        print(pos, pacbioResult.name, queryResult.name)

        #genotype is agreed
        if var.agree_genotype(pacbioResult) and var.agree_genotype(queryResult):
            #consensus
            if (pacbioResult == var.Result.AGREE or pacbioResult == var.Result.AGREE_UNPHASED) and \
                queryResult == var.Result.AGREE:
                candidate.set_confidence(var.Confidence.CONSENSUS)
                continue
            
            #phasing disagreements / uncalled phase
            candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
            candidate.add_issue(var.Issues.QUESTIONABLE_PHASE)
            continue

        #some level of genotype disagreement
        if var.disagree_genotype(pacbioResult) or var.disagree_genotype(queryResult):
            candidate.set_confidence(var.Confidence.CANDIDATE)
            candidate.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)
            continue
            
        #missing a call
        if queryResult == var.Result.MISSING:
            if pacbioResult == var.Result.AGREE:
                candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
                candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)
            elif pacbioResult == var.Result.AGREE_UNPHASED or pacbioResult == var.Result.INCONSISTENT_PHASE:
                candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
                candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)
                candidate.add_issue(var.Issues.QUESTIONABLE_PHASE)
            elif pacbioResult == var.Result.MISSING:
                candidate.set_confidence(var.Confidence.UNLIKELY_CANDIDATE)
                candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)
            continue
        
        #missing a call
        if pacbioResult == var.Result.MISSING:
            if queryResult == var.Result.AGREE:
                candidate.set_confidence(var.Confidence.CALL)
            elif queryResult == var.Result.AGREE_UNPHASED or queryResult == var.Result.INCONSISTENT_PHASE:
                candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
                candidate.add_issue(var.Issues.QUESTIONABLE_PHASE)
            continue
            
        print("unhandled case 1")


    for pos in variantSetPolish.get_positions():
        
        polishedVariant = variantSetPolish.get(pos) #.pop(pos)
        queryVariant = variantSetQuery.get(pos) #.pop(pos)
        
        candidate = var.VariantCandidate()
        candidate.add_call(LONGSHOT, None)
        candidate.add_call(ARROW, polishedVariant)
        candidate.add_call(QUERYTIG, queryVariant)
        candidates.append(candidate)

        result = var.compare(polishedVariant, queryVariant)

        print(pos, result.name)

        #not called by longshot!
        if polishedVariant.is_snp():
            candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)

        if result == var.Result.AGREE:
            if polishedVariant.is_snp():
                candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
            else:
                candidate.set_confidence(var.Confidence.CONSENSUS)
            continue
        
        if result == var.Result.AGREE_UNPHASED or pacbioResult == var.Result.INCONSISTENT_PHASE:
            candidate.set_confidence(var.Confidence.LIKELY_CANDIDATE)
            candidate.add_issue(var.Issues.QUESTIONABLE_PHASE)
            continue
        if result == var.Result.MISSING:
            if polishedVariant.is_snp():
                candidate.set_confidence(var.Confidence.UNLIKELY_CANDIDATE)
            else:
                candidate.set_confidence(var.Confidence.CANDIDATE)
                candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)
            continue
        if var.disagree_genotype(result):
            candidate.set_confidence(var.Confidence.CANDIDATE)
            candidate.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)
            continue

        print("unhandled case 2")

        for pos in variantSetQuery.get_positions():
        
            queryVariant = variantSetPolish.get(pos) #.pop(pos)
            
            candidate = var.VariantCandidate()
            candidate.add_call(LONGSHOT, None)
            candidate.add_call(ARROW, None)
            candidate.add_call(QUERYTIG, queryVariant)
            candidates.append(candidate)
    
            candidate.set_confidence(var.Confidence.CANDIDATE)
            candidate.add_issue(var.Issues.INCOMPLETE_SUPPORT)

    return extract_confident_calls(candidates)

        
        
'''

variantSetLS = var.vcf_to_variantcallset(longshotVCF, region)

variantSetPolishA = var.vcf_to_variantcallset(vcfA, region)
variantSetPolishB = var.vcf_to_variantcallset(vcfB, region) 
variantSetPolish = var.combine_as_genotypes(variantSetPolishA, variantSetPolishB, phased=True)

variantSetPacbio = helper.combine_variantsets_pacbio(variantSetLS, variantSetPolish)

variantSetQueryA = paf.parse_cs_string(alnA["cs"], alnA["rid"], alnA["rstart"], rSeq, noComplex=True)
variantSetQueryB = paf.parse_cs_string(alnB["cs"], alnB["rid"], alnB["rstart"], rSeq, noComplex=True)          
variantSetQuery = var.combine_as_genotypes(variantSetQueryA, variantSetQueryB, phased=True)

'''
def combine_variantsets_pacbio_query(variantSetPacbio, variantSetQuery, bothHaplotypes=True):

    combinedSet = var.VariantCallSet()

    confidence_upgrade = { 
        var.Confidence.UNKNOWN            : var.Confidence.UNKNOWN,
        var.Confidence.CONSENSUS          : var.Confidence.CONSENSUS,
        var.Confidence.LIKELY_CANDIDATE   : var.Confidence.CONSENSUS,
        var.Confidence.CANDIDATE          : var.Confidence.CONSENSUS,
        var.Confidence.UNLIKELY_CANDIDATE : var.Confidence.CALL,
        var.Confidence.SPURIOUS           : var.Confidence.CANDIDATE,
        var.Confidence.CALL               : var.Confidence.CONSENSUS
        }

    confidence_downgrade = { 
        var.Confidence.UNKNOWN            : var.Confidence.UNKNOWN,
        var.Confidence.CONSENSUS          : var.Confidence.LIKELY_CANDIDATE,
        var.Confidence.LIKELY_CANDIDATE   : var.Confidence.CANDIDATE,
        var.Confidence.CANDIDATE          : var.Confidence.UNLIKELY_CANDIDATE,
        var.Confidence.UNLIKELY_CANDIDATE : var.Confidence.UNLIKELY_CANDIDATE,
        var.Confidence.SPURIOUS           : var.Confidence.SPURIOUS,
        var.Confidence.CALL               : var.Confidence.CANDIDATE
        }
    
    for pos in variantSetPacbio.get_positions():
        variant = variantSetPacbio.get(pos) #.pop(pos)
        otherVariant = variantSetQuery.get(pos) #.pop(pos)

        result = var.compare(variant, otherVariant)
        print(result.name)
        
        if result == var.Result.AGREE:
            variant.set_confidence(confidence_upgrade[variant.get_confidence()])
            
        elif result == var.Result.AGREE_UNPHASED or result == var.Result.INCONSISTENT_PHASE:
            newCall.set_confidence(var.Confidence.LIKELY_CANDIDATE)
            newCall.add_issue(var.Issues.QUESTIONABLE_PHASE)
        
        elif result == var.Result.DIFFERENT_ALLELE or result == var.Result.POTENTIAL_BIALLELIC:
            #todo: need to see examples of this
            #tig00000446_pilon_pilon:20,690-20,730
            print("TODO: RESOLVE DIFFERENT ALLELE ISSUE")
            print(snp)
            print(otherVariant)

            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)

        elif result == var.Result.UNCALLED:
            #this should not occur in this loop
            print("WARNING: Variant reported as UNCALLED when it should exist")
            continue
        
        elif result == var.Result.MISSING:
            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.INCOMPLETE_SUPPORT)

        elif result == var.Result.POTENTIAL_HET or result == var.Result.POTENTIAL_HOMALT:
            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)

        combinedSet.add(newCall)

    for pos in variantSetPolish.get_positions():
        variant = variantSetPolish.pop(pos)
        alleleA, alleleB = variant.get_alleleA(clean=True), variant.get_alleleB(clean=True)
        newCall = var.VariantCall(alleleA, alleleB, phased=True)
        if variant.is_snp():
            newCall.set_confidence(var.Confidence.UNLIKELY_CANDIDATE)
        else:
            newCall.set_confidence(var.Confidence.CANDIDATE)

        #todo: assess arrow confidence in genotype and # of phased reads at every pos
        newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)
        newCall.add_issue(var.Issues.QUESTIONABLE_PHASE)
        combinedSet.add(newCall)

    return combinedSet
        




'''
def combine_variantsets(variantSetLS, variantSetPolish, variantSetQuery, queryHaplotypes=2):

    combinedSet = var.VariantCallSet()
    
    for pos in variantSetLS.get_positions():
        
        snp = variantSetLS.get(pos) #.pop(pos)
        polishedVariant = variantSetPolish.get(pos) #.pop(pos)
        queryVariant = variantSetQuery.get(pos) #.pop(pos)
        
        result = var.compare(snp, polishedVariant)
        
        alleleA, alleleB = snp.get_alleleA(clean=True), snp.get_alleleB(clean=True)
        phased, ps = snp.is_phased(), snp.get_phase_set()
        newCall = var.VariantCall(alleleA, alleleB, phaseSet=ps, phased=phased)

        if result == var.Result.AGREE:
            newCall.set_confidence(var.Confidence.CALL)
            
        elif result == var.Result.AGREE_UNPHASED or result == var.Result.INCONSISTENT_PHASE:
            newCall.set_confidence(var.Confidence.LIKELY_CANDIDATE)
            newCall.add_issue(var.Issues.QUESTIONABLE_PHASE)
        
        elif result == var.Result.DIFFERENT_ALLELE or result == var.Result.POTENTIAL_BIALLELIC:
            #todo: need to see examples of this
            #tig00000446_pilon_pilon:20,690-20,730
            print("TODO: RESOLVE DIFFERENT ALLELE ISSUE")
            print(snp)
            print(polishedVariant)

            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)

        elif result == var.Result.UNCALLED:
            #this should not occur in this loop
            print("WARNING: Variant reported as UNCALLED when it should exist")
            continue
        
        elif result == var.Result.MISSING:
            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.INCOMPLETE_SUPPORT)

        elif result == var.Result.POTENTIAL_HET or result == var.Result.POTENTIAL_HOMALT:
            newCall.set_confidence(var.Confidence.CANDIDATE)
            newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)

        combinedSet.add(newCall)

    for pos in variantSetPolish.get_positions():
        variant = variantSetPolish.pop(pos)
        alleleA, alleleB = variant.get_alleleA(clean=True), variant.get_alleleB(clean=True)
        newCall = var.VariantCall(alleleA, alleleB, phased=True)
        if variant.is_snp():
            newCall.set_confidence(var.Confidence.UNLIKELY_CANDIDATE)
        else:
            newCall.set_confidence(var.Confidence.CANDIDATE)

        #todo: assess arrow confidence in genotype and # of phased reads at every pos
        newCall.add_issue(var.Issues.GENOTYPE_UNCERTAINTY)
        newCall.add_issue(var.Issues.QUESTIONABLE_PHASE)
        combinedSet.add(newCall)

    return combinedSet
'''


