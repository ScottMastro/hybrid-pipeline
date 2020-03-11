import sys
sys.path.append("./analysis")
import external_tools as tools
import paf_helper as paf
import variants as var
from regions import SimpleRegion, PhasedRegion, region_difference
import helper


def resolve_qvariants(alnA, alnB, region, seqData):
    
    rid = alnA["rid"]
    rseq = seqData[rid]
            
    variantSetQueryA = paf.parse_cs_string(alnA["cs"], rid, alnA["rstart"], rseq, noComplex=True)
    
    if alnB is None:
        variantSetQueryB = None
    else:
        variantSetQueryB = paf.parse_cs_string(alnB["cs"], alnB["rid"], alnB["rstart"], rseq, noComplex=True)
        #variantSetQuery = var.combine_as_genotypes(variantSetQueryA, variantSetQueryB, phased=True)
    
    result = PhasedRegion(region)
    result.add_variants(variantSetQueryA, variantSetQueryB)

    return result



'''
alignments = aligndf[(aligndf["rid"] == "tig00000446_pilon_pilon") & ((aligndf["qid"] == "1175") | (aligndf["qid"] == "1176"))]
alnA = alignments.iloc[0] ; alnB = alignments.iloc[1]
'''
def resolve_variants_2haplotigs(alnA, alnB, region, seqData, param, dummy=False):
    
    rid = alnA["rid"]
    rSeq = seqData[rid]
    prefix = helper.file_prefix(region, param)

    if dummy:
        print("Resolving with 2 supernova haplotigs:", region)
        return
    
    variantSetLS, variantSetPolish = resolve_variants_pacbio(region, seqData, param)
    
    variantSetQueryA = paf.parse_cs_string(alnA["cs"], alnA["rid"], alnA["rstart"], rSeq, noComplex=True)
    variantSetQueryB = paf.parse_cs_string(alnB["cs"], alnB["rid"], alnB["rstart"], rSeq, noComplex=True)          
    variantSetQuery = var.combine_as_genotypes(variantSetQueryA, variantSetQueryB, phased=True)

    calls, candidates = helper.combine_variantsets_2haplotigs(variantSetLS, variantSetPolish, variantSetQuery)
    return (calls, candidates)

def resolve_variants_1haplotig(alnA, region, seqData, param, dummy=False):
    '''
    Idea: call variants from haplotig, phase by hetsnps, polish with arrow
    
    '''
    
    rid = alnA["rid"]
    rSeq = seqData[rid]
    prefix = helper.file_prefix(region, param)
    
    if dummy:
        print("Resolving with 1 supernova haplotigs:", region)
        return

    #todo:
    '''
    subsetBam = tools.samtools_subset(param.REF_ALIGNED_READS, region, prefix)
       
    variantSet = paf.parse_cs_string(alnA["cs"], alnA["rid"], alnA["rstart"], seqData[rid], noComplex=True)
    variantSetHaploid2 = genotype_variants(rid, variantSet, subsetBam, param)
    variantSetHaploid2.filter_nonhet()
    variantSetHaploid2.filter_nonsnp()
    variantSetHaploid3 = phase_variantset(rid, variantSetHaploid2, subsetBam, param)
    hpA, hpB, unphased = phase_reads(rid, variantSetHaploid3, subsetBam, param)

    #todo: robustly determine if phasing worked

    variantSetPolish = polish_region(rid, rstart, rend, hpA, hpB, unphased, seqData, param, assignUnphased=True)



    #var.set_homopolymer_status(variantSetHaploid2, seqData, minRepeat=2)
    #variantSetHaploid2.filter_homopolymer()
    #variantSetHaploid2.filter_nonhet()
    #variantSetHaploid2.filter_nonsnp()
    #variantSetHaploid2.filter_by_tag("GQ", 100, gt=True)
    #variantSetHaploid2.unfilter_snps()

    
    #todo: robustly determine if phasing worked
    '''

def query_correct_1haplotig(region, aln, seqData, param, dummy=False):
    '''
    Idea: call variants from haplotig, phase by hetsnps, polish with arrow
    '''
    
    
    
    
    rid = alnA["rid"]
    rSeq = seqData[rid]
    prefix = helper.file_prefix(region, param)
    
    if dummy:
        print("Resolving with 1 supernova haplotigs:", region)
        return

    
# region = SimpleRegion("tig00000446_pilon_pilon", 84575, 84769, lengthData)
# simple phaseable gap:
# region = SimpleRegion("tig00000446_pilon_pilon", 445891, 449220, lengthData)
# strand bias insertion:
# region = SimpleRegion("tig00000446_pilon_pilon", 1723113, 1724572, lengthData)
# tough example low cov:
# region =  SimpleRegion("tig00000446_pilon_pilon", 2605537, 2606528, lengthData)   
# region = SimpleRegion("tig00000446_pilon_pilon", 3537856, 3538159)
def resolve_gap_method1(region, seqData, param, buffer=3000):
    '''
    '''
    ba=False

    # Isolate gap + buffer
    gapRegion = region.extend(buffer)
    prefix = helper.file_prefix(gapRegion, param)
    regionFa, regionBam = helper.isolate_region(gapRegion, seqData, param)

    # Polish with all reads
    polishedFa = iterative_polish(regionFa, regionBam, None, prefix +"_allreads", \
                                  niter=2, useUnphased=False, backAlign=ba)

    # Phase on polished sequence
    alignedBam = tools.align_pacbio(polishedFa, regionBam, prefix + "_polished")
    tools.longshot_genotype(alignedBam, polishedFa, prefix, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)

    # Check if phasing is good
    phasedReads = tools.samtools_fetch(bamA) + tools.samtools_fetch(bamB)
    unphasedReads = tools.samtools_fetch(bamUnphased)
    
    if len(phasedReads) < 0.1 * len(unphasedReads) :
        # Reads were not phased so polish with all of them again
        polishedFa = iterative_polish(polishedFa, regionBam, None, prefix +"_allreads_homozygous", \
                                      niter=3, useUnphased=False, backAlign=ba)
        #todo: call SVs?
        #tools.find_sv(alignedBam, polishedFa, prefix + "_polished")    

    else:
        # Reads were phased, polish each haplotype separately
        polishedFaA = iterative_polish(polishedFa, bamA, bamUnphased, \
                                   prefix +"_Areads", niter=3, useUnphased=True, backAlign=ba)
        polishedFaB = iterative_polish(polishedFa, bamB, bamUnphased, \
                                       prefix +"_Breads", niter=3, useUnphased=True, backAlign=ba)
    
       
        #todo: return something



# region = SimpleRegion("tig00000446_pilon_pilon", 84575, 84769, lengthData)
# simple phaseable gap:
# region = SimpleRegion("tig00000446_pilon_pilon", 445891, 449220, lengthData)
# strand bias insertion:
# region = SimpleRegion("tig00000446_pilon_pilon", 1723113, 1724572, lengthData)
#weird phase difference I can't figure out:
# region =  SimpleRegion("tig00000446_pilon_pilon", 2605537, 2606528, lengthData)   
# region = SimpleRegion("tig00000446_pilon_pilon", 3537856, 3538159)
# buffer = 15000

# buffer = 15000
def resolve_gap(region, seqData, param, buffer1=30000, buffer2=5000):
    '''
    Idea: 
        1) Isolate gap and polish with all reads.
        2) Realign to new consensus sequence.


    #todo: are reads from the same PS????????????

    vcfA = tools.pacbio_polish(bamA, reference, prefix + ".A", region, useMapFilter=useMapFilter)
    vcfB = tools.pacbio_polish(bamB, reference, prefix + ".B", region, useMapFilter=useMapFilter)

    return [vcfA, vcfB]

        1) Polish gap & realign reads.
        2) Attempt to phase
        3) Verify phasing
        4) Polish haplotypes         self.unitigMapSecondary = map_unitig_contig_custom(canuDir)
        (if they exist)
        5) Verify phased variants
        5b) Mark uncertain variants for 10x inspection??
        6) Output region as 1 or 2 fa files
        7) Reintegrate into hybrid???
    '''

    # 1) Isolate gap + buffer1, phase reads in this region.
    phaseRegion = region.extend(buffer1)
    prefix = helper.file_prefix(phaseRegion, param)
    #todo: what if extendedRegion does not include full buffer?
    tools.longshot_genotype(param.REF_ALIGNED_READS, param.REF_FA, prefix, region=phaseRegion, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    
    # 2) Subset reads to gap + buffer2 region.
    gapRegion = region.extend(buffer2)
    prefix = helper.file_prefix(gapRegion, param)
    bamA = tools.samtools_subset(bamA, gapRegion, prefix + ".A")
    bamB = tools.samtools_subset(bamB, gapRegion, prefix + ".B")
    bamUnphased = tools.samtools_subset(bamUnphased, gapRegion, prefix + ".unphased")
    bamA = tools.unalign_bam(bamA, prefix + ".uA")
    bamB = tools.unalign_bam(bamB, prefix + ".uB")
    bamUnphased = tools.unalign_bam(bamUnphased, prefix + ".uUnphased")

    # 3) Isolate region and polish with all reads.
    regionFa, regionBam = helper.isolate_region(gapRegion, seqData, param)
    gapConsensusFa = regionFa
    
    polishedFa = iterative_polish(gapConsensusFa, regionBam, None, \
                               prefix +"_allreads", niter=5, useUnphased=False, backAlign=True)

    alignedBam = tools.align_pacbio(polishedFa, regionBam, prefix)
    tools.longshot_genotype(alignedBam, polishedFa, prefix, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)


    # tools.pacbio_polish(regionBam, regionFa, prefix, outputFasta=True, useMapFilter=False)

    #todo?: polish again with all reads?

    # 4) Iteratively polish with both haplotypes.

    #todo: use uphased reads? when should I??
    
    polishedFaA = iterative_polish(polishedFa, bamA, bamUnphased, \
                                   prefix +"_Areads", niter=10, useUnphased=False, backAlign=True)
    polishedFaB = iterative_polish(polishedFa, bamB, bamUnphased, \
                                   prefix +"_Breads", niter=10, useUnphased=False, backAlign=True)

"""

def genotype_variants(tigId, variantSet, subsetBam, param):
    
    start, end = variantSet.get_range()
    prefix = file_prefix(tigId, start, end, param)

    vsVCF = prefix + ".variantset.vcf"
    vsVCFgz = variantSet.write_vcf(vsVCF, hetsOnly=False, compressed=True)
   
    whVCF = prefix + ".whatshap.gt.vcf"
    tools.whatshap_genotype(whVCF, param.REF_FA, vsVCFgz, subsetBam)                

    return var.vcf_to_variantset(whVCF)

def phase_reads(tigId, variantSet, subsetBam, param):

    start, end = variantSet.get_range()
    prefix = file_prefix(tigId, start, end, param)

    inVCF = prefix + ".variantset.vcf"
    inVCFgz = variantSet.write_vcf(inVCF, hetsOnly=False, compressed=True)
   
    outBAM = prefix + ".whatshap.bam"
    tools.whatshap_haplotag(outBAM, param.REF_FA, inVCFgz, subsetBam)                

    alignments = tools.samtools_fetch(outBAM, tigId, start, end)
    return sam.split_haplotype(alignments)

def phase_variantset(tigId, variantSet, subsetBam, param):

    start, end = variantSet.get_range()
    prefix = file_prefix(tigId, start, end, param)

    inVCF = prefix + ".variantset.vcf"
    variantSet.write_vcf(inVCF, hetsOnly=False)
    outVCF = prefix + ".whatshap.vcf"

    #todo use phase info in vcf?
    inputs = [subsetBam]
    tools.whatshap_phase(outVCF, param.REF_FA, inVCF, inputs, genotype=True, indels=True)
    
    #read in whatshap results
    variantSetWH = var.vcf_to_variantset(outVCF)
    return variantSetWH

def polish_region(region, reference, hpA, hpB, unphased, seqData, param, assignUnphased=False, useMapFilter=True):
    
    prefix = file_prefix(region, param)

    unphasedA = []
    unphasedB = []
    
    if assignUnphased:
        random.shuffle(unphased)
        half = math.ceil(len(unphased)/2)
        unphasedA = unphased[:half]
        unphasedB = unphased[half:]

    if len(hpA) > 1:
        print("Multiple phase blocks in this region detected.....")

    readsA, readsB = unphasedA, unphasedB
    
    for PS in hpA:
        readsA.extend(hpA[PS])
        readsB.extend(hpB[PS])
        
    readsA = sorted(readsA, key=lambda x: x.pos)
    readsB = sorted(readsB, key=lambda x: x.pos)

    bamA = prefix + ".A.bam"
    bamB = prefix + ".B.bam"
    tools.samtools_write(bamA, readsA, param.REF_ALIGNED_READS)
    tools.samtools_write(bamB, readsB, param.REF_ALIGNED_READS)
    vcfA = tools.pacbio_polish(bamA, reference, prefix + ".A", region, useMapFilter=useMapFilter)
    vcfB = tools.pacbio_polish(bamB, reference, prefix + ".B", region, useMapFilter=useMapFilter)


    '''
    for PS in hpA:
        bamA = prefix + "_PS" + str(PS) + ".A.bam"
        tools.samtools_write(bamA, hpA[PS], param.REF_ALIGNED_READS)
        vcfA = tools.pacbio_polish(bamA, param.REF_FA, prefix + ".A", tigId, start, end)
    
        bamB = prefix + "_PS" + str(PS) + ".B.bam"
        tools.samtools_write(bamB, hpB[PS], param.REF_ALIGNED_READS)
        vcfB = tools.pacbio_polish(bamB, param.REF_FA, prefix + ".B", tigId, start, end)
    '''
    
    rSeq = seqData[region.chrom]

    variantSetA = var.vcf_to_variantset(vcfA)
    variantSetB = var.vcf_to_variantset(vcfB)     
    variantSet = var.combine_as_genotypes(variantSetA, variantSetB, rSeq)
    
    return variantSet




def find_variation(alnA, alnB, seqData, param):
    #two query haplotigs paf alignments to reference

    rid = alnA["rid"]
    rseq = seqData[rid]
rid
    #identify variants supernova
    variantSetA = paf.parse_cs_string(alnA["cs"], rid, alnA["rstart"], rseq)
    variantSetB = paf.parse_cs_string(alnB["cs"], rid, alnB["rstart"], rseq)            
    variantSet = var.combine_as_genotypes(variantSetA, variantSetB, rseq)
    
    #identify variants canu
    start = alnA["rstart"]
    end = alnA["rend"]
    prefix = param.OUTPUT_DIR + "/" + "_".join([rid, str(start), str(end)])

    vsVCF = prefix + ".variantset.vcf"
    vsVCFgz = variantSet.write_vcf(vsVCF, hetsOnly=False, compressed=True)
   
    subsetBam = tools.samtools_subset(param.REF_ALIGNED_READS, param.OUTPUT_DIR,
                    chrom=rid, startPos=start, endPos=end)

    whBAM = prefix + ".whatshap.bam"
    tools.whatshap_haplotag(whBAM, param.REF_FA, vsVCFgz, subsetBam)                

    alignments = tools.samtools_fetch(whBAM, rid, start, end)
    hpA, hpB, unphased = sam.split_haplotype(alignments)

    if len(hpA) > 1:
        print("Multiple phase blocks in this region detected.....")

    for PS in hpA:
        bamA = prefix + "_PS" + str(PS) + ".A.bam"
        tools.samtools_write(bamA, hpA[PS], whBAM)

        vcfA = tools.pacbio_polish(bamA, param.REF_FA, prefix + ".A", rid, start, end)
    

        bamB = prefix + "_PS" + str(PS) + ".B.bam"
        tools.samtools_write(bamB, hpB[PS], whBAM)

    
    
    '''
    #isolate reads in this region
    rstart = min(alnA["rstart"], alnB["rstart"])
    rend = max(alnA["rend"], alnB["rend"])
 
    subsetBam = tools.samtools_subset(param.REF_ALIGNED_READS, param.OUTPUT_DIR,
                    chrom=rid, startPos=rstart, endPos=rend)

    #phase variants / reads with whatshap
    whIn = param.OUTPUT_DIR + "/whatshap.in.vcf"
    whOut = param.OUTPUT_DIR + "/whatshap.out.vcf"
    variantSet.write_vcf(whIn, hetsOnly=False)
    
    inputs = [whIn, subsetBam]
    tools.whatshap_phase(whOut, param.REF_FA, whIn, inputs, genotype=True, indels=True)

    rstart = min(alnA["rstart"], alnB["rstart"])
    rend = max(alnA["rend"], alnB["rend"])

    #read in whatshap results
    variantSetWH = var.vcf_to_variantset(whOut)
    '''
    
    region = Region(rid, alnA, alnB=alnB)
    region.add_supernova_variantset(variantSet)
    #region.add_whatshap_variantset(variantSetWH)
    return [region]



    
def local_assemble(tigId, start, end, readsA, readsB, param):

    nameA =  "_".join([tigId, str(start), str(end), "A"])
    outdirA = param.OUTPUT_DIR + "/" + nameA
    
    readsAfasta = tools.bam2fasta(readsA)
    tools.canu_correct(readsAfasta, nameA, outdirA, end-start)
    
    nameB =  "_".join([tigId, str(start), str(end), "B"])
    outdirB = param.OUTPUT_DIR + "/" + nameB

    readsBfasta = tools.bam2fasta(readsB)
    tools.canu_correct(readsBfasta, nameB, outdirB, end-start)

    name =  "_".join([tigId, str(start), str(end)])
    outdir = param.OUTPUT_DIR + "/" + name
    
    correctedReadsA = outdirA + "/" + nameA + ".correctedReads.fasta.gz"
    correctedReadsB = outdirB + "/" + nameB + ".correctedReads.fasta.gz"
    tools.canu_assemble([correctedReadsA, correctedReadsB], name, outdir, end-start)


    #todo: unphased?
    
    


    return (variantSetSVA, variantSetSVB, variantSetSV)


def match_quality(qid, aligndf):
    alignments = aligndf[aligndf["qid"] == qid]
    rids = set(list(alignments["rid"]))
    
    qlen = list(alignments["qlen"])[0]
    
    for rid in rids:

        rAlignments = alignments[alignments["rid"] == rid]
        rlen = list(rAlignments["rlen"])[0]
        aln = sum(list(rAlignments["alnlen"]))
        divergance = np.sum(np.array(rAlignments["dv"]) * np.array(rAlignments["alnlen"])) / aln
        print(rid, "\t",
              str(round(100*aln/rlen,1)) + " of ref\t", 
              str(round(100*aln/qlen,2)) + " of query\t",
              "dv: " + str(round(divergance,4)))
    
    
def fill_gap(region, seqData, param, buffer=15000):
    '''
    Idea: 
        1) Isolate gap and polish with all reads.
        2) Realign to new consensus sequence.


    #todo: are reads from the same PS????????????

    vcfA = tools.pacbio_polish(bamA, reference, prefix + ".A", region, useMapFilter=useMapFilter)
    vcfB = tools.pacbio_polish(bamB, reference, prefix + ".B", region, useMapFilter=useMapFilter)

    return [vcfA, vcfB]

        1) Polish gap & realign reads.
        2) Attempt to phase
        3) Verify phasing
        4) Polish haplotypes (if they exist)
        5) Verify phased variants
        5b) Mark uncertain variants for 10x inspection??
        6) Output region as 1 or 2 fa files
        7) Reintegrate into hybrid???
    '''


    # 1) Isolate gap and polish with all reads.
    extendedRegion = region.extend(buffer)
    buff=buffer-2000
    #todo: what if extendedRegion does not include full buffer?
    regionFa, regionBam = helper.isolate_region(extendedRegion, seqData, param)

    prefix = helper.file_prefix(extendedRegion, param)
    unalignedBam = tools.unalign_bam(regionBam, prefix)

    fid = "_".join([extendedRegion.chrom, str(extendedRegion.start), str(extendedRegion.end)])
    gapRegion = SimpleRegion(fid, buff, len(extendedRegion) - buff)
    
    gapConsensusFa = tools.pacbio_polish(regionBam, regionFa, prefix, outputFasta=True, useMapFilter=False, region=gapRegion)
    gapLen = helper.get_fasta_len(gapConsensusFa)

    # only polished region is output so we reappend the unpolished flanking region
    polishedFa = helper.append_flanking_region(regionFa, gapConsensusFa, gapRegion, prefix + ".polished")
    gapRegion = SimpleRegion(fid, buff, buff + gapLen)

    # 2) Realign to new consensus sequence.
    realignedBam = tools.align_pacbio(polishedFa, unalignedBam, prefix + ".realigned")    
    
    # 3) Phase reads if possible.
    longshotVCF = tools.longshot_genotype(realignedBam, polishedFa, prefix, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    bamA = tools.unalign_bam(bamA, prefix + ".A")
    bamB = tools.unalign_bam(bamB, prefix + ".B")
    bamUnphased = tools.unalign_bam(bamUnphased, prefix + ".unphased")

    def repolish(fa, bam, bamUnphased, prefix):
        global gapRegion
        
        '''
        reads = tools.samtools_fetch(bam)
        
        if bamUnphased is not None:
            unphased = tools.samtools_fetch(bamUnphased)
            random.shuffle(unphased)
            half = math.ceil(len(unphased)/2)
            reads = reads + unphased[:half]
            reads = sorted(reads, key=lambda x: x.pos)
            
        readsBam = tools.samtools_write(reads, prefix, bam)
        '''

        alignedBam = tools.align_pacbio(fa, bam, prefix)
        gapFa = tools.pacbio_polish(alignedBam, fa, prefix, outputFasta=True, useMapFilter=False, region=gapRegion)
        gapLen = helper.get_fasta_len(gapFa)
        gapRegion = SimpleRegion(fid, buff, buff + gapLen)
        return helper.append_flanking_region(fa, gapFa, gapRegion, prefix + ".polished")
        
    iterations = 0
    
    faA, faB = polishedFa, polishedFa
    
    while iterations < 20:
        
        unphased = bamUnphased
        if iterations > 5: unphased = None
        faA = repolish(faA, bamA, unphased, prefix + "_A_" + str(iterations))
        faB = repolish(faB, bamB, unphased, prefix + "_B_" + str(iterations))
        
        alignedBam = tools.align_pacbio(regionFa, faA, param.OUTPUT_DIR + "/A_" + str(iterations))
        alignedBam = tools.align_pacbio(regionFa, faB, param.OUTPUT_DIR + "/B_" + str(iterations))

        
        iterations +=1
       
"""