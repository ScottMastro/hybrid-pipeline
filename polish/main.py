import sys
sys.path.append("./analysis")
import log.log as logger
import parameters
import file_handler as io
import external_tools as tools
import supernova_helper as supernova
import canu_helper as canu
import paf_helper as paf
import variants as var
import regions as rgn
import contig_connector as contigger
import sam_helper as sam
import resolve_region as resolver
import helper
import polish_helper as polisher
import assessment as assesser

def find_structural_variation(region, param):
    svVCF = "/media/scott/HDD/sickkids/CF062/pacbio_to_canu.pbsv.vcf"
    variantSetSV = var.vcf_to_variantcallset(svVCF, region)

    subsetBam = tools.samtools_subset(param.REF_ALIGNED_READS, param.OUTPUT_DIR,
                    chrom=tigId, startPos=start, endPos=end)

    prefix = param.OUTPUT_DIR + "/" + "_".join([tigId, str(start), str(end)])

    vsVCF = prefix + ".variantset.vcf"
    vsVCFgz = variantSet.write_vcf(vsVCF, hetsOnly=False, compressed=True)

    whBAM = prefix + ".whatshap.bam"
    tools.whatshap_haplotag(whBAM, param.REF_FA, vsVCFgz, subsetBam)                

    alignments = tools.samtools_fetch(whBAM, tigId, start, end)
    hpA, hpB, unphased = sam.split_haplotype(alignments)

    if len(hpA) > 1:
        print("Multiple phase blocks in this region detected.....")

    for PS in hpA:
        bamA = prefix + "_PS" + str(PS) + ".A.bam"
        tools.samtools_write(bamA, hpA[PS], whBAM)

        svPrefix = prefix + ".A.pbsv"
        svFile = tools.find_sv(bamA, param.REF_FA, svPrefix)
        variantSetSVA = var.vcf_to_variantset(svFile)
                
        bamB = prefix + "_PS" + str(PS) + ".B.bam"
        tools.samtools_write(bamB, hpB[PS], whBAM)

        svPrefix = prefix + ".B.pbsv"
        svFile = tools.find_sv(bamB, param.REF_FA, svPrefix)
        variantSetSVB = var.vcf_to_variantset(svFile)
    
        svPrefix = prefix + ".pbsv"
        svFile = tools.find_sv(subsetBam, param.REF_FA, svPrefix)
        variantSetSV = var.vcf_to_variantset(svFile)
    
    return (variantSetSVA, variantSetSVB, variantSetSV)

#    regionList = get_query_regions(rid, queryGraph, aligndf, seqData, param)
def get_query_regions(rid, queryGraph, aligndf, seqData, param):

    regionList = []
    
    alignments = aligndf[(aligndf["rid"] == rid)]
    alignments = alignments.sort_values(by=["rstart"]).reset_index()

    minSize=350

    for _,alnA in (alignments).iterrows():
        
        regionA = rgn.SimpleRegion(rid, alnA["rstart"], alnA["rend"])   
        hapId = supernova.find_haplotig(alnA["qid"], queryGraph)
    
        print(alnA["qid"], hapId)
        
        #todo: resolve non-hap overlaps
        #ex: tig00000446_pilon_pilon:6121909-6125469
        #supernova contigs 131183,126652

        #todo: multi-mapping supernova contigs, use divergance?????

        if hapId is not None:
            alnHap = alignments[alignments["qid"] == hapId]

            for _, alnB in (alnHap).iterrows():
                regionB = rgn.SimpleRegion(rid, alnB["rstart"], alnB["rend"])
                intersection = regionA & regionB
                if len(intersection) > 0 and len(intersection[0]) > minSize:
                    regions = rgn.region_difference(intersection[0], regionList, minSize=minSize)
                    
                    for region in regions:
                        #calls, candidates = resolve_variants_2haplotigs(alnA, alnB, seqData, param)
                        qregion = resolver.resolve_qvariants(alnA, alnB, region, seqData)
                        regionList.append(qregion)      
            
        singleHapRegions = rgn.region_difference(regionA, regionList, minSize=minSize)
        for region in singleHapRegions:
            qregion = resolver.resolve_qvariants(alnA, None, region, seqData)
            regionList.append(qregion)

    return regionList

'''
rid="tig00000446_pilon_pilon"
'''
def resolve_contig(rid, refGraph, queryGraph, aligndf, seqData, lengthData, param):
    dummyRun = True
  
    helper.contig_fasta(rid, seqData, param)
    rRegion = helper.contig_region(rid, lengthData, param)

    currentContigId = rid
    containedContigs, altPaths = refGraph.get_connected_contigs(currentContigId)
    
    svRegionList = []
    
    # 1) resolve SV contigs
    for svContig in containedContigs:              
        mainContig = currentContigId
        resolvedRegion = resolver.resolve_structural_haplotigs(mainContig, svContig, seqData, lengthData, param, dummy=dummyRun)
        svRegionList.append(resolvedRegion)


    # i=0 ; svRegion, qAlignment = svRegionList[i], containedContigs[i]
    for svRegion, svContig in zip(svRegionList, containedContigs):
            
        print(get_query_regions(svContig, queryGraph, aligndf, seqData, param))

        overlaps = rgn.region_overlap(svRegion, regionList)
        






    
    regionGaps = rgn.region_difference(rRegion, regionList)
    firstGap = regionGaps[0]
    lastGap = regionGaps[-1]
    
    for gap in regionGaps:
        #todo: resolve first and last "gap"
        fill_gap(gap, seqData, param, buffer1=30000, buffer2=1500)

    
    
    
    
    unitigId = refGraph.get_closest_unitig(rid, 0)
    ustart, uend = refGraph.unitig_pos(rid, unitigId)
    unitigRegion = SimpleRegion(rid, ustart, uend)
    
    #todo: utig coordinates are a bit off.....
    

    #       idx=3; alnA = alignments.iloc[idx]

    
        

    traversed = dict()
    





    while(True):
        print("------------------------------\n", "current unitig:", currentUnitigId)
        
        leftTigs, rightTigs = refGraph.get_unitig_neighbours(currentUnitigId)
        nextTigs = leftTigs if direction == LEFT else rightTigs
        traversed[currentUnitigId] = True
        #todo: update contigId

        traversed = backtrack(currentUnitigId, traversed, refGraph, direction)

        if len(nextTigs) == 0:
            results = zero_path(currentContigId, aligndf, traversed, refGraph, queryGraph, direction)
            currentContigId, currentUnitigId, direction = results

        elif len(nextTigs) == 1:
            currentUnitigId = one_path(currentUnitigId, rightTigs[0])
            
        elif len(nextTigs) == 2:
            results = two_path(currentUnitigId, nextTigs, aligndf, traversed, refGraph, direction)
            currentUnitigId, traversed = results
            
        else:
            three_path(currentContigId, currentUnitigId, nextTigs, aligndf, traversed, refGraph, queryGraph, direction)



        
'''   
rid="tig00000446_pilon_pilon"
alnA = alignments.iloc[0] ; alnB = alignments.iloc[1]

'''
def resolve_contigs(rid, refGraph, queryGraph, aligndf, seqData, param):
    
    unitigId = refGraph.get_closest_unitig(rid, 0)
    ustart, uend = refGraph.unitig_pos(rid, unitigId)
    unitigRegion = rgn.SimpleRegion(rid, ustart, uend)
    
    #todo: utig coordinates are a bit off.....
    
    alignments = aligndf[(aligndf["rid"] == rid) & 
                         (aligndf["rend"] >= ustart) & 
                         (aligndf["rstart"] <= uend)]
    alignments = alignments.sort_values(by=["rstart"]).reset_index()

    #       idx=3; alnA = alignments.iloc[idx]

    minSize=350
    regionList = []

    for _,alnA in (alignments).iterrows():
        
        regionA = rgn.SimpleRegion(rid, alnA["rstart"], alnA["rend"])   

        hapId = supernova.find_haplotig(alnA["qid"], queryGraph)
        
        print(alnA["qid"], hapId)
        
        #todo: resolve non-hap overlaps
        #ex: tig00000446_pilon_pilon:6121909-6125469
        #supernova contigs 131183,126652

        #todo: multi-mapping supernova contigs, use divergance?????


        if hapId is not None:
            alnHap = alignments[alignments["qid"] == hapId]

            for _, alnB in (alnHap).iterrows():
                regionB = rgn.SimpleRegion(rid, alnB["rstart"], alnB["rend"])
                intersection = regionA & regionB
                if len(intersection) > 0 and len(intersection[0]) > minSize:
                    regions = region_difference(intersection[0], regionList, minSize=minSize)
                    
                    for region in regions:
                        print("twohap", region)
                        regionList.append(region)
                        #calls, candidates = resolve_variants_2haplotigs(alnA, alnB, seqData, param)
            
        singleHapRegions = region_difference(regionA, regionList, minSize=minSize)
        for region in singleHapRegions:
            print("onehap", region)
            regionList.append(region)
            #calls, candidates = resolve_variants_2haplotigs(alnA, alnB, seqData, param)

    
    regionGaps = region_difference(unitigRegion, regionList)
    firstGap = regionGaps[0]
    lastGap = regionGaps[-1]
    
    for gap in regionGaps:
        #todo: resolve first and last "gap"
        fill_gap(gap, seqData, param, buffer1=30000, buffer2=1500)

    
def main_human_reference_polish():
    
    param = parameters.get_parameters_reference_polish("CF001")
    
    print("Parsing parameters...")
    seqData = dict()

    print("Reading reference fasta...")
    refData = io.read_fasta(param.REF_FA)
    seqData.update(refData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}

    backAlign = True
    plot = True

    #isolate region of interest
    region = rgn.region_from_string(param.REF_REGION, lengthData)

    for CFID in [ \
            #"CF001", "CF002", "CF003", "CF004", "CF006", "CF007", "CF010", "CF011", "CF013", "CF014", "CF016", 
            #"CF022", "CF024", 
            "CF045", "CF047", "CF049", "CF052", "CF060", "CF062", "CF063", "CF066", "CF067",
                 "CF071", "CF072", "CF073", "CF075", "CF076", "CF077", "CF078"]:

        param = parameters.get_parameters_reference_polish(CFID)    
        faDict, faDictOrder, bamDict = dict(), [], dict()
    
        prefix = helper.file_prefix(region, param)
        regionFa, alignedBam = helper.isolate_region(seqData, param.REF_ALIGNED_READS, region, prefix)
    
        faDict["HG38"] = regionFa
        faDictOrder.append("HG38")
    
        #polish with all reads
        polishedFa = polisher.iterative_polish(regionFa, alignedBam, prefix, niter=2, chunkSize=5000)
        polishedBam = tools.align_pacbio(polishedFa, alignedBam, prefix + "_polished")
        
        faDict["Polished"] = polishedFa
        faDictOrder.append("Polished")
    
        #phase reads
        vcf = tools.longshot_genotype(polishedBam, polishedFa, prefix, writeBams=True)
        bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
        
        bamDict["HapA"] = bamA
        bamDict["HapB"] = bamB
        bamDict["Unphased"] = bamUnphased
    
        fid, flen = helper.get_fasta_id(polishedFa), helper.get_fasta_len(polishedFa)
        tempRegion = rgn.SimpleRegion(fid, 0, flen-1)
        lengthData[fid] = flen
        
        phasedSegments, unphasedSegments = polisher.phaseblock_split(vcf, tempRegion, lengthData, param, collapseBuffer=10000)
    
        polishedFaA = polisher.iterative_polish(polishedFa, bamA, prefix + "hapA", niter=2, bamUnphased=bamUnphased, chunkSize=5000)
        polishedFaA2 = polisher.iterative_polish(polishedFaA, bamA, prefix + "hapA_final", niter=2, bamUnphased=bamUnphased)
    
        faDict["Haplotype A"] = polishedFaA2
        faDictOrder.append("Haplotype A")
    
        polishedFaB = polisher.iterative_polish(polishedFa, bamB, prefix + "hapB", niter=2, bamUnphased=bamUnphased, chunkSize=5000)
        polishedFaB2 = polisher.iterative_polish(polishedFaB, bamB, prefix + "hapB_final", niter=2, bamUnphased=bamUnphased)
    
        faDict["Haplotype B"] = polishedFaB2
        faDictOrder.append("Haplotype B")
    
        if backAlign:
            backFa = param.REF_FA
            #backFa = regionFa
            tools.align_pacbio(backFa, polishedFa, prefix + "_backaligned_polishAll")
            tools.align_pacbio(backFa, polishedFaA, prefix + "_backaligned_polishA1")
            tools.align_pacbio(backFa, polishedFaA2, prefix + "_backaligned_polishA2")
            tools.align_pacbio(backFa, polishedFaB, prefix + "_backaligned_polishB1")
            tools.align_pacbio(backFa, polishedFaB2, prefix + "_backaligned_polishB2")
    
        
        if plot:
            assesser.plot_alignments(faDict, bamDict, prefix, faDictOrder=faDictOrder, image="png")
            
        
        
    


    '''
    svBamA, svBamB, svBamUnphased = helper.phase_sv_reads(polishedFaA, bamA, polishedFaB, bamB, prefix, param, readsBamUnphased=bamUnphased)
    tools.align_pacbio(polishedFa, svBamA, prefix + "_A_svphased")
    tools.align_pacbio(polishedFa, svBamB, prefix + "_B_svphasedd")
    tools.align_pacbio(polishedFa, svBamUnphased, prefix + "_U_svphased")

    polishedBam = tools.align_pacbio(polishedFa, alignedBam, prefix + "_polished_haps")
    
    polishedSegments = []
    
    vcf = tools.longshot_genotype(param.REF_ALIGNED_READS, param.REF_FA, prefix, region=region, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    
    for phasedSegment in phasedSegments:
        polishedPhasedSegment = polisher.polish_region_haplotypes(phasedSegment, bamA, bamB, bamUnphased, seqData, lengthData, param, niter=2)
        polishedPhasedSegment = polisher.inject_pbsv_calls(polishedPhasedSegment, prefix, doA=True, doB=True)

        if backAlign:
            faA, faB = polishedPhasedSegment.write_seq(prefix)
            tools.align_pacbio(param.REF_FA, faA, prefix + "_A_polish_backaligned")
            tools.align_pacbio(param.REF_FA, faB, prefix + "_B_polish_backaligned")

        faA, faB = polishedPhasedSegment.write_seq(prefix)
        bamA, bamB, bamUnphased = helper.phase_sv_reads(faA, bamA, faB, bamB, prefix, param, readsBamUnphased=bamUnphased)
 

        if backAlign:
            faA, faB = polishedPhasedSegment.write_seq(prefix)
            tools.align_pacbio(param.REF_FA, faA, prefix + "_A_polish_backaligned")
            tools.align_pacbio(param.REF_FA, faB, prefix + "_B_polish_backaligned")

        polishedSegments.append(polishedPhasedSegment)

        
    for unphasedSegment in unphasedSegments:
        polishedUnphasedSegment = polisher.polish_region_unphased(unphasedSegment, bamA, bamB, bamUnphased, seqData, lengthData, param, dummy=dummyRun)
        polishedSegments.append(polishedUnphasedSegment)



    polishedSegments.sort(key=lambda x: x.start)
    i=0
    while True:
        if i < len(polishedSegments) -1:
            polish_junction(polishedSegments[i], polishedSegments[i+1], bamA, bamB, bamUnphased, seqData, lengthData, param)

    '''

    


def main():
    humanReferencePolish = True
    
    
    print("Parsing parameters...")
    param = parameters.get_parameters()
    seqData = dict()

    print("Reading reference fasta...")
    refData = io.read_fasta(param.REF_FA)
    seqData.update(refData)
    rids = list(refData.keys())

    if not humanReferencePolish:
        print("Reading Supernova fasta...")
        queryData, queryDescription = io.read_fasta(param.QUERY_FA, description=True)
        queryGraph = supernova.create_graph_df(queryDescription)
        refGraph = canu.CanuGraph(param.CANU_DIR)

        aligndf = paf.parse_paf(param.ALIGN)
        seqData.update(queryData)
        qids = list(refData.keys())

    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    rids.sort(key=lambda x: -lengthData[x])

    if not humanReferencePolish:
        qids.sort(key=lambda x: -lengthData[x])


    '''
    for rid in rids:

        contigger.resolve_path(rid, refGraph, queryGraph, aligndf, seqData, lengthData, param)






    #tig00983269_pilon_pilon   chr2:229,073,284-229,123,907
    #tig00000446_pilon_pilon   chr2:16,032,439-16,335,095
    
    for rid in rids:
        rid="tig00000446_pilon_pilon"
        logger.out("Doing contig: " + rid, 1, param)
                
        done = dict()
        regions = []
        
        aligndfSubset = aligndf[aligndf["rid"] == rid]
        for idx, row in (aligndfSubset).iterrows():
           
            qid, qstart, qend = row["qid"], row["qstart"], row["qend"]
            rid, rstart, rrFaend = row["rid"], row["rstart"], row["rend"]
            
            if qid in done: continue
            
            print(qid)
            if(qid == "290"):
                break

            hapId = supernova.find_haplotig(qid, queryGraph)
            
            #find easiest examples, 2 haplotigs + 1 to 1 alignments
            if hapId is None: continue
            alignments = aligndf[(aligndf["qid"] == qid) | (aligndf["qid"] == hapId)]
            
            rmatches = set(alignments["rid"])
            if len(rmatches) > 1: continue
            
            print(qid, " is good")
            done[qid] = True
            done[hapid] = True
            
            alnA = alignments[alignments["qid"] == qid].iloc[0]
            alnB = alignments[alignments["qid"] == hapid].iloc[0]

            region = find_variation(alnA, alnB, seqData, param)
            
            variantSet = region[0].get_supernova_variantset()
            variantSetSVA, variantSetSVB, variantSetSV = find_structural_variation(rid, rstart, rend, variantSet, param)
            #todo merge SV into small variant vcf
            
            
            if len(variantSetSVB) > 0 or len(variantSetSVA) > 0:
                print("!")
                input()
            
            
            print(variantSetSV)
            #todo: make sv for whole ref genome!
            regions.extend(region)
            
            
            
            
            
            
            
            
            
            
            
            if row["qid"] == "292":
                break
            
            

            
            
            
            
            
            
            #todo: logic to properly pair q and r contigs
            
            

            
            
            
            
            
            
            
            overlapping = aligndf[((aligndf["rstart"] <= rend) & (aligndf["rend"] >= rstart)) | \
                                  ((aligndf["rend"] <= rstart) & (aligndf["rstart"] >= rend))]
            
            haplotigID = supernova.find_haplotig(qid, queryGraph)
            if haplotigID is not None:
                haplotigAlignments = aligndfSubset[aligndfSubset["qid"] == haplotigID]
                haploRow = haplotigAlignments.iloc[0, :]
            else:
                haploRow = None
                
            #todo: order multiple alignments to haplotigs properly
                
            variantSetA = paf.parse_cs_string(row["cs"], rid, rstart, refData[rid])
            variantSetB = paf.parse_cs_string(haploRow["cs"], rid, haploRow["rstart"], refData[rid])            
            variantSet = var.combine_as_genotypes(variantSetA, variantSetB, refData[rid])
            
            whInput= param.OUTPUT_DIR + "/whatshap.in.vcf"
            
            variantSet.write_vcf(whInput, hetsOnly=True)

            
            subsetBam = tools.samtools_subset(param.REF_ALIGNED_READS, param.OUTPUT_DIR,
                            chrom=rid, startPos=rstart, endPos=rend)
            
            whVCF = param.OUTPUT_DIR + "/whatshap.out.vcf"
            inputs = [whInput, subsetBam]
            tools.whatshap_phase(whVCF, param.REF_FA, whInput, inputs, genotype=True, indels=True)

            variantSetWH = var.vcf_to_variantset(whVCF)
            
            
            whVCF = tools.bgzip(whVCF)

            whBAM = param.OUTPUT_DIR + "/whatshap.out.bam"            
            tools.whatshap_haplotag(whBAM, param.REF_FA, whVCF, subsetBam)                





            whOutgz = tools.bgzip(whOut)
            whBam = param.OUTPUT_DIR + "/whatshap.out.bam"
           
            leftBam = tools.samtools_subset(param.QUERY_ALIGNED_READS, param.OUTPUT_DIR,
                            chrom=leftTig, startPos=(lengthData[leftTig] - buffer), endPos=lengthData[leftTig])
        
            leftReads = tools.samtools_fetch(param.QUERY_ALIGNED_READS, leftTig, 
                           startPos=(lengthData[leftTig] - buffer), endPos=lengthData[leftTig])
        
            rightReads = tools.samtools_fetch(param.QUERY_ALIGNED_READS, rightTig, 
                           startPos=0, endPos=buffer)
            
            intersectingIds = set([a.qname for a in leftReads]).intersection([a.qname for a in rightReads])  
            '''
    

















if __name__== "__main__":
  main_human_reference_polish()
  print("done")
  #exit()