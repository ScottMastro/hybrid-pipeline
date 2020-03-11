import parameters
import file_handler as io
import canu_helper as canu
import regions as rgn
import helper
import polish_helper as polisher

'''
rid="tig00000446_pilon_pilon"
rid="tig00007781_pilon_pilon"

'''
def polish_contig(rid, refGraph, seqData, lengthData, param):
    dummyRun = True
  
    helper.contig_fasta(rid, seqData, param)
    contigRegion = helper.contig_region(rid, lengthData, param)

    currentContigId = rid
    containedContigs, altPaths = refGraph.get_connected_contigs(currentContigId)
    
    svSegments = []
    
    # 1) resolve SV contigs
    for svContig in containedContigs:
        mainContig = currentContigId
        svSegment = polisher.polish_structural_haplotigs(mainContig, svContig, seqData, lengthData, param, dummy=dummyRun)
        svSegments.append(svSegment)
        
    regions = rgn.region_difference(contigRegion, svSegments)

    # 2) polish rest of contig
    polishedSegments = []
    polishedSegments.extend(svSegments)
    
    for region in regions:
        
        vcf,bamA,bamB,bamUnphased = polisher.phase_region(param.REF_FA, param.REF_ALIGNED_READS, region, param)
        phasedSegments, unphasedSegments = polisher.phaseblock_split(vcf, region, lengthData, param)
        
        for phasedSegment in phasedSegments:
            polishedPhasedSegment = polisher.polish_region_haplotypes(phasedSegment, bamA, bamB, bamUnphased, seqData, lengthData, param, dummy=dummyRun)
            polishedSegments.append(polishedPhasedSegment)
            
        for unphasedSegment in unphasedSegments:
            polishedUnphasedSegment = polisher.polish_region_unphased(unphasedSegment, bamA, bamB, bamUnphased, seqData, lengthData, param, dummy=dummyRun)
            polishedSegments.append(polishedUnphasedSegment)

    polishedSegments.sort(key=lambda x: x.start)
    i=0
    while True:
        if i < len(polishedSegments) -1:
            polish_junction(polishedSegments[i], polishedSegments[i+1], bamA, bamB, bamUnphased, seqData, lengthData, param)

        
        
    
        
        
        
def polish_main():
    print("Parsing parameters...")
    param = parameters.get_parameters()

    print("Reading Canu data...")
    refData = io.read_fasta(param.REF_FA)
    refGraph = canu.CanuGraph(param.CANU_DIR)
   
    seqData = dict()
    seqData.update(refData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    rids = list(refData.keys())
    rids.sort(key=lambda x: -lengthData[x])

    for rid in rids:
        polish_contig(rid, refGraph, seqData, lengthData, param)

