#various attemts to close this gap but it appears near impossible....
'''
# refTig = "tig00000871_pilon_pilon" ; scaffoldTig = "tig00000446_pilon_pilon" ; queryTig = "1175"
def scaffold_left(refTig, queryTig, refGraph, queryGraph, seqData, lengthData, param):
    
    #todo: resolve strand / positions
    
    unitigId = refGraph.get_closest_unitig(refTig, 0)
    leftUnitigs, rightUnitigs = refGraph.get_unitig_neighbours(unitigId)
    
    if len(leftUnitigs) > 0:
        print("different case, todo")
    
    #if no unitig to scaffold to... look to supernova
    aligndf[aligndf["qid"] == queryTig]
    
    #todo: logic to determine unambiguous connection between q and r
    #scaffoldTig = ....  
    
    buffer=5000
    
    #todo: resolve strand / positions using supernova alignment??? can we trust it??? todoodododod
    leftRegion = SimpleRegion(refTig, 0, buffer, lengthData)
    rightRegion = SimpleRegion(scaffoldTig, 0, buffer, lengthData)
    
        
    leftReads = tools.samtools_fetch(param.REF_ALIGNED_READS, leftRegion)
    rightReads = tools.samtools_fetch(param.REF_ALIGNED_READS, rightRegion)
    
    
    leftReadIds = set([read.qname for read in leftReads])
    rightReadIds = set([read.qname for read in leftReads])
    
    overlap = leftReadIds.intersection(rightReadIds)
    print("Left reads overlapping = %", 100*len(overlap)/len(leftReads))
    print("Right reads overlapping = %", 100*len(overlap)/len(rightReadIds))
    
    len(rightReads)/len(overlap)
    
    leftReadsSubset = [read for read in leftReads if read.qname in overlap]
    rightReadsSubset = [read for read in rightReads if read.qname in overlap]

    readDict = dict()
    for read in leftReads + rightReads:
        readDict[read.qname] = read
    

    gapFa = param.OUTPUT_DIR + "/gap_sequence.fasta"
    
    leftFlank = seqData[refTig][:buffer]
    leftFlank = reverse_complement(leftFlank)
    rightFlank = seqData[scaffoldTig][:buffer]
    
    writer = open(gapFa, 'w+')
    writer.write(">GAP\n" + leftFlank + rightFlank)
    writer.close()
    tools.samtools_faidx(gapFa)

    allReads = [readDict[qname] for qname in readDict]
    allReads = sorted(allReads, key=lambda x: (x.reference_name, x.pos))
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    alignedReads = tools.samtools_write(allReads, prefix, param.REF_ALIGNED_READS)
    realignedBam = tools.align_pacbio(gapFa, alignedReads, prefix)

    gapFa = param.OUTPUT_DIR + "/gap_sequence.fasta"

    tools.longshot_genotype(realignedBam, gapFa, prefix, writeBams=True)
    hpA, hpB, unaligned = tools.get_longshot_phased_reads(prefix)
    readsA = tools.samtools_fetch(hpA)
    readsB = tools.samtools_fetch(hpB)
    fastaFile = tools.reads2fasta(readsA, prefix)
    
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig
    
    #todo: estimate size better????
    genomeSize = buffer*2 + 5000
    tools.canu_correct(fastaFile, prefix, genomeSize, rawErrorRate=0.5, trim=False)
    #tools.canu_assemble(fastaFile, prefix, genomeSize)

    correctedReadsFasta = prefix +  "/" + refTig + "_" + scaffoldTig + ".correctedReads.fasta.gz"
    alignedCorrectedReadsBam  = tools.align_pacbio(gapFa, correctedReadsFasta, prefix + ".correctedB")
    alignedCorrectedReads = tools.samtools_fetch(alignedCorrectedReadsBam)

    #todo: pick best read.
    alignLen = {read.qname : [] for read in alignedCorrectedReads}
    readDict = {read.qname : read for read in alignedCorrectedReads}

    for read in alignedCorrectedReads:
        alignLen[read.qname].append(read.alen)
        
    alignLenList = [(readDict[qname], alignLen[qname]) for qname in alignLen]
    alignLenList = sorted(alignLenList, reverse=True, key=lambda x: (len(x[1]) == 2, len(x[1]), sum(x[1])))
    
    bestRead = alignLenList[0][0]
    
    gapFillFa = param.OUTPUT_DIR + "/gap_fill.fasta"

    writer = open(gapFillFa, 'w+')
    writer.write(">gap_fill" + "\n" + bestRead.seq + "\n")
    writer.close()

    prefix = param.OUTPUT_DIR + "/gap_fill.reads"
    readsBBam = tools.samtools_write(readsB,  prefix + ".gapfillreadsB", param.REF_ALIGNED_READS)

    numIter=10
    iteration = 0
    consensusFa = gapFillFa
    while iteration <= numIter:
        prefix = param.OUTPUT_DIR + "/gap_consensus" + str(iteration)
    
        realignedBamB = tools.align_pacbio(consensusFa, readsBBam, prefix)
        iteration += 1
        consensusFa = tools.pacbio_polish(realignedBamB, consensusFa, prefix, outputGff=True, outputFasta=True, useMapFilter=False)



    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    alignedReads = tools.samtools_write(readsB, prefix, param.REF_ALIGNED_READS)

    numIter=10
    iteration = 0
    consensusFa = gapFa
    while iteration <= numIter:
        prefix = param.OUTPUT_DIR + "/gap_fill_cons" + str(iteration)

        realignedBam = tools.align_pacbio(consensusFa, alignedReads, prefix)
        pbsvVCF = tools.find_sv(realignedBam, consensusFa, param.OUTPUT_DIR + "/pbsv" +str(iteration), threads=6)
        variantSetSV = var.vcf_to_variantset(pbsvVCF)
        
        consensusSeq = tools.fasta2dict(consensusFa, toUpper=True).popitem()[1]
        for variant in variantSetSV.get_variants():
            pos = variant.pos
            if "DUP" in variant.alt: continue
            consensusSeq[pos]
            consensusSeq = consensusSeq[:pos] + variant.alt + consensusSeq[pos+1:]
            
        consensusFa = prefix + ".pbsv.fasta"
        writer = open(consensusFa, 'w+')
        writer.write(">consensus" + str(iteration) + "\n" + consensusSeq + "\n")
        writer.close()
        realignedBam = tools.align_pacbio(consensusFa, alignedReads, prefix)
        iteration += 1

        consensusFa = tools.pacbio_polish(realignedBam, consensusFa, prefix, outputFasta=True, useMapFilter=False)




    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    alignedReads = tools.samtools_write(allReads, prefix, param.REF_ALIGNED_READS)
    realignedBam = tools.align_pacbio(gapFa, alignedReads, prefix)
    tools.find_sv(realignedBam, gapFa, param.OUTPUT_DIR + "/pbsv", threads=6)
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    fastaFile = tools.reads2fasta(allReads, prefix)
    
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig
    
    #todo: estimate size better????
    genomeSize = buffer*2 + 50000
    tools.canu_correct(fastaFile, prefix, genomeSize)
    #tools.canu_assemble(fastaFile, prefix, genomeSize)

    correctedReadsFa = prefix +  "/" + refTig + "_" + scaffoldTig + ".correctedReads.fasta.gz"
    realignedBam = tools.align_pacbio(consensusFa, alignedReads, prefix)

    flanksFa = param.OUTPUT_DIR + "/flanks.fasta"
    
    leftFlank = seqData[refTig][:buffer]
    rightFlank = seqData[scaffoldTig][:buffer]

    LEFT_FLANK = "leftflank"
    RIGHT_FLANK = "rightflank"

    writer = open(flanksFa, 'w+')
    writer.write(">" + LEFT_FLANK + "\n" + leftFlank + "\n")
    writer.write(">" + RIGHT_FLANK + "\n" + rightFlank)

    writer.close()

    alignedReads = tools.align_pacbio(flanksFa, correctedReadsFa, prefix)
    correctedReads = tools.samtools_fetch(alignedReads)

    #todo: pick best read.
    alignLen = {read.qname : [] for read in correctedReads}
    readDict = {read.qname : read for read in correctedReads}

    for read in correctedReads:
        alignLen[read.qname].append(read.alen)
        
    alignLenList = [(readDict[qname], alignLen[qname]) for qname in alignLen]
    alignLenList = sorted(alignLenList, reverse=True, key=lambda x: (len(x[1]) == 2, len(x[1]), sum(x[1])))
    
    bestRead = alignLenList[0][0]
    
    
    
    
    
    gapFa = param.OUTPUT_DIR + "/gap_fill.fasta"

    writer = open(gapFa, 'w+')
    writer.write(">gap_fill" + "\n" + bestRead.seq + "\n")
    writer.close()

    prefix = param.OUTPUT_DIR + "/gap_fill.reads"
    
    readDict = dict()
    for read in leftReads + rightReads:
        readDict[read.qname] = read
    
    readDict.pop(bestReadId)
        
    allReads = [readDict[qname] for qname in readDict]
    allReads = sorted(allReads, key=lambda x: (x.reference_name, x.pos))

    readsBam = tools.samtools_write(allReads, prefix, param.REF_ALIGNED_READS)
    
    prefix = param.OUTPUT_DIR + "/gap_fill.aligned"
    alignedReads = tools.align_pacbio(gapFa, readsBam, prefix, minConcordance=0)
    
    numIter=10
    iteration = 1

    while iteration <= numIter:
        prefix = param.OUTPUT_DIR + "/gap_consensus" + str(iteration)
    
        consensusFa = tools.pacbio_polish(alignedReads, gapFa, prefix, outputFasta=True, useMapFilter=False)
        realignedBam = tools.align_pacbio(consensusFa, alignedReads, prefix)
        iteration += 1
    
    
    
    










    prefix = param.OUTPUT_DIR + "/gap_fill.reads"
    

    
    
    
    
    
    #todo: pick best read.
    bestRead = sorted(leftReadsSubset, reverse=True, key=lambda x: x.query_length)[0]
    bestReadId = bestRead.qname

    bestLeftAlignment = [read for read in leftReads if read.qname == bestReadId][0]
    bestLeftPos = (bestLeftAlignment.qstart, bestLeftAlignment.qend)
    if bestLeftAlignment.is_reverse:
        length = bestLeftAlignment.query_length
        bestLeftPos = (length - bestLeftAlignment.qend, length - bestLeftAlignment.qstart)
    
    bestRightAlignment = [read for read in rightReads if read.qname == bestReadId][0]
    bestRightPos = (bestRightAlignment.qstart, bestRightAlignment.qend)
    if bestRightAlignment.is_reverse:
        length = bestLeftAlignment.query_length
        bestRightPos = (length - bestRightAlignment.qend, length - bestRightAlignment.qstart)

    subreadPos = (min(max(bestRightPos), max(bestLeftPos)) - 100, \
                  max(min(bestRightPos), min(bestLeftPos)) + 100)
    subread = bestLeftAlignment.get_forward_sequence()[subreadPos[0]:subreadPos[1]]
    
    
    #todo: get right strand
    #leftFlank = seqData[refTig][:buffer]
    #leftFlank = reverse_complement(q1)
    #rightFlank = seqData[scaffoldTig][:buffer]    
    
    
    gapFa = tools.reads2fasta([bestRead], prefix)
    
    prefix = param.OUTPUT_DIR + "/gap_fill.reads"
    readsBam = tools.samtools_write(rightReads, prefix, param.REF_ALIGNED_READS)
    
    prefix = param.OUTPUT_DIR + "/gap_fill.aligned"
    alignedReads = tools.align_pacbio(gapFa, readsBam, prefix, minConcordance=0)
    
    
    
    
    
    
    
    #todo: is this enough to do assembly??
   

    
    
    
    
    
    
    newFa  = param.OUTPUT_DIR + "/trying_this.fasta"
    writer = open(newFa, 'w+')
    q1 = seqData[refTig][:buffer]
    q1rc = reverse_complement(q1)
    q2 = seqData[scaffoldTig][:buffer]


    writer.write(">testing" + "\n" + q1+q2)
    writer.close()
    
    leftReadsSubset = sorted(leftReadsSubset, key=lambda x: x.pos)
    rightReadsSubset = sorted(rightReadsSubset, key=lambda x: x.pos)
   
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    bam = tools.samtools_write(rightReads, file_prefix(rightRegion, param), param.REF_ALIGNED_READS)
    
    tools.align_pacbio(newFa, bam, param.OUTPUT_DIR + "/trying_this.align" )
    


    for alignment in alignments:
        fastaLines[alignment.qname] = alignment.seq
    for readId in fastaLines:
        writer.write(">" + readId + "\n" + fastaLines[readId] + "\n")
    writer.close()

    
    
    leftReadsSubset = [read for read in leftReads if read.qname in overlap]
    rightReadsSubset = [read for read in rightReads if read.qname in overlap]
    leftReadsSubset = sorted(leftReadsSubset, key=lambda x: x.pos)
    rightReadsSubset = sorted(rightReadsSubset, key=lambda x: x.pos)
   
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig + ".reads"
    fastaFile = tools.reads2fasta(leftReads + rightReads, prefix)
    
    prefix = param.OUTPUT_DIR + "/" + refTig + "_" + scaffoldTig

    #todo: estimate size better????
    genomeSize = buffer*2 + 50000
    tools.canu_correct(fastaFile, prefix, genomeSize)
    tools.canu_assemble(fastaFile, prefix, genomeSize)






    leftBam = tools.samtools_write(leftReadsSubset, file_prefix(leftRegion, param), param.REF_ALIGNED_READS)
    rightBam = tools.samtools_write(rightReadsSubset, file_prefix(rightRegion, param), param.REF_ALIGNED_READS)




    
   startPos=(lengthData[leftTig] - buffer), endPos=lengthData[leftTig])

    
    leftBam = tools.samtools_subset(param.REF_ALIGNED_READS, leftRegion, leftPrefix)
    
                    chrom=leftTig, startPos=(lengthData[leftTig] - buffer), endPos=lengthData[leftTig])

    leftReads = tools.samtools_fetch(param.REF_ALIGNED_READS, leftTig, 
                   startPos=(lengthData[leftTig] - buffer), endPos=lengthData[leftTig])

    rightReads = tools.samtools_fetch(param.QUERY_ALIGNED_READS, rightTig, 
                   startPos=0, endPos=buffer)
    
    intersectingIds = set([a.qname for a in leftReads]).intersection([a.qname for a in rightReads])  




'''