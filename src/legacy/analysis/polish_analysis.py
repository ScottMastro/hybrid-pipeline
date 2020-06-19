def reference_analysis(faFile, refReads, outdir, param):
    
    outName="analysis"
    refBam = tools.align_pacbio(faFile, refReads, outdir + outName)
    alignments = tools.samtools_fetch(refBam)


    alignments = tools.samtools_fetch(param.REF_ALIGNED_READS, rRegion)
    analyzer

    return refBam, queryBam

def analyze_contig(refReads, queryReads, unpolishedFa, consensusFa, hap1Fa, hap2Fa, outDir, param):
    
    d="/media/scott/Zapdos/out/scaff139_/"
    outDir=d
    consensusFa=outDir+"consensus.fasta"     
    unpolishedFa=outDir+"unpolished.fasta" 
    consensusFa=outDir+"consensus.fasta" 
    hap1Fa=outDir+"hap1.fasta" 
    hap2Fa=outDir+"hap2.fasta" 

    refReads=outDir+"pacbio_reads.bam" 
    queryReads=(outDir+"10x_reads/bamtofastq_S1_L000_R1_001.fastq.gz", outDir+"10x_reads/bamtofastq_S1_L000_R2_001.fastq.gz")
    faFile=unpolishedFa

    phasedBam=outDir+"ref_to_consensus.pbmm2.haplotag.bam"
    phasedAlns = tools.samtools_fetch(phasedBam)
    h1,h2,unknown = [],[],[]
    for aln in phasedAlns:
        try:
            hp=str(aln.get_tag("HP"))
            if hp == "1":
                h1.append(aln)
            elif hp == "2":
                h2.append(aln)
            else:
                unknown.append(aln)
        except:
            unknown.append(aln)

    tempDir = outDir + ("/" if outDir[-1] != "/" else "") + "temp/"
    io.make_dir(tempDir)

    h1RefBam = tools.samtools_write(h1, tempDir + "h1_ref", phasedBam, index=True)
    h2RefBam = tools.samtools_write(h2, tempDir + "h2_ref", phasedBam, index=True)
    h0RefBam = tools.samtools_write(unknown, tempDir + "h0_ref", phasedBam, index=True)

    files = [unpolishedFa, consensusFa, hap1Fa, hap2Fa]
    names = ["Unpolished", "Consensus", "Hap1", "Hap2"]
    
    refBamFiles, queryBamFiles = [], []
    
    concordanceList, concordanceListName = [], []
    for faFile, name in zip(files, names):
        haplotypes = ["1","2","0"]
        bamFiles = []
        for haplotypeBam, hp in zip([h1RefBam, h2RefBam, h0RefBam], haplotypes):
    
            if not os.path.isfile(tempDir + name + "." + hp + ".pbmm2.bam"):
                rBam = tools.align_pacbio(faFile, haplotypeBam, tempDir + name + "." +hp)
                #analyzer.plot_coverage(rBam, title=None)
            else:
                rBam = tempDir + name + "." + hp + ".pbmm2.bam"
            refBamFiles.append(rBam)
            concordanceList.append(analyzer.get_concordance(rBam, hap=hp))
            concordanceListName.append(name)
            bamFiles.append(rBam)
            
        #analyzer.get_coverage(bamFiles, haplotypes, outDir)

            
            
    analyzer.plot_concordance(concordanceList, concordanceListName, outDir)

