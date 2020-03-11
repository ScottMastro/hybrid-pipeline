import parameters
import file_handler as io
import external_tools as tools
import regions as rgn
import helper
import paf_helper
import polish_helper as polisher
import assessment as assesser
import variants as var

def make_haplo_vcfs():
    
    param = parameters.get_parameters_reference_polish()

    print("Parsing parameters...")
    seqData = dict()

    print("Reading reference fasta...")
    refData = io.read_fasta(param.REF_FA)
    seqData.update(refData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}

    #isolate region of interest
    region = rgn.region_from_string(param.REF_REGION, lengthData)


    for CFID in [ \
            "CF001", "CF002", "CF003", "CF004", "CF006", "CF007", "CF010", "CF011", "CF013", "CF014", "CF016", 
            "CF022", "CF024", "CF045", "CF047", "CF049", "CF052", "CF060", "CF062", "CF063", "CF066", "CF067",
            "CF071", "CF072", "CF073", "CF075", "CF076", "CF077", "CF078"]:

        param = parameters.get_parameters_reference_polish(CFID)        
        prefix = helper.file_prefix(region, param)
        regionFa, alignedBam = helper.isolate_region(seqData, param.REF_ALIGNED_READS, region, prefix)
        hapAFa = prefix + "hapA_final.polished.fasta"
        hapBFa = prefix + "hapB_final.polished.fasta"
            
        pafFileA = tools.align_paf_very_lenient(regionFa, hapAFa, prefix + "_hapA")
        pafFileB = tools.align_paf_very_lenient(regionFa, hapBFa, prefix + "_hapB")

        variantSetA = paf_helper.get_variantset(pafFileA, regionFa, hapAFa)
        variantSetB = paf_helper.get_variantset(pafFileB, regionFa, hapBFa)
        
        variantSet = var.combine_as_genotypes(variantSetA, variantSetB, phased=True)
        vcf = variantSet.write_vcf(prefix + "_haps.vcf")
        variantSet.shift_all(region.start, region.chrom)
        vcfRef = variantSet.write_vcf(prefix + "_hg38_haps.vcf")


        phasedVcf = tools.whatshap_phase([alignedBam], vcf, regionFa, prefix, indels=False)
        variantSet = var.vcf_to_variantcallset(phasedVcf)

        variantSet.shift_all(region.start, region.chrom)
        variantSet.write_vcf(phasedVcf)

        '''
        vcfA = variantSetA.write_vcf(prefix + "_hapA.vcf")
        vcfB = variantSetB.write_vcf(prefix + "_hapB.vcf")

        vcfAPhased = tools.whatshap_phase([alignedBam], vcfA, regionFa, prefix + "_A", indels=False)
        vcfBPhased = tools.whatshap_phase([alignedBam], vcfB, regionFa, prefix + "_B", indels=False)
        
        whVariantSetA = var.vcf_to_variantcallset(vcfAPhased)
        whVariantSetB = var.vcf_to_variantcallset(vcfBPhased)

        whVariantSetA.shift_all(region.start, region.chrom)
        whVariantSetB.shift_all(region.start, region.chrom)

        vcfA = whVariantSetA.write_vcf(prefix + "_A.whatshap.phase.vcf")
        vcfB = whVariantSetB.write_vcf(prefix + "_B.whatshap.phase.vcf")
        '''

def main_human_reference_polish():
    
    param = parameters.get_parameters_reference_polish()
    
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
            "CF001", "CF002", "CF003", "CF004", "CF006", "CF007", "CF010", "CF011", "CF013", "CF014", "CF016", 
            "CF022", "CF024", "CF045", "CF047", "CF049", "CF052", "CF060", "CF062", "CF063", "CF066", "CF067",
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
            
        
        
    

if __name__== "__main__":
  make_haplo_vcfs() 
  #main_human_reference_polish()
  print("done")
  #exit()