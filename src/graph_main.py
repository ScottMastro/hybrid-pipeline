import sys
import os
import pyfaidx
import random, copy

import utils.log as logger
import utils.file_handler as io
import utils.fasta_handler as fasta
import utils.external_tools as tools
  
def extract_from_fasta(fasta, region=None):
    
    if region is not None:
        try:
            fa = fasta.faidx_fetch(fasta, region, toUpper=True)
            return {str(region) : fa}
        except:
            print("Could not parse region subset (\"" + region + "\"). Format as chr:start-end")
            sys.exit()    
    else:
        return io.read_fasta(fasta, toUpper=True)

def extract_from_paf(paf, x, y, fasta, outdir, sampleID):
    regions = []    
    for line in open(paf, "r"):
        split = line.rstrip().split("\t")
        chrom, start, end = split[0], int(split[2]), int(split[3])
        
        #filter small alignments
        if abs(end - start) < x: continue
        regions.append((chrom, start, end))
        
    if len(regions) > 1:

        regions = sorted(regions, key=lambda tup: tup[2])
        rdict = {}
        for rgn in regions:
            if rgn[0] not in rdict: rdict[rgn[0]] = []
            rdict[rgn[0]].append(rgn)
            
        regions = []
        for chrom in rdict:
            if len(rdict[chrom]) == 1:
                regions.append(rdict[chrom][0])
            else:
                prevRgn = rdict[chrom][0]
                for rgn in rdict[chrom][1:]:
                    
                    #collapse close alignments
                    if rgn[1] - prevRgn[2] < y:
                        prevRgn = (chrom, min(rgn[1], prevRgn[1]), max(rgn[2], prevRgn[2]))
                    else:
                        regions.append(prevRgn)
                        prevRgn = rgn
                regions.append(prevRgn)
                    
        #==============================
        
        if len(regions) < 1: return None
        
        fa = pyfaidx.Faidx(fasta)
        fastaDict = {}

        for rgn in regions:
            name = "_".join([sampleID]+[str(r) for r in rgn])
            fastaDict[name] = str(fa.fetch(rgn[0], rgn[1]+1, rgn[2]+1))
                        
        faPath = outdir + sampleID + ".fasta"
        return fasta.write_fasta(faPath, fastaDict)
              
def build_graph(fasta, paf, prefix, param):
       
    gfa = tools.seqwish_graph(fasta, paf, prefix, minMatchLen=16)
    
    #for visualization tools:
    vg = tools.vg_normalize_gfa(gfa, prefix, generateGfa=True) 
    xg = tools.vg_index(vg, xg=True) 
    tools.vg_path_to_gam(xg, prefix)     
    
def main(param):
    
    #==================================================
    # Extract target region 
    #==================================================

    outdir = param.OUTPUT_DIR + ("" if param.OUTPUT_DIR[-1] == "/" else "/")
    logger.Logger(clean=True, level=param.VERBOSE)
    logger.FileLogger(clean=True, outdir=outdir)

    try:
        targetRegion = param.TARGET_REGION
    except:
        targetRegion = None

    if not param.SKIP_EXTRACTION:
        
        logger.log("Extracting target sequence(s)...")
        io.make_dir(outdir)
        
        targetFa = outdir + "target.fasta"
        sequences = extract_from_fasta(param.REF_GENOME, targetRegion)
        fasta.write_fasta(targetFa, sequences)
        
        ePafDir = outdir + "extract_paf/"    
        io.reset_directory(ePafDir)
        fastaDir = outdir + "fasta/"    
        io.reset_directory(fastaDir)
        
        logger.log("Extracting target sequence(s)...")

        extractedFasta = []
        for line in open(param.TSV_INPUT, "r"):
            
            fid, fdir = line.rstrip().split("\t")
            logger.log("Aligning target to " + fid + ".", logger.LOG_DETAILS)

            prefix = ePafDir + fid
            #minimap2 -cx asm5 --secondary=no $TARGET_FA $QUERY_FA > $OUT_PAF          
            pafOut = tools.align_paf(targetFa, fdir, prefix, asm="asm5", secondary=False)
                
            extract = extract_from_paf(pafOut, param.FILTER_X, param.COMBINE_Y, fdir, fastaDir, fid)
            if extract is not None:
                extractedFasta.append(extract)
    
    #==================================================
    # Align input sequences
    #==================================================

        cPafDir = outdir + "cross_paf/"    
        io.reset_directory(cPafDir)
    
        shuffledFasta = copy.deepcopy(extractedFasta)
        random.seed(428)
        random.shuffle(shuffledFasta)
    
        while len(shuffledFasta) > 1:
            extractedFa = shuffledFasta.pop()
            fid = os.path.splitext(os.path.basename(extractedFa))[0]
            print(fid, "v", "ref")
              
            prefix = cPafDir + "ref_" + fid
            #minimap2 -cx asm5 --cs --secondary=no $TARGET_FA $QUERY_FA > $OUT_PAF
            pafOut = tools.align_paf(targetFa, extractedFa, prefix, asm="asm5", secondary=False)

    
            for otherFasta in shuffledFasta:
                otherFid = os.path.splitext(os.path.basename(otherFasta))[0]
                print(fid, "v", otherFid)
    
                prefix = cPafDir + otherFid + "_" + fid
                #minimap2 -cx asm5 --cs --secondary=no $TARGET_FA $QUERY_FA > $OUT_PAF
                pafOut = tools.align_paf(otherFasta, extractedFa, prefix, asm="asm5", secondary=False)
    
    #==================================================
    # Generate graph
    #==================================================

    allFasta = outdir + "all.fasta"
    allPaf = outdir + "all.paf"
    
    io.cat(fastaDir, ".fa*", allFasta, others=[targetFa])
    io.cat(cPafDir, ".paf", allPaf)

    graphPrefix = outdir + "graph"
    build_graph(allFasta, allPaf, graphPrefix, param)