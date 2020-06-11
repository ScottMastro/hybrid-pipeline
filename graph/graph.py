import os, sys, re
import shutil, glob
import gzip
import argparse
import pyfaidx
from Bio import SeqIO
import subprocess 
import random, copy

def reset_directory(directory):
    try:
        shutil.rmtree(directory)
    except:
        pass
    
    try:
        os.mkdir(directory)
    except:
        pass

def cat(directory, extension, outfile, others=[]):
    with open(outfile, 'wb') as out:
        for other in others:
            with open(other, 'rb') as readfile:
                shutil.copyfileobj(readfile, out)
                
        for filename in glob.glob(directory + "*" + extension):
            if filename == outfile:
                # don't want to copy the output into the output
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, out)
                
    return outfile
    
def read_fasta(fasta, region=None):
    
    if region is not None:
        try:
            fa = pyfaidx.Faidx(fasta)
            chrom = region.split(":")[0]
            start = int(region.split(":")[1].split("-")[0])
            end = int(region.split(":")[1].split("-")[1])
        except:
            print("Could not parse region subset (\"" + region + "\"). Format as chr:start-end")
            sys.exit()
        
        faDict = {region : str(fa.fetch(chrom, start, end))}
    
    else:
        if(fasta[-2:] == "gz"):
            with gzip.open(fasta, "rt") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            records = list(SeqIO.parse(fasta, "fasta"))
        
        faDict = dict(zip([str(r.id) for r in records], [str(r.seq).upper() for r in records]))

    return faDict


def write_fasta(faPath, fastaDict):

    if not faPath.endswith(".fa") and not faPath.endswith(".fasta"):
        faPath = faPath + ".fasta"
        
    writer = open(faPath, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")
            
    writer.close()
    return faPath


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
        return write_fasta(faPath, fastaDict)
                            
def main():

    # default parameters
    INPUT = "/home/scott/Dropbox/hybrid-pipeline/graph/input.txt"
    OUTDIR = "/media/scott/Zapdos/out/all/out"
    REF = "/media/scott/Zapdos/out/all/chr5_polish/chr5_393462_677667.fasta"
    REGION=None
    #REGION="chr5_393462_677667:1-100"
    X = 200
    Y = 20000
    SKIP = False

    parser = argparse.ArgumentParser(description="Local Graph Constructor")

    parser.add_argument("-i", "--input", type=str, default=INPUT,
            help="Tab-separated list of input FASTA. First column is sample ID, second column is FASTA directory")
    parser.add_argument("-o", "--outdir", type=str, default=OUTDIR,
                help="Directory where output will be written." )
    parser.add_argument("-r", "--ref", type=str, default=REF,
                help="FASTA file with sequence to target." )
    parser.add_argument("-s", "--subset", type=str, default=REGION,
            help="Extract a subsequence from ref FASTA as target (optional). Format as chr:start-end." )

    parser.add_argument("-x", "--filter", type=str, default=X,
            help="Filter alignments smaller than this distance in bp." )
    parser.add_argument("-y", "--collapse", type=str, default=Y,
            help="Combine alignments that are separated by a distance smaller than this threshold in bp." )

    parser.add_argument("-k", "--skip", type=bool, default=SKIP,
            help="Skip over FASTA extraction step and build graph from existing PAF files." )

    args = parser.parse_args()
    #set parameters from user input
    INPUT = args.input
    OUTDIR = args.outdir
    OUTDIR = OUTDIR + ("" if OUTDIR[-1] == "/" else "/")
    REF = args.ref
    SKIP = args.skip

    X = args.filter
    Y = args.collapse

    try:
        REGION = args.subset
    except:
        REGION = None

    if not SKIP:

        if not os.path.exists(OUTDIR): os.mkdir(OUTDIR)
        
        targetFa = OUTDIR + "target.fasta"
        sequences = read_fasta(REF, REGION)
        write_fasta(targetFa, sequences)
        
        E_PAFDIR=OUTDIR + "extract_paf/"    
        reset_directory(E_PAFDIR)
        FASTADIR=OUTDIR + "fasta/"    
        reset_directory(FASTADIR)
    
        extractedFasta = []
        for line in open(INPUT, "r"):
            fid, fdir = line.rstrip().split("\t")
            
            pafOut = E_PAFDIR + fid + ".paf"
            subprocess.run(["bash", "extract_minimap.sh", targetFa, fdir, pafOut])
    
            extract = extract_from_paf(pafOut, X, Y, fdir, FASTADIR, fid)
            if extract is not None:
                extractedFasta.append(extract)
    
        C_PAFDIR=OUTDIR + "cross_paf/"    
        reset_directory(C_PAFDIR)
    
        shuffledFasta = copy.deepcopy(extractedFasta)
        random.seed(776)
        random.shuffle(shuffledFasta)
    
        while len(shuffledFasta) > 1:
            extractedFa = shuffledFasta.pop()
            fid = os.path.splitext(os.path.basename(extractedFa))[0]
            print(fid, "v", "ref")
    
            pafOut = C_PAFDIR + "ref_" + fid + ".paf"
            subprocess.run(["bash", "cross_minimap.sh", targetFa, extractedFa, pafOut])
    
            for otherFasta in shuffledFasta:
                otherFid = os.path.splitext(os.path.basename(otherFasta))[0]
                print(fid, "v", otherFid)
    
                pafOut = C_PAFDIR + otherFid + "_" + fid + ".paf"
                subprocess.run(["bash", "cross_minimap.sh", otherFasta, extractedFa, pafOut])
    
    allFasta = OUTDIR + "all.fasta"
    allPaf = OUTDIR + "all.paf"
    
    cat(FASTADIR, ".fa*", allFasta, others=[targetFa])
    cat(C_PAFDIR, ".paf", allPaf)

    graphPrefix = OUTDIR + "graph"
    subprocess.run(["bash", "build_graph.sh", allFasta, allPaf, graphPrefix])


if __name__== "__main__":
  main()