import os, sys, re
import shutil, glob
import gzip
import pyfaidx
from Bio import SeqIO

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
                            
