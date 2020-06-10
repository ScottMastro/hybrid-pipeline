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


'''
import argparse





TARGET_FA=$1
OUTDIR="${2:-"./extract"}"


clean () {
# 1=in, 2=out, 3=new fasta id

  if [ -z "$var" ]
    then
      awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $1 > $2
    else
      awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $1 | \
           sed 's/>.*/${3}/' > $2
  fi
}


print_something Mars
print_something Jupiter








echo $OUTDIR

exit 1

cfid=CF011
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${cfid}/hybrid_assembly.fasta | sed 's/>.*/&_CF011/' > all/${cfid}.fasta

#normalize lines
gawk 'BEGIN { RS=">";FS="\n" } NR > 1 { print ">"$1;print gensub(/.{65}/,"&\n","g",gensub(/\n/,"","g",substr($0,length($1)+1))) }' fasta 

minimap2 -cx asm20 --secondary=no CF011.fasta CF019.fasta > 11_19.paf
cat CF011.fasta CF019.fasta > 11_19.fasta
seqwish --repeat-max=1 -s 11_19.fasta -p 11_19.paf -g 11_19.gfa
vg view -Fv all.gfa > all_unmod.vg

#-------------------
cfid=CF052
minimap2 -cx asm20 --secondary=no ../${cfid}.fasta chr5_393462_677667.fasta > chr5_${cfid}.paf
#-------------------
minimap2 -cx asm5 --secondary=no chr5_393462_677667.fasta cf011.fasta > cf011.paf
minimap2 -cx asm5 --secondary=no chr5_393462_677667.fasta cf019.fasta > cf019.paf
minimap2 -cx asm5 --secondary=no chr5_393462_677667.fasta cf045.fasta > cf045.paf
minimap2 -cx asm5 --secondary=no chr5_393462_677667.fasta cf062.fasta > cf062.paf

#19
# -

#11
minimap2 -cx asm5 --secondary=no cf011.fasta cf019.fasta > cf011_cf019.paf

#45
minimap2 -cx asm5 --secondary=no cf045.fasta cf011.fasta > cf045_cf011.paf
minimap2 -cx asm5 --secondary=no cf045.fasta cf019.fasta > cf045_cf019.paf

#62
minimap2 -cx asm5 --secondary=no cf062.fasta cf011.fasta > cf062_cf011.paf
minimap2 -cx asm5 --secondary=no cf062.fasta cf019.fasta > cf062_cf019.paf
minimap2 -cx asm5 --secondary=no cf062.fasta cf045.fasta > cf062_cf045.paf

#52
minimap2 -cx asm5 --secondary=no cf052.fasta cf011.fasta > cf052_cf011.paf
minimap2 -cx asm5 --secondary=no cf052.fasta cf019.fasta > cf052_cf019.paf
minimap2 -cx asm5 --secondary=no cf052.fasta cf045.fasta > cf052_cf045.paf
minimap2 -cx asm5 --secondary=no cf052.fasta cf062.fasta > cf052_cf062.paf

#72
minimap2 -cx asm5 --secondary=no cf072.fasta cf011.fasta > cf072_cf011.paf
minimap2 -cx asm5 --secondary=no cf072.fasta cf019.fasta > cf072_cf019.paf
minimap2 -cx asm5 --secondary=no cf072.fasta cf045.fasta > cf072_cf045.paf
minimap2 -cx asm5 --secondary=no cf072.fasta cf062.fasta > cf072_cf062.paf
minimap2 -cx asm5 --secondary=no cf072.fasta cf052.fasta > cf072_cf052.paf

mkdir -p pafs
list=()
for d in ./sequence/* ; do  
    cfid=`basename $d`
    filename="${cfid%.*}"
    echo $filename
    #echo "${list[@]}"
    
    for e in "${list[@]}" ; do
        eid=`basename $e`
        ename="${eid%.*}"
	minimap2 -cx asm5 --secondary=no --cs $e $d > pafs/${ename}_${filename}.paf
	#echo "minimap2 -cx asm5 --secondary=no $e $d > pafs/${ename}_${filename}.paf"
    done
    
    list+=($d)
    minimap2 -cx asm5 --secondary=no --cs chr5_393462_677667.fasta $d > pafs/ref_${filename}.paf
    
done

cat pafs/*.paf > all.paf
rm all.fasta
cat *.fasta > all.fasta

seqwish -k 16 -s all.fasta -p all.paf -g all.gfa
vg view -Fv all.gfa | vg mod -X 256 - | vg mod -n - | vg ids -c -| vg mod -X 256 - > all.vg
vg stats -z all.vg  ; vg paths -L -v all.vg
vg view -Vg all.vg > normalized.gfa

vg index -x all.xg all.vg
vg paths -X -x all.xg > all.gam

#vg find -x all.xg -n 9105 -c 3 | vg view -dp - | dot -Tsvg -o subgraph.svg

vg find -N nodes.txt -x all.xg | vg paths -F -v -



#-------------------

mkdir -p vcfs
list=()
for d in ./pafs/ref* ; do  
    cfid=`basename $d`
    filename="${cfid%.*}"
    echo $filename

    sort -k6,6 -k8,8n $d | \
       paftools call -s $filename -f chr5_393462_677667.fasta - > ./vcfs/${filename}.vcf
    bgzip -f ./vcfs/${filename}.vcf
    tabix ./vcfs/${filename}.vcf.gz
done

bcftools merge -0 ./vcfs/*.gz > all.ref.vcf
vg construct -r chr5_393462_677667.fasta -v all.ref.vcf > all.ref.vg
vg mod -X 256 all.ref.vg | vg mod -n - | vg ids -c -| vg mod -X 256 - > normalized.ref.vg
vg stats -z normalized.ref.vg  ; vg paths -L -v normalized.ref.vg
vg view -Vg normalized.ref.vg > normalized.ref.gfa

#-------------------

#who aligns to who???
minimap2 -cx asm20 --secondary=no chr3_195721919_195989697.fasta scaff286_CF011.fasta > aln.paf
minimap2 -cx asm20 --secondary=no chr3_195721919_195989697.fasta scaff208_CF019.fasta > aln2.paf

minimap2 -cx asm20 --secondary=no scaff286_CF011.fasta chr3_195721919_195989697.fasta > aln.paf
minimap2 -cx asm20 --secondary=no scaff208_CF019.fasta chr3_195721919_195989697.fasta  > aln2.paf


cat chr3_195721919_195989697.fasta scaff286_CF011.fasta > aln.fasta
seqwish --repeat-max=1 -s aln.fasta -p aln.paf -g aln.gfa
vg view -Fv aln.gfa > aln_unmod.vg
vg mod -X 256 aln_unmod.vg > aln.vg

cat chr3_195721919_195989697.fasta scaff208_CF019.fasta > aln2.fasta
seqwish --repeat-max=1 -s aln2.fasta -p aln2.paf -g aln2.gfa
vg view -Fv aln2.gfa > aln2_unmod.vg
vg mod -X 256 aln2_unmod.vg > aln2.vg

vg stats -z aln.vg  ; vg paths -L -v aln.vg
vg stats -z aln2.vg ; vg paths -L -v aln2.vg

vg deconstruct -p chr3:195721919-195989697 -e -A scaff aln.vg > aln.vcf
vg deconstruct -p chr3:195721919-195989697 -e -A scaff aln2.vg > aln2.vcf

# build the graph
cat aln.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > aln.sorted.vcf
bgzip aln.sorted.vcf ; tabix aln.sorted.vcf.gz
vg construct -a -r chr3_195721919_195989697.fasta -v aln.sorted.vcf.gz > final.vg
vg index -x final.xg -G final.gbwt -v aln.sorted.vcf.gz final.vg
(cat final.vg ; vg paths -g final.gbwt -x final.xg -T -V ) >x+.vg

vg paths -x final.xg -g final.gbwt -V -q scaff286_CF011 > aln_paths.vg

vg index -x final.xg -G final.gbwt -v aln.vcf -g final.gcsa final.vg
# extract the threads from the GBWT as vg Paths



#clean fastas
#------------------------------
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' consensus.fasta > uconsensus.fasta
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' hap1.fasta > uhap1.fasta
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' hap2.fasta > uhap2.fasta

#ref alignment
#------------------------------
#minimap2 -x asm20 -d ../../reference/hg38.mmi ../../reference/hg38.fa

minimap2 -cx asm20 --secondary=no ../../reference/hg38.mmi uconsensus.fasta > cref.paf
minimap2 -cx asm20 --secondary=no ../../reference/hg38.mmi uhap1.fasta > ref1.paf
minimap2 -cx asm20 --secondary=no ../../reference/hg38.mmi uhap2.fasta > ref2.paf

cat cref.paf ref1.paf ref2.paf > ref.paf

#graph alignment
#------------------------------

minimap2 -cx asm20 --secondary=no uconsensus.fasta uhap1.fasta > c1.paf
minimap2 -cx asm20 --secondary=no uconsensus.fasta uhap2.fasta > c2.paf
cat c1.paf c2.paf ref.paf > all.paf
rm c1.paf ; rm c2.paf

cat uconsensus.fasta uhap1.fasta uhap2.fasta > all.fasta

#build graph
#------------------------------
seqwish --repeat-max=1 -s all.fasta -p all.paf -g all.gfa

vg view -Fv all.gfa > all_unmod.vg
vg mod -X 256 all_unmod.vg > all.vg

vg stats -z all.vg
vg paths -L -v all.vg

#index and annotate
#------------------------------
vg index -x all.xg all.vg

BED=test.bed

vg annotate -p -x all.xg -b $BED > haps.gam
#vg mod -i test.gam all_mod.vg > all_mod_v2.vg
vg augment -B -i all.vg haps.gam > annotated.vg
vg stats -z annotated.vg
vg paths -L -v annotated.vg

#rm all.vg

#remove path
#------------------------------
vg paths -d -Q consensus -v annotated.vg > haps_temp.vg
vg mod -N haps_temp.vg > haps.vg
rm haps_temp.vg

vg stats -z haps.vg
vg paths -L -v haps.vg











cat uconsensus.fasta uhap1.fasta > 1c.fasta
cat uconsensus.fasta uhap2.fasta > 2c.fasta

minimap2 -cx asm20 --secondary=no uconsensus.fasta uhap1.fasta > 1c.paf
minimap2 -cx asm20 --secondary=no uconsensus.fasta uhap2.fasta > 2c.paf

seqwish --repeat-max=1 -s 1c.fasta -p 1c.paf -g 1c.gfa
seqwish --repeat-max=1 -s 2c.fasta -p 2c.paf -g 2c.gfa

gzip 1c.gfa
gzip 2c.gfa

zcat 1c.gfa.gz | vg view -Fv - > 1c.vg
zcat 2c.gfa.gz | vg view -Fv - > 2c.vg

vg mod -X 256 1c.vg > 1c_mod.vg
vg mod -X 256 2c.vg > 2c_mod.vg

echo ">hap1_subset" > part_hap1.fasta
tail -109779 uhap1.fasta | head -50000 >> part_hap1.fasta

vg index -x 1c_mod.xg 1c_mod.vg
vg index -g 1c_mod.gcsa 1c_mod.vg
vg index -x 2c_mod.xg 2c_mod.vg
vg index -g 2c_mod.gcsa 2c_mod.vg

vg map -x 1c_mod.xg -g 1c_mod.gcsa -F part_hap1.fasta > part.gam







#https://github.com/vgteam/vg/wiki/Index-Construction
vg index -x all_mod.xg all_mod.vg
vg prune all_mod.vg > all_mod.pruned.vg
mkdir -p tmp
vg index --temp-dir ./tmp -g all_mod.gcsa all_mod.pruned.vg
rm -f all_mod.pruned.vg


vg map -x all_mod.xg -g all_mod.gcsa -F part_hap1.fasta > part.gam

#remove path
#------------------------------
vg paths -d consensus -v all.vg > haps.vg
vg mod -N all_mod.vg





vg deconstruct -p consensus all_mod.vg -e hap1 > hap1_.vcf





















for d in ../../../slc9a3_polish/out/CF* ; do  
    cfid=`basename $d`
    echo $cfid
    
    cp $d/out/chr*67hapA.polished.fasta ${cfid}.hapA.fasta
    cp $d/out/chr*67hapB.polished.fasta ${cfid}.hapB.fasta
done


for d in * ; do  
    cfid=`basename $d`
    echo $cfid
    
    sed -i "s/>chr5_393462_677667|arrow|arrow|arrow|arrow|arrow|arrow/>${cfid}/" $d 
    
done




'''