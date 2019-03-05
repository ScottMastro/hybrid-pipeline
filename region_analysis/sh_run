#!/bin/bash

#PBS -l mem=64g
#PBS -l vmem=64g
#PBS -l gres=localhd:10
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=16
#PBS -o /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/region_analysis/jobout/
#PBS -e /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/region_analysis/jobout/
#PBS -N region_analysis

cd /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/region_analysis

#name=$1
#region=$2

#cftr = chr7:117469963-117678664
#slc9a3 = chr5:461986-534434
#slc26a9 = chr1:205903000-206003600

name=cftr
chr=7
start=117469963
end=117678664
region=chr${chr}:${start}-${end}
minlen=2000
minid=95

module load samtools
module load mummer
module load fig2dev
module load bedops
module load bedtools

mkdir -p regions

fasta=./regions/$name.ref.fa
samtools faidx ./hg38/hg38.fa.gz $region > $fasta
sed -i "1s/.*/\>$name $region/" $fasta

echo -e "${chr}\t${start}\t${end}" > $name.region.bed
bedtools intersect -a ./hg38/gene_annotations.gtf -b $name.region.bed | 
 awk -v r="$name" -F $'\t' 'BEGIN {OFS = FS}{$1=r}1' > $name.temp.gff
bedtools shift -s -$start -i $name.temp.gff -g /dev/null > ./regions/$name.gff
rm $name.region.bed
rm $name.temp.gff

mkdir -p out

for filename in ./assemblies/*.fasta; do	
    base=${filename##*/}.$name
    echo $base
    nucmer -t 16 -p ./out/$base ./regions/$name.ref.fa $filename
    delta-filter -i $minid -l $minlen ./out/$base.delta > ./out/$base.filtered.delta
    show-coords -B -r ./out/$base.filtered.delta > ./out/$base.coords
    show-snps -TH ./out/$base.filtered.delta > ./out/$base.snps
    #/hpf/tools/centos6/mummer/3.23/mapview -I -Ir -n 1 -p ./out/$base -f pdf ./out/$base.coords ./regions/$name.gff
done

