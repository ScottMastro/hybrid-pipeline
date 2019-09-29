#!/bin/bash

#PBS -l mem=16g
#PBS -l vmem=16g
#PBS -l gres=localhd:10
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -o /hpf/largeprojects/struglis/angela/log
#PBS -e /hpf/largeprojects/struglis/angela/log

module load python/3.7.1

# imported variables: 
OUT_DIR=$OUT_DIR
WD=$WD
REF=$REF
CONTIG_NAME=$CONTIG_NAME
N=$N


cd $OUT_DIR
date
echo $PWD

set -o x trace

PSLX=$OUT_DIR"/all.pslx"

# concatenate all the pslx file into one file
cat *.pslx > all.pslx

echo "Concatenated all psl files into all.pslx"

# parse the pslx file to select the top hit for each "chop"
# output summary file, filtered pslx with only the top hit, and a vep input file
$WD"/pslx_parse.py" -pslx $PSLX -n $N

HITS="../all.pslx_tophits.pslx"
VEP="./vep_input.txt"

# remove the first blank line
sed -i '1d' $VEP 

# VEP SNPs from top BLAT hits for each contig
qsub -v DIR=$OUT_DIR,VEP=$VEP,CONTIG_NAME=$CONTIG_NAME \
-N "vep_BLAT_hits:"$CONTIG_NAME \
$WD"/vep_rsid.sh"


# create a new directory with all the contig plots
mkdir -p contig_images
cd contig_images

# plot the contigs and their hits to different chromosomes 
$WD"/plot_contig.py" -hits $HITS


echo "Mapping contig onto the genome complete"

	
