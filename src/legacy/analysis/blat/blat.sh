#!/bin/bash

#PBS -l mem=16g
#PBS -l vmem=16g
#PBS -l gres=localhd:10
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=4
#PBS -o /hpf/largeprojects/struglis/angela/log
#PBS -e /hpf/largeprojects/struglis/angela/log

module load blat

date
cd $OUT_DIR
echo $PWD

blat -fastMap -noHead -t=dna -q=dna -out=pslx $REF $QUERY $PSLX

echo "Done running BLAT for"${QUERY}
