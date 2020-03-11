#!/bin/bash

#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -l gres=localhd:10
#PBS -l walltime=1:59:00
#PBS -l nodes=1:ppn=2
#PBS -o /hpf/largeprojects/struglis/angela/log
#PBS -e /hpf/largeprojects/struglis/angela/log
#PBS -N contig_segmentation

module load python/3.7.1

WD=$WD #where the scripts are located
REF=$REF
CONTIG=$CONTIG

# checking if the N option was specified 
if [ -z "$N" ]
then
	echo "Option N missing, default N=1000"
    N=1000
else
    N=$N
fi

REF="/hpf/largeprojects/struglis/angela/lib/hg38/GRCh38.main_chr.fa"


# name of the contig file 
BASE_CONTIG_NAME="$(basename -- $CONTIG)"
CONTIG_NAME="${BASE_CONTIG_NAME%.*}"

# move to working directory
OUT="/hpf/largeprojects/struglis/angela/blat_out"
mkdir -p $OUT
OUT_DIR=$OUT"/blat_output_"$CONTIG_NAME"_"$N 
mkdir -p $OUT_DIR
date

set -o x trace

# chop up the contigs into random N segments with 10% coverage (assuming no overlap)
# save "chops" in new fasta files in OUTDIR
$WD"/contig_chop_random.py" -contig $CONTIG -n $N -outdir $OUT_DIR

cd $OUT_DIR


job_ids=""
# BLAT each "chop" file
for chops in *_chops*.fa; do
chop="${chops%.*}"
job_id=$(qsub -v OUT_DIR=$OUT_DIR,REF=$REF,QUERY=$chops,PSLX=${chop}".pslx" \
-N BLAT_contigs_${chop} \
$WD"/blat.sh")
job_ids="$job_ids:$job_id"
done


job_ids="${job_ids:1}"
echo "BLAT Job IDs:"
echo $job_ids

# Job on hold until all BLAT processes are complete
# Then concatenate all the output pslx files 
qsub -W depend=afterok:$job_ids -v \
WD=$WD,OUT_DIR=$OUT_DIR,REF=$REF,N=$N,CONTIG_NAME=$CONTIG_NAME \
-N "concatenate_pslx:"$CONTIG_NAME \
$WD"/concatenate_pslx.sh"

echo "Main script complete"
