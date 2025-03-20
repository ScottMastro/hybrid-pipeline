module load samtools/1.9

FASTA=$1
READS=$2
PBMM2=$3
OUTDIR="${4:-"."}"
CORES="${5:-"16"}"
MEM="${6:-"64"}"

CWD=`pwd`
cd $OUTDIR

echo "Running pbmm2"
# =========================================
FILENAME=${READS##*/}
echo $FILENAME

BAM=${OUTDIR}/${FILENAME}.pbmm2.bam
$PBMM2 align -j $CORES --sort $FASTA $READS > $BAM

echo "Reheader BAM to be compliant with PB BAM standards"
# =========================================

REHEADER_BAM=${OUTDIR}/${FILENAME}.reheader.bam
samtools view -H $BAM | sed "s/SO:coordinate/pb:3.0.4\tSO:coordinate/" | samtools reheader - $BAM > $REHEADER_BAM
mv $REHEADER_BAM $BAM

samtools index $BAM
cd $CWD
