module load java/1.8.0_91
module load samtools/1.9

SAMPLENAME=$1
FASTA=$2
READSDIR=$3
LONGRANGER=$4
PICARD=$5

OUTDIR="${6:-"."}"
CORES="${7:-"32"}"
MEM="${8:-"121"}"

CWD=`pwd`
cd $OUTDIR
OUTDIR=`pwd`

echo "Creating indexes"
# =========================================

FILENAME=`basename $FASTA`
PREFIX=${FILENAME%%.*}

REF_INDEX=${OUTDIR}/refdata-${PREFIX}
echo $REF_INDEX

$LONGRANGER mkref $FASTA 
java -jar $PICARD CreateSequenceDictionary R=${REF_INDEX}/fasta/genome.fa O=${REF_INDEX}/fasta/genome.dict

echo "Running Long Ranger"
# =========================================

$LONGRANGER align --id=${SAMPLENAME} --fastqs=`ls -d ${READSDIR}` \
 --reference=${REF_INDEX} --jobmode=local --localcores=${CORES} --localmem=${MEM}

mv ${SAMPLENAME}/outs/possorted_bam.bam ${SAMPLENAME}.longranger.bam
mv ${SAMPLENAME}/outs/possorted_bam.bam.bai ${SAMPLENAME}.longranger.bam.bai

cd $CWD
