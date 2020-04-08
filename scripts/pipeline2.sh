CFID=$1
BASEDIR="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid"
HG38="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/hg38.fa.gz"
HG38_ANNOTATIONS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/gene_annotations.gff"
QUEUE_NAME="tcagdenovo"
PYTHON="${BASEDIR}/tools/python"

if [ -z "$QUEUE_NAME" ]; then QUEUE="" ; else QUEUE="-q $QUEUE_NAME"; fi

LONGRANGER_ALN=${BASEDIR}/hybrid-pipeline/scripts/run_longranger.sh
PBMM2_ALN=${BASEDIR}/hybrid-pipeline/scripts/run_pbmm2.sh
POLISH_SCRIPT=${BASEDIR}/hybrid-pipeline/src/main_polish.sh

JOBOUT=${BASEDIR}/jobout/${CFID}
mkdir -p $JOBOUT

echo "STEP 5: REALIGNMENT"
# =========================================
HYBRIDDIR=${BASEDIR}/${CFID}/hybrid
HYBRID_FA=${HYBRIDDIR}/hybrid_assembly.fasta

REALIGN_DIR=${BASEDIR}/${CFID}/hybrid_align
LONGRANGER_DIR=$REALIGN_DIR/10x
mkdir -p $LONGRANGER_DIR
READSDIR_10X=${BASEDIR}/${CFID}/10x/reads
LONGRANGER=${BASEDIR}/tools/longranger
PICARD=${BASEDIR}/tools/picard.jar

PBMM2_DIR=${REALIGN_DIR}/pacbio
PBMM2_TEMP_DIR=${PBMM2_DIR}/alns
mkdir -p $PBMM2_TEMP_DIR
READSDIR_PB=${BASEDIR}/${CFID}/pacbio/reads
PBMM2=${BASEDIR}/tools/pbmm2

# EXPECTED RESULT OF STEP 5
PBMM2_BAM=${PBMM2_DIR}/${CFID}.pbmm2.bam
LONGRANGER_BAM=${LONGRANGER_DIR}/${CFID}.longranger.bam

if [ -f $LONGRANGER_BAM ]; then
   echo "ALIGNMENT FILE FOUND, SKIPPING: $LONGRANGER_BAM"
   LR_JID=""
else
   CORES=32 
   MEM=121
   JOB="bash $LONGRANGER_ALN $CFID $HYBRID_FA $READSDIR_10X $LONGRANGER $PICARD $LONGRANGER_DIR"
   LR_JID=$(echo $JOB | qsub $QUEUE -l nodes=1:ppn=$CORES -l mem=${MEM}g -l vmem=${MEM}g -l walltime=239:00:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N ${CFID}_longranger_aln "-")
   LR_JID=:$LR_JID
fi

if [ -f $PBMM2_BAM ]; then
   echo "ALIGNMENT FILE FOUND, SKIPPING: $PBMM2_BAM"
   PBMM2_JID=""
else
   BAMLIST=${PBMM2_DIR}/bamlist.txt
   echo -n > $BAMLIST

   JID_LIST=""

   for READS in $READSDIR_PB/*; do
      CORES=16
      MEM=64

      JOB="bash $PBMM2_ALN $HYBRID_FA $READS $PBMM2 $PBMM2_TEMP_DIR"
      JID=$(echo $JOB | qsub $QUEUE -l nodes=1:ppn=$CORES -l mem=${MEM}g -l vmem=${MEM}g -l walltime=16:00:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N ${CFID}_pbmm2_aln "-")
      JID_LIST=${JID_LIST}:${JID}
   
      echo ${PBMM2_TEMP_DIR}/${READS##*/}.pbmm2.bam >> $BAMLIST
   done

   JOB="module load samtools/1.9 ; samtools merge --threads 16 -b $BAMLIST $PBMM2_BAM $BAMLIST ; samtools index $PBMM2_BAM"
   PBMM2_JID=$(echo $JOB | qsub $QUEUE -W depend=afterok${JID_LIST} -l nodes=1:ppn=16 -l mem=24g -l vmem=24g -l walltime=41:59:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N ${CFID}_pbmm2_merge "-")
   PBMM2_JID=:$PBMM2_JID
fi

if [ ! -z $LR_JID ] || [ ! -z $PBMM2_JID ] ; then
   DEPEND_4="-W depend=afterok${LR_JID}${PBMM2_JID}"
else
   DEPEND_4=""
fi

echo "STEP 6: POLISHING"
# =========================================



