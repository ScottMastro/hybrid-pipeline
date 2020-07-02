CFID=$1
BASEDIR="${2:-"."}"
QUEUE_NAME="${3:-""}"
if [ -z "$QUEUE_NAME" ]; then QUEUE="" ; else QUEUE="-q $QUEUE_NAME"; fi

PURGE_DUPS=${BASEDIR}/tools/purge_dups
READDIR=${BASEDIR}/${CFID}/pacbio/reads
ASSEM=`echo ${BASEDIR}/${CFID}/pacbio/*.contigs.fasta`
OUTDIR=${BASEDIR}/${CFID}/pacbio/purge

FQDIR="${OUTDIR}/fa_reads"
ALNDIR="${OUTDIR}/alignments"
READSLIST=${OUTDIR}/reads.txt

mkdir -p $OUTDIR
mkdir -p ${OUTDIR}/jobout
mkdir -p $ALNDIR
mkdir -p $FQDIR

XOUTDIR=${BASEDIR}/${CFID}/10x/purge
XASSEM=`echo ${BASEDIR}/${CFID}/10x/*pseudohap.fasta*`

XALNDIR="${XOUTDIR}/alignments"

mkdir -p $XOUTDIR
mkdir -p ${XOUTDIR}/jobout
mkdir -p $XALNDIR

XASSEM_PREFIX=`basename $XASSEM`
XSPLIT="${XOUTDIR}/${XASSEM_PREFIX}.split"
XSELFALN="${XALNDIR}/${XASSEM_PREFIX}.split.self.paf.gz"

CWD=`pwd`
cd $OUTDIR

ASSEM_PREFIX=`basename $ASSEM`
SPLIT="${OUTDIR}/${ASSEM_PREFIX}.split"
SELFALN="${ALNDIR}/${ASSEM_PREFIX}.split.self.paf.gz"

#---------------------------------------------------

echo "Set up BAM->FASTQ jobs"

JID_LIST_1=""
for BAM in ${READDIR}/*bam; do

   NAME=$(basename $BAM)
   PREFIX="${NAME%%.bam}"
   echo $PREFIX

   FQ=${FQDIR}/${PREFIX}.fastq
   FASTQ=`echo "$(cd "$(dirname "$FQ")"; pwd)/$(basename "$FQ")"`
   JOB="module load samtools/1.9 ; samtools bam2fq $BAM > $FASTQ"
   JID=$(echo $JOB | \
         qsub $QUEUE -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -l walltime=11:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N bam2fastq_${PREFIX} "-")
   JID_LIST_1="${JID_LIST_1}:${JID}"
   echo $FASTQ >> $READSLIST
done

#---------------------------------------------------

echo "Set up read alignment jobs"

JID_LIST_2=""
while read pb
   do
      PREFIX=`basename $pb`
      echo "aligning $PREFIX"

          JOB="module load minimap2 ; minimap2 -xmap-pb $ASSEM $pb | gzip -c - > ${ALNDIR}/${PREFIX}.paf.gz"
          JID=$(echo $JOB | \
                qsub $QUEUE -W depend=afterok${JID_LIST_1} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=23:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_align_${CFID}_${PREFIX} "-")
          JID_LIST_2=${JID_LIST_2}":${JID}"

done < $READSLIST

#---------------------------------------------------

JOB="${PURGE_DUPS}/bin/pbcstat ${ALNDIR}/*.paf.gz ; ${PURGE_DUPS}/bin/calcuts ${OUTDIR}/PB.stat > ${OUTDIR}/cutoffs 2> ${OUTDIR}/calcults.log"
JID_1=$(echo $JOB | qsub $QUEUE -W depend=afterok${JID_LIST_2} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=3:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N pbcstat_${CFID} "-")
echo "Produces PB.base.cov and PB.stat files"
echo $JOB

JOB="module load minimap2 ; ${PURGE_DUPS}/bin/split_fa $ASSEM > $SPLIT ; minimap2 -xasm5 -DP $SPLIT $SPLIT | gzip -c - > $SELFALN"
JID_2=$(echo $JOB | qsub $QUEUE -W depend=afterok:${JID_1} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_split_${CFID} "-")
echo "Split assembly and do a self-self alignment"
echo $JOB

JOB="${PURGE_DUPS}/bin/purge_dups -2 -T ${OUTDIR}/cutoffs -c ${OUTDIR}/PB.base.cov ${SELFALN} > ${OUTDIR}/dups.bed 2> ${OUTDIR}/purge_dups.log"
JID_3=$(echo $JOB | qsub $QUEUE -W depend=afterok:${JID_2} -l nodes=1:ppn=1 -l mem=64g -l vmem=64g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_dups_${CFID} "-")
echo "Purge haplotigs and overlaps"
echo $JOB

JOB="${PURGE_DUPS}/bin/get_seqs ${OUTDIR}/dups.bed $ASSEM"
JID_4=$(echo $JOB | qsub $QUEUE -W depend=afterok:${JID_3} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_out_${CFID} "-")
echo "Get purged primary and haplotig sequences"
echo $JOB


cd $XOUTDIR

#---------------------------------------------------

echo "Set up read alignment jobs - 10x"

JID_LIST_3=""
while read pb
   do
      PREFIX=`basename $pb`
      echo "aligning $PREFIX"

          JOB="module load minimap2 ; minimap2 -xmap-pb $XASSEM $pb | gzip -c - > ${XALNDIR}/${PREFIX}.paf.gz"
          JID=$(echo $JOB | \
                qsub $QUEUE -W depend=afterok${JID_LIST_1} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=23:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N xpurge_align_${CFID}_${PREFIX} "-")
          JID_LIST_3=${JID_LIST_3}":${JID}"

done < $READSLIST

#---------------------------------------------------

JOB="${PURGE_DUPS}/bin/pbcstat ${XALNDIR}/*.paf.gz ; ${PURGE_DUPS}/bin/calcuts ${XOUTDIR}/PB.stat > ${XOUTDIR}/cutoffs 2> ${XOUTDIR}/calcults.log"
XJID_1=$(echo $JOB | qsub $QUEUE -W depend=afterok${JID_LIST_3} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=3:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N xpbcstat_${CFID} "-")
echo "Produces PB.base.cov and PB.stat files"
echo $JOB

JOB="module load minimap2 ; ${PURGE_DUPS}/bin/split_fa $XASSEM > $XSPLIT ; minimap2 -xasm5 -DP $XSPLIT $XSPLIT | gzip -c - > $XSELFALN"
XJID_2=$(echo $JOB | qsub $QUEUE -W depend=afterok:${XJID_1} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N xpurge_split_${CFID} "-")
echo "Split assembly and do a self-self alignment"
echo $JOB

JOB="${PURGE_DUPS}/bin/purge_dups -2 -T ${XOUTDIR}/cutoffs -c ${XOUTDIR}/PB.base.cov ${XSELFALN} > ${XOUTDIR}/dups.bed 2> ${XOUTDIR}/purge_dups.log"
XJID_3=$(echo $JOB | qsub $QUEUE -W depend=afterok:${XJID_2} -l nodes=1:ppn=1 -l mem=64g -l vmem=64g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N xpurge_dups_${CFID} "-")
echo "Purge haplotigs and overlaps"
echo $JOB

JOB="${PURGE_DUPS}/bin/get_seqs ${XOUTDIR}/dups.bed $XASSEM"
XJID_4=$(echo $JOB | qsub $QUEUE -W depend=afterok:${XJID_3} -l nodes=1:ppn=1 -l mem=36g -l vmem=36g -l walltime=6:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N xpurge_out_${CFID} "-")
echo "Get purged primary and haplotig sequences"
echo $JOB

#---------------------------------------------------


JOB="rm -r $ALNDIR ; rm -r $XALNDIR ; rm $SPLIT ; rm $XSPLIT ;\
     mv ${OUTDIR}/purged.fa ${OUTDIR}/${CFID}.canu.purged.fa ; \
     mv ${OUTDIR}/hap.fa ${OUTDIR}/${CFID}.canu.hap.fa ; \
     mv ${XOUTDIR}/purged.fa ${XOUTDIR}/${CFID}.supernova.purged.fa ; \
     mv ${XOUTDIR}/hap.fa ${XOUTDIR}/${CFID}.supernova.hap.fa ; \
	 rm -r $READSLIST ; rm -r $FQDIR"
FINAL_JID=$(echo $JOB | qsub $QUEUE -W depend=afterany:${JID_4}:${XJID_4} -l nodes=1:ppn=1 -l mem=4g -l vmem=4g -l walltime=1:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_clean_${CFID} "-")
echo "Clean up job"
echo $JOB

cd $CWD
echo $FINAL_JID

#cd ${PURGE_DUPS}/runner && python setup.py install --user
#cd $CWD

#JOB="module load python/3.3.5 ; ${PURGE_DUPS}/scripts/pd_config.py ${ASSEM} ${OUTDIR}/reads.txt"
#JID_1=$(echo $JOB | qsub $QUEUE -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -l walltime=1:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_setup_${PREFIX} "-")
#echo "Set up purge config"


#run="module load python/3.3.5 ; module load minimap2 ; ${PURGE_DUPS}/scripts/run_purge_dups.py config.json ${PURGE_DUPS}/bin/ human -p bash"
#mkdir -p jobout
#echo $run |  qsub -q tcagdenovo -l nodes=1:ppn=6 -l mem=224g -l vmem=224g -l walltime=42:36:00 -o ./jobout/ -e ./jobout/ -d `pwd` -N purge_dups_${CFID} "-"

