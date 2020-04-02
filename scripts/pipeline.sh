
CFID=$1
BASEDIR="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid"
HG38="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/hg38.fa.gz"
HG38_ANNOTATIONS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/gene_annotations.gff"
QUEUE_NAME="tcagdenovo"
PYTHON="/hpf/largeprojects/struglis/scott/anaconda3/bin/python"

if [ -z "$QUEUE_NAME" ]; then QUEUE="" ; else QUEUE="-q $QUEUE_NAME"; fi

PURGE=${BASEDIR}/hybrid-pipeline/scripts/purge.sh
ALIGN_TIGS=${BASEDIR}/hybrid-pipeline/scripts/align_chunks/align_contigs.sh
#SUMMARIZE_HITS=${BASEDIR}/hybrid-pipeline/scripts/align_chunks/summarize_hits.sh

JOBOUT=${BASEDIR}/jobout/${CFID}
mkdir -p $JOBOUT

echo "STEP 1: PURGE DUPLICATES"
# =========================================

# EXPECTED RESULT OF STEP 1
REF_FA=${BASEDIR}/${CFID}/pacbio/purge/${CFID}.canu.purged.fa
QUERY_FA=${BASEDIR}/${CFID}/10x/purge/${CFID}.supernova.purged.fa

if [ -f "$REF_FA" ] && [ -f "$QUERY_FA" ]; then
   DEPEND_1=""
   echo "FOUND FILES: $REF_FA $QUERY_FA"
else
   PURGE_JID=$(bash $PURGE $CFID $BASEDIR $QUEUE_NAME | tail -1)
   DEPEND_1="-W depend=afterok:${PURGE_JID}"
fi

echo "STEP 2: ALIGN ASSEMBLIES & BUILD BLOCKS"
# =========================================
BLOCKDIR=${BASEDIR}/${CFID}/block
SUMMARY_PREFIX=${CFID}_supernova_vs_canu

# EXPECTED RESULT OF STEP 2
SUMMARY=${BASEDIR}.id${ID}.cov${COV}.summary.txt

if [ -f "$SUMMARY" ]; then
   DEPEND_2=""
else
   JOB="bash $ALIGN_TIGS $REF_FA $QUERY_FA $BLOCKDIR $SUMMARY_PREFIX $QUEUE_NAME"
   BLOCK_JID=$(echo $JOB | qsub $QUEUE $DEPEND_1 -l nodes=1:ppn=1 -l mem=16g -l vmem=16g -l walltime=6:00:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N align_tigs_${CFID} "-")

   DEPEND_2="-W depend=afterok:${BLOCK_JID}"
fi

echo "STEP 3: MAKE DRAFT HYBRID"
# =========================================

HYBRIDDIR=${BASEDIR}/${CFID}/hybrid
mkdir -p $HYBRIDDIR

JOB="$PYTHON ${BASEDIR}/hybrid-pipeline/src/main.py $SUMMARY $QUERY_FA $REF_FA -o $HYBRIDDIR"
HYBRID_JID=$(echo $JOB | qsub $QUEUE $DEPEND_2 -l nodes=1:ppn=1 -l mem=32g -l vmem=32g -l walltime=6:00:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N hybrid_${CFID} "-")


echo "STEP 4: RUN QUAST"
# =========================================

QUERY_FA_RAW=`echo ${BASEDIR}/${CFID}/10x/*pseudohap.fasta*`
REF_FA_RAW=`echo ${BASEDIR}/${CFID}/pacbio/*.contigs.fasta`

QUASTDIR=${BASEDIR}/${CFID}/quast
mkdir -p $QUASTDIR

JOB="module load quast ; \
     quast-lg.py $REF_FA_RAW $REF_FA $QUERY_FA_RAW $QUERY_FA $hybrid \
                 -o $QUASTDIR -r $HG38 -g $HG38_ANNOTATIONS"
	
#Options:
#-o  --output-dir  <dirname>       Directory to store all result files [default: quast_results/results_<datetime>]
#-r                <filename>      Reference genome file
#-g  --features [type:]<filename>  File with genomic feature coordinates in the reference (GFF, BED, NCBI or TXT)
#                                  Optional 'type' can be specified for extracting only a specific feature type from GFF
#-m  --min-contig  <int>           Lower threshold for contig length [default: 3000]
#-t  --threads     <int>           Maximum number of threads [default: 25% of CPUs]

QUAST_JID=$(echo $JOB | qsub $QUEUE -W depend=afterok:${HYBRID_JID} -l nodes=1:ppn=8 -l mem=64g -l vmem=64g -l walltime=71:59:00 -o $JOBOUT -e $JOBOUT -d `pwd` -N quast_${CFID} "-")

