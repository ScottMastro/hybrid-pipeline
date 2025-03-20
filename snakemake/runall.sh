ID=${1:-RUN_ALL}
BASEDIR=$2
TXASSEM=$3
PBASSEM=$4
PBUNITIGS=$5

SNAKEFILE_DIR=1_snakemake/

# ---- RUN SNAKEMAKE IF SINGLE SAMPLE PROVIDED, ELSE LOOP THROUGH ALL SAMPLES AND RUN THIS SCRIPT AS QUEUED JOB ----

if [ ! "$ID" = "RUN_ALL" ]; then
  echo "RUNNING SAMPLE $ID"

  #load correct python version with snakemake
  DIR_=/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/anaconda
  __conda_setup="$(${DIR_}'/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
  if [ $? -eq 0 ]; then
    eval "$__conda_setup"
  else
    if [ -f "${DIR_}/etc/profile.d/conda.sh" ]; then
    . "${DIR_}/etc/profile.d/conda.sh"
    else
        export PATH="${DIR_}/bin:$PATH"
    fi
  fi
  unset __conda_setup
  echo "Active python version: " `which python`
  
  mkdir -p ./jobout

  DIR=${BASEDIR}/${ID}

  HYBRIDFA=${DIR}/hybrid/hybrid_assembly.fasta.gz

  check=$(ls ${HYBRIDFA} | wc -l)
  if [ ! "$check" = "1" ]; then
     
     ENV_FILE=${SNAKEFILE_DIR}set_env.sh
     #snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_hybrid --config env=${ENV_FILE} r=$PBASSEM q=$TXASSEM out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS scripts=${BASEDIR}/../hybrid-pipeline/ --unlock
     #snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_hybrid --config env=${ENV_FILE} r=$PBASSEM q=$TXASSEM out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS  \
     #   scripts=${BASEDIR}/../hybrid-pipeline/  --cluster "sbatch -N 1 --tmp 100G --mem {resources.mem}G -t 52:00:00 -c {threads} -o ./jobout/slurm-%j.out -e ./jobout/slurm-%j.out --job-name ${ID}_hybrid_task"

     echo """snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_hybrid --config env=${ENV_FILE} r=$PBASSEM q=$TXASSEM out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS  \
        scripts=${BASEDIR}/../hybrid-pipeline/  --cluster "sbatch -N 1 --tmp 100G --mem {resources.mem}G -t 52:00:00 -c {threads} -o ./jobout/slurm-%j.out -e ./jobout/slurm-%j.out --job-name ${ID}_hybrid_task""""


     JOBOUT=./jobout
     JOBSH=./jobsh
     mkdir -p $JOBOUT
     mkdir -p $JOBSH

     LOG=${JOBOUT}/${ID}_hybrid-%j.out

     # unlock
     snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_hybrid --config env=${ENV_FILE} r=$PBASSEM q=$TXASSEM out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS scripts=${BASEDIR}/../hybrid-pipeline/ --unlock

     JOB="snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_hybrid --config env=${ENV_FILE} r=$PBASSEM q=$TXASSEM out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS  \
       scripts=${BASEDIR}/../hybrid-pipeline/ --rerun-incomplete --cluster \"sbatch -N 1 --tmp 100G --mem {resources.mem}G -t 52:00:00 -c {threads} -o ./jobout/slurm-%j.out -e ./jobout/slurm-%j.out --job-name ${ID}_hybrid_task\""

     echo $JOB

     echo "#!/bin/bash" > ${JOBSH}/job.${ID}.hybrid.sh
     echo $JOB >> ${JOBSH}/job.${ID}.hybrid.sh

     sbatch -N 1 --mem 4G -t 23:00:00 -o $LOG -e $LOG -D $PWD --job-name ${ID}_hybrid ${JOBSH}/job.${ID}.hybrid.sh

     exit 1
fi

  #--cluster "qsub -q tcagdenovo -l nodes=1:ppn={threads} -l mem={resources.mem}g -l walltime=16:00:00 -o ./jobout -e ./jobout"

  echo "---"

  #ENV_FILE=${SNAKEFILE_DIR}set_env_stats.sh
  #echo snakemake --jobs 20 -s ${SNAKEFILE_DIR}Snakefile_stats --config env=${ENV_FILE} out=${DIR}/stats hybrid=${DIR}/hybrid/ r=$PBASSEM q=$TXASSEM --unlock
  #echo snakemake --jobs 20 -s ${SNAKEFILE_DIR}Snakefile_stats --config env=${ENV_FILE} out=${DIR}/stats hybrid=${DIR}/hybrid/ r=$PBASSEM q=$TXASSEM  \
  #   --cluster "qsub -q tcagdenovo -l nodes=1:ppn={threads} -l mem={resources.mem}g -l walltime=82:00:00 -o ./jobout -e ./jobout"

  #echo "---"

  #ENV_SDA_FILE=${SNAKEFILE_DIR}set_env_sda.sh
  #ENV_FILE=${SNAKEFILE_DIR}set_env.sh
  #snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_prep --config env=${ENV_FILE} env_sda=$ENV_SDA_FILE hybrid=${DIR}/hybrid rreads=${DIR}/pacbio/reads qreads=${DIR}/10x/reads --unlock

  #echo snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_prep --config env=${ENV_FILE} env_sda=$ENV_SDA_FILE hybrid=${DIR}/hybrid rreads=${DIR}/pacbio/reads qreads=${DIR}/10x/reads  \
  #   --cluster "qsub -q tcagdenovo -l nodes=1:ppn={threads} -l mem={resources.mem}g -l walltime=52:00:00 -o ./jobout -e ./jobout"
  #snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_prep --config env=${ENV_FILE} env_sda=$ENV_SDA_FILE hybrid=${DIR}/hybrid rreads=${DIR}/pacbio/reads qreads=${DIR}/10x/reads  \
  #   --cluster "sbatch -N 1 --mem {resources.mem}G -t 52:00:00 -c {threads} -o ./jobout/slurm-%j.out -e ./jobout/slurm-%j.out --job-name ${ID}_${TIG}_polish"

  #echo "---"
  
  JOBOUT=./jobout
  mkdir -p $JOBOUT
  LOG=${JOBOUT}/${ID}_prep-%j.out
  ENV_SDA_FILE=${SNAKEFILE_DIR}set_env_sda.sh
  ENV_FILE=${SNAKEFILE_DIR}set_env.sh
  JOBSH=./jobsh
  mkdir -p $JOBSH
    
  PBALN=${DIR}/hybrid/hybrid_assembly.pbmm2.bam
  TXALN=${DIR}/hybrid/hybrid_assembly.longranger.bam
  SDADONE=${DIR}/hybrid/sda/sda.done

  check1=$(ls ${PBALN} | wc -l)
  check2=$(ls ${TXALN} | wc -l)
  check3=$(ls ${SDADONE} | wc -l)

  if [ "$check1 $check2 $check3" = "1 1 1" ]; then
      echo "Prep done."

  else
  
     SCRIPT=${JOBSH}/job.${ID}.prep.sh

     UNLOCK="""snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_prep --config env=${ENV_FILE} env_sda=$ENV_SDA_FILE hybrid=${DIR}/hybrid rreads=${DIR}/pacbio/reads qreads=${DIR}/10x/reads --unlock"""
     JOB="""snakemake --jobs 55 -s ${SNAKEFILE_DIR}Snakefile_prep --config env=${ENV_FILE} env_sda=$ENV_SDA_FILE hybrid=${DIR}/hybrid rreads=${DIR}/pacbio/reads qreads=${DIR}/10x/reads  \
        --cluster \"sbatch -N 1 --tmp 1000G --mem {resources.mem}G -t 52:00:00 -c {threads} -o ./jobout/slurm-%j.out -e ./jobout/slurm-%j.out --job-name ${ID}_prep_task\""""

     echo "---"
  
     echo $JOB
     echo "#!/bin/bash" > $SCRIPT
     echo $UNLOCK >> $SCRIPT
     echo $JOB >> $SCRIPT

     sbatch -N 1 --mem 4G -t 123:00:00 -o $LOG -e $LOG -D $PWD --job-name ${ID}_prep $SCRIPT
     exit 1

  fi

  JOB="source $ENV_FILE ; snakemake --jobs 3 -s ${SNAKEFILE_DIR}Snakefile_polish --config sample=${ID} env=${ENV_FILE} out=${DIR}/hybrid/polish hybridfa=$HYBRIDFA sda=${DIR}/hybrid/sda r=$PBALN q=$TXALN iters=2 tig=scaff388 \
      --latency-wait 20 --cluster \"qsub -q tcagdenovo -l nodes=1:ppn={threads} -l mem={resources.mem}g -l walltime=52:00:00 -o $JOBOUT -e $JOBOUT\""
  echo $JOB

  JOBS="./jobs.txt"
  squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me > $JOBS
  QUEUE="-q tcagdenovo"


  JOBS="./jobs.txt"
  squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me > $JOBS


  cat ${DIR}/hybrid/${ID}.scaff.t2t.txt | grep "chr7" | cut -f1 | while read -r TIG ; do

    if [ -f "${DIR}/hybrid/polish/${TIG}/.done_compress" ]; then
      echo "$TIG done."
      continue

    else

      RUNNING=`cat $JOBS | grep ${ID}_${TIG}_ | wc -l`

      if [ "$RUNNING" -ne "0" ]; then

        echo "$TIG running."
        continue

      fi


      echo "doing $TIG"

      JOBOUT=${DIR}/hybrid/polish/${TIG}/jobout
      mkdir -p ${DIR}/hybrid/polish/${TIG}
      mkdir -p $JOBOUT

      ENV_FILE=${SNAKEFILE_DIR}set_env.sh
      LOG=${JOBOUT}/${ID}_${TIG}_task-%j.out

      # unlock
      source $ENV_FILE ; snakemake --jobs 3 -s ${SNAKEFILE_DIR}Snakefile_polish --config sample=${ID} env=${ENV_FILE} out=${DIR}/hybrid/polish tig=$TIG \
         hybridfa=$HYBRIDFA sda=${DIR}/hybrid/sda r=$PBALN q=$TXALN iters=2 --latency-wait 120 --unlock

      #snakemake --jobs 3 -s ${SNAKEFILE_DIR}Snakefile_polish --config sample=${ID} env=${ENV_FILE} out=${DIR}/hybrid/polish tig=$TIG \
      #   hybridfa=$HYBRIDFA sda=${DIR}/hybrid/sda r=$PBALN q=$TXALN iters=2 --touch

      JOB="""snakemake --jobs 3 -s ${SNAKEFILE_DIR}Snakefile_polish --config sample=${ID} env=${ENV_FILE} out=${DIR}/hybrid/polish tig=$TIG \
         hybridfa=$HYBRIDFA sda=${DIR}/hybrid/sda r=$PBALN q=$TXALN iters=2 --latency-wait 120 \
         --cluster \"sbatch -N 1 -c {threads} --tmp 100G --mem {resources.mem}G -t 52:00:00 -o $LOG -e $LOG\""

      echo $JOB

      JOBSH=./jobsh
      mkdir -p $JOBSH
      echo "#!/bin/bash" > ${JOBSH}/job.${ID}.${TIG}.sh
      echo $JOB >> ${JOBSH}/job.${ID}.${TIG}.sh

      sbatch -N 1 --mem 4G -t 23:00:00 -o $LOG -e $LOG -D $PWD --job-name ${ID}_${TIG}_polish ${JOBSH}/job.${ID}.${TIG}.sh

    fi

  done

else

  BASEDIR=`pwd`/data
  for d in $BASEDIR/*; do
    ID=`basename $d`

    echo "#------"
    echo '#'$ID

    TXSUFF="_supernova.pseudohap.fasta.gz.irods"
    TXID=`basename ${d}/10x/*$TXSUFF`
    TXID=${TXID%"$TXSUFF"}
  
    TXASSEM=${d}/10x/${TXID}${TXSUFF}
    TXASSEM=${TXASSEM%".irods"}

    PBSUFF=".contigs.fasta.irods"
    PBID=`basename ${d}/pacbio/*$PBSUFF`
    PBID=${PBID%"$PBSUFF"}

    PBASSEM=${d}/pacbio/${PBID}${PBSUFF}
    PBASSEM=${PBASSEM%".irods"}

    PBUNITIGS=${d}/pacbio/${PBID}.unitigs.bed


    JOB="bash runall.sh $ID $BASEDIR $TXASSEM $PBASSEM $PBUNITIGS"
    echo $JOB 

  done
fi


