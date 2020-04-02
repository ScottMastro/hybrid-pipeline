#!/bin/bash

DATA=$1
XREADDIR="/hpf/largeprojects/tcagstor/projects_snapshort/illumina/basecalling_10x"

while read LINE; do
   #CFIT_ID	PB_SEQUENCE_ID	GMS_ID	10X_SEQUENCE_ID	NOTES	10X_RUN_ID	LONGRANGER	SUPERNOVA	CANU				

   CFID=`echo "${LINE}" | cut -f 1`
   PBID=`echo "${LINE}" | cut -f 2`
   GMSID=`echo "${LINE}" | cut -f 3`
   XID=`echo "${LINE}" | cut -f 4`
   NOTES=`echo "${LINE}" | cut -f 5`
   RUNID=`echo "${LINE}" | cut -f 6`
   LRDIR=`echo "${LINE}" | cut -f 7`
   SNDIR=`echo "${LINE}" | cut -f 8`
   PBDIR=`echo "${LINE}" | cut -f 9`
   CNDIR=`echo "${LINE}" | cut -f 10`
   
   echo "-------------- $CFID --------------"
   
   if [ -z "$PBID" ]; then
	  echo "Skipping..."
      continue
   fi
   if [ -z "$SNDIR" ]; then
	  echo "Skipping..."
      continue
   fi
   if [ -z "$CNDIR" ]; then
	  echo "Skipping..."
      continue
   fi
   if [ -z "$RUNID" ]; then
	  echo "Skipping..."
      continue
   fi

   
   mkdir -p ${CFID}/pacbio/reads
   
   # link 10x reads
   #---------------------------------------------
   export IFS=","
   for RUN in $RUNID; do
      echo "Linking reads from ${RUN}"
	  READS=`echo ${XREADDIR}/${RUN}*bc/output/${RUN}/${XID}`
	  
      if [ ! -d "$READS" ]; then
		 echo "READ DIRECTORY DOES NOT EXIST!! $READS"
		 exit 1
	  else
	     mkdir -p ${CFID}/10x/reads
	     ln -s ${READS}/*.f*q.gz ${CFID}/10x/reads/
	  fi
	  
   done
   
   # link 10x assembly
   #---------------------------------------------
   echo "Linking 10x assembly"
   ASSEMBLY=`echo ${SNDIR}/${XID}_supernova_mkoutput/*.pseudohap.f*a.gz`

   if [ ! -f "$ASSEMBLY" ]; then
      echo "ASSEMBLY DOES NOT EXIST!! $ASSEMBLY"
	  exit 1

   else
      mkdir -p ${CFID}/10x
	  ln -s ${ASSEMBLY} ${CFID}/10x

   fi
  
   # link PacBio reads
   #---------------------------------------------
   echo "Linking PacBio reads"
   mkdir -p ${CFID}/pacbio/reads
   ln -s ${PBDIR}/*/*/*.subreads.bam ${CFID}/pacbio/reads/

   
   # link PacBio assembly and unitigs
   #---------------------------------------------
   echo "Linking PacBio assembly"
   ASSEMBLY=`echo ${CNDIR}/*.contigs.f*a`

   if [ ! -f "$ASSEMBLY" ]; then
      echo "ASSEMBLY DOES NOT EXIST!! $ASSEMBLY"
	  exit 1
   else
      mkdir -p ${CFID}/pacbio
	  ln -s ${ASSEMBLY} ${CFID}/pacbio
   fi
   
   UNITIGS=`echo ${CNDIR}/*.unitigs.bed`

   if [ ! -f "$UNITIGS" ]; then
      echo "UNITIGS BED FILE DOES NOT EXIST!! $UNITIGS"
	  exit 1
   else
      mkdir -p ${CFID}/pacbio
	  ln -s ${UNITIGS} ${CFID}/pacbio
   fi

done < $DATA
