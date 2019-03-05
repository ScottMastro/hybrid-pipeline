#!/usr/bin/env bash

#note: requires gzipped and indexed fa files
canu=$1
supernova=$2

#gzip -k $canu
#gzip -k $supernova

module load samtools

samtools faidx $canu
samtools faidx $supernova
