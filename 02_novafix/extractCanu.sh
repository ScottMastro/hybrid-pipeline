#!/bin/bash

dir=$1
contig=$2
left=$3
right=$4
len="$(($right-$left+1))"

#extract region, remove header, remove newlines
tig="$(samtools faidx ${dir} ${contig}_pilon_pilon:${left}-${right} | tail -n +2 | tr -d '\n')"

echo $tig
