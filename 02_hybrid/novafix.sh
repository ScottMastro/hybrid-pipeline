#!/usr/bin/env bash

#note: requires indexed fa files
canu=$1
supernova=$2

input=$3
: ${input:=../blocks}

output=$4
: ${output:=./output}

mkdir -p $output

block_file=$input/blocks.txt
summary_file=$input/summary.txt

echo "Starting novafix"

module load python/2.7.9
echo "construct blocks"
python python/novafix.py $block_file $summary_file $canu $supernova $output

