#!/bin/bash

contig=$1
left=$2
right=$3

tig="$(sed -n "/>${contig} edges*/,/>/{/>${contig}/b;/>/b;p}" CF003B2D_318510_supernova.pseudohap2.1.fasta)"

echo $(echo ${tig:$left:$right} | tr -d ' ')

