#!/bin/bash

FASTA=$1
PAF=$2
OUT_PREFIX=$3

SEQWISH="/home/scott/bin/seqwish/bin/seqwish"
VG="/home/scott/bin/vg"

$SEQWISH -k 16 -s $FASTA -p $PAF -g ${OUT_PREFIX}.gfa

$VG view -Fv ${OUT_PREFIX}.gfa | vg mod -X 256 - | vg mod -n - | vg ids -c -| vg mod -X 256 - > ${OUT_PREFIX}.vg
$VG stats -z ${OUT_PREFIX}.vg  ; vg paths -L -v ${OUT_PREFIX}.vg
$VG view -Vg ${OUT_PREFIX}.vg > ${OUT_PREFIX}_normalized.gfa

$VG index -x ${OUT_PREFIX}.xg ${OUT_PREFIX}.vg
$VG paths -X -x ${OUT_PREFIX}.xg > ${OUT_PREFIX}.gam
