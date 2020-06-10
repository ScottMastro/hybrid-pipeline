#!/bin/bash

TARGET_FA=$1
QUERY_FA=$2
OUT_PAF=$3

/home/scott/bin/minimap2-2.17_x64-linux/minimap2 -cx asm5 --secondary=no $TARGET_FA $QUERY_FA > $OUT_PAF