#!/bin/bash
#create a file of file names for input: ls -d `pwd`/*/transposed_report.tsv > reports.txt
#create a file of file names for input: ls -d CF* > ids.txt

REPORTS=$1
IDS=$2
OUT="quast_stats.txt"

HEADER="id,assembly,filename,total_contigs,contigs_1kb,contigs_5kb,contigs_10kb,contigs_25kb,contigs_50kb,total_span,span_1kb,span_5kb,span_10kb,span_25kb,span_50kb,contigs,largest,span,ref_span,gc,ref_gc,n50,ng50,n75,ng75,l50,lg50,l75,lg75,misassembly,misassembly_contig,misassembly_span,local_misassembly,scaffold_gap,possible_ mge,unaligned_misassembly,unaligned_contigs,unaligned_span,genome_fraction,dup_ratio,n_100kb,mismatches_100kb,indel_100kb,genomic_features,longest_alignment,total_alignment,na50,nga50,na75,nga75,la50,lga50,la75,lga75"
echo $HEADER > $OUT



IFS=$'\n' read -d '' -r -a ID < $IDS

COLS="___columns___.txt"

i=0
while IFS= read -r REPORT; do

   echo "${ID[$i]},canu" > $COLS
   echo "${ID[$i]},supernova" >> $COLS
   echo "${ID[$i]},hybrid" >> $COLS
   echo "${ID[$i]},leftovers" >> $COLS
   i=$((i+1))

   TEMP="___temp___.txt"
   echo "$REPORT"
   tail -n +2 $REPORT | sed 's/\t/,/g' >> $TEMP
   paste -d "," $COLS $TEMP >> $OUT
   rm $TEMP


done < "$REPORTS"

rm $COLS
