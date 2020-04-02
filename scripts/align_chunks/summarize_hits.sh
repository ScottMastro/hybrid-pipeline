#!/usr/bin/env bash

DIR=$1
PREFIX_NAME="${2:-"query_vs_ref"}"    # prefix for output file 

ID="${3:-"95"}"          # filter: minimum identity required
COV="${4:-"40"}"         # filter: minimum percentage of query coverage required

PREFIX=$DIR/$PREFIX_NAME

echo "combining alignment files"

BLAST_TABULAR_OUTPUT=$DIR/blastn.vs_ref.txt

cat $DIR/alignments/*.out > $BLAST_TABULAR_OUTPUT
#rm $DIR/alignments/*.out


#echo "filter blast results"
cat $BLAST_TABULAR_OUTPUT | perl -naF'\t' -e 'BEGIN{ $id='$ID'; $cov='$COV' } $p=(($F[7]-$F[6]+1)/$F[12])*100; if($F[2]>=$id && $p>=$cov){ chomp($F[1]); @s=split(":",$F[0]); $n=$s[3]; $n=~s/(exon|part)//g; printf("$s[1]\t$s[2]\t$n\t$F[1]\t$F[8]\t$F[9]\t$F[12]\t$F[13]\t%.2f\t$_",$p); }' > $PREFIX.id$ID.cov$COV.filter.txt

#echo "sort filtered hits"
cat $PREFIX.id$ID.cov$COV.filter.txt | sort -k1,1d -k3,3n -k4,4n -k5,5n -k6,6nr > $PREFIX.id$ID.cov$COV.filter.sort.txt

echo "summarize information"
#cat $PREFIX.id$ID.cov$COV.filter.sort.txt | perl -nae '$k="$F[0]\t$F[1]\t$F[3]\t$F[7]"; if(! exists($h{$k})){ $h{$k}=$F[2]; $o[$i]=$k; $hs{$k}=$F[4]; $he{$k}=$F[5]; $hqs{$k}=$F[15]; $hqe{$k}=$F[16]; $hql{$k}=$F[21]; $i++; } else{ $h{$k}.=",".$F[2]; $hs{$k}.=",".$F[4]; $he{$k}.=",".$F[5]; $hqs{$k}.=",".$F[15]; $hqe{$k}.=",".$F[16]; $hql{$k}.=",".$F[21]; } END{ foreach(@o){ print "$_\t$h{$_}\t$hs{$_}\t$he{$_}\t$hqs{$_}\t$hqe{$_}\t$hql{$_}\n" } }' | perl -naF'\t' -e 'chomp($F[9]); @s=split(",",$F[4]); chomp; print join("\t",@F[0..8])."\t$s[0]\t".@s."\t$F[9]\n";' | sort -k1,1d -k10,10n -k11,11nr -k4,4nr > $PREFIX.id$ID.cov$COV.filter.sort.summarize.txt

#added percent identity column
cat $PREFIX.id$ID.cov$COV.filter.sort.txt | perl -nae '$k="$F[0]\t$F[1]\t$F[3]\t$F[7]"; if(! exists($h{$k})){ $h{$k}=$F[2]; $o[$i]=$k; $hs{$k}=$F[4]; $he{$k}=$F[5]; $hqs{$k}=$F[15]; $hqe{$k}=$F[16]; $hql{$k}=$F[21]; $hid{$k}=$F[11]; $i++; } else{ $h{$k}.=",".$F[2]; $hs{$k}.=",".$F[4]; $he{$k}.=",".$F[5]; $hqs{$k}.=",".$F[15]; $hqe{$k}.=",".$F[16]; $hql{$k}.=",".$F[21]; $hid{$k}.=",".$F[11]; } END{ foreach(@o){ print "$_\t$h{$_}\t$hs{$_}\t$he{$_}\t$hqs{$_}\t$hqe{$_}\t$hql{$_}\t$hid{$_}\n" } }' | perl -naF'\t' -e 'chomp($F[10]); @s=split(",",$F[4]); chomp; print join("\t",@F[0..8])."\t$s[0]\t".@s."\t$F[9]\t$F[10]\n";' | sort -k1,1d -k10,10n -k11,11nr -k4,4nr > $PREFIX.id$ID.cov$COV.summary.txt

#remove intermediate directories
rm $PREFIX.id$ID.cov$COV.filter.txt
rm $PREFIX.id$ID.cov$COV.filter.sort.txt
rm $BLAST_TABULAR_OUTPUT
