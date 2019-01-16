#!/usr/bin/env bash

#PBS -l mem=64gb
#PBS -l vmem=64gb
#PBS -l gres=localhd:10
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -o /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/pipeline/01_align_contigs/out/jobout
#PBS -e /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/pipeline/01_align_contigs/out/jobout
#PBS -N blockmaker



#dir=$1
#: ${dir:=./out}

dir=./out
cd /hpf/largeprojects/cfcentre/strug/10XG_fromTCAG/pipeline/01_align_contigs

echo "combining alignment files"
#cat $dir/*.out > $/dir/blastn.vs_canu.txt
#rm $dir/*.out

# all options required
BLAST_TABULAR_OUTPUT=$dir/alignments/blastn.vs_canu.txt

ID=95				# filter: minimum identity required 
COV=40				# filter: minimum percentage of query coverage required
PREFIX=$dir/nova_vs_canu	# prefix for output file (eg: output_dir/blast.query_vs_reference)

echo "filter blast results"
#cat $BLAST_TABULAR_OUTPUT | perl -naF'\t' -e 'BEGIN{ $id='$ID'; $cov='$COV' } $p=(($F[7]-$F[6]+1)/$F[12])*100; if($F[2]>=$id && $p>=$cov){ chomp($F[1]); @s=split(":",$F[0]); $n=$s[3]; $n=~s/(exon|part)//g; printf("$s[1]\t$s[2]\t$n\t$F[1]\t$F[8]\t$F[9]\t$F[12]\t$F[13]\t%.2f\t$_",$p); }' > $PREFIX.id$ID.cov$COV.filter.txt

echo "sort filtered hits"
#cat $PREFIX.id$ID.cov$COV.filter.txt | sort -k1,1d -k3,3n -k4,4n -k5,5n -k6,6nr > $PREFIX.id$ID.cov$COV.filter.sort.txt

echo "summarize information"
#cat $PREFIX.id$ID.cov$COV.filter.sort.txt | perl -nae '$k="$F[0]\t$F[1]\t$F[3]\t$F[7]"; if(! exists($h{$k})){ $h{$k}=$F[2]; $o[$i]=$k; $hs{$k}=$F[4]; $he{$k}=$F[5]; $hqs{$k}=$F[15]; $hqe{$k}=$F[16]; $hql{$k}=$F[21]; $i++; } else{ $h{$k}.=",".$F[2]; $hs{$k}.=",".$F[4]; $he{$k}.=",".$F[5]; $hqs{$k}.=",".$F[15]; $hqe{$k}.=",".$F[16]; $hql{$k}.=",".$F[21]; } END{ foreach(@o){ print "$_\t$h{$_}\t$hs{$_}\t$he{$_}\t$hqs{$_}\t$hqe{$_}\t$hql{$_}\n" } }' | perl -naF'\t' -e 'chomp($F[9]); @s=split(",",$F[4]); chomp; print join("\t",@F[0..8])."\t$s[0]\t".@s."\t$F[9]\n";' | sort -k1,1d -k10,10n -k11,11nr -k4,4nr > $PREFIX.id$ID.cov$COV.filter.sort.summarize.txt


module load python/2.7.9
echo "construct blocks"
python python/main.py $PREFIX.id$ID.cov$COV.filter.sort.summarize.txt $PREFIX.id$ID.cov$COV.blocks.txt
echo "contigs"
python python/contigs.py $PREFIX.id$ID.cov$COV.blocks.txt $PREFIX.id$ID.cov$COV.blocks_contigs.txt
echo "summarize blocks"
python python/BlockStats.py $PREFIX.id$ID.cov$COV.blocks.txt > $PREFIX.id$ID.cov$COV.genomecov.txt

mkdir -p ../blocks
cp $PREFIX.id$ID.cov$COV.filter.sort.summarize.txt ../blocks/summary.txt
cp $PREFIX.id$ID.cov$COV.blocks.txt ../blocks/blocks.txt

