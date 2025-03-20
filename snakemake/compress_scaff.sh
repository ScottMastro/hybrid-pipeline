#!/bin/bash

SAMPLE=$1
SCAFF=$2

NAME=${SAMPLE}.${SCAFF}
echo $NAME

if [ ! -f ".done" ]; then
    echo "$NAME - The .done file does not exist. Exiting script."
    exit 1
fi

mkdir -p iterpolish
mkdir -p consensus
mkdir -p reads

zip -j -r jobout.zip jobout
rm -r jobout

sed -i "1s/.*/>${NAME}.unpolished/" unpolished.fasta
gzip unpolished.fasta
mv unpolished.fasta.gz unpolished.fa.gz
rm unpolished.fasta.fai

rm *.t2t.paf
rm *.t2t.chr.txt

zcat query_reads.I1.fastq.gz | awk 'NR%4==1' > query.readids.txt
gzip query.readids.txt

rm query_reads.R*
rm query_reads.I*
rm -r refdata-consensus

samtools view ref_reads.bam | cut -f1 | gzip > reads/ref.readids.txt.gz
gzip query.readids.txt ; mv query.readids.txt.gz reads/
rm ref_reads.bam

mv sda.phased.readids.txt reads/
gzip reads/sda.phased.readids.txt

sed -i "1s/.*/>${NAME}.refpolished/" unpolished/refpolished.fasta
mv unpolished/refpolished.fasta consensus/refpolished.fa
gzip consensus/refpolished.fa
rm -r unpolished

rm consensus/high_confidence_hets.phased.vcf.gz.tbi

sed -i "1s/.*/>${NAME}.consensus/" consensus/consensus.fasta
gzip consensus/consensus.fasta
rm consensus/consensus.fasta.fai
mv consensus/consensus.fasta.gz consensus/consensus.fa.gz

zip -j -r consensus/pilon.zip consensus/pilon/
rm -r consensus/pilon/

rm consensus/*bam
rm consensus/*bai


samtools view reads/ref_hap1.bam | cut -f1 | gzip > reads/ref.hap1.readids.txt.gz
samtools view reads/ref_hap2.bam | cut -f1 | gzip > reads/ref.hap2.readids.txt.gz
samtools view reads/ref_unassigned.bam | cut -f1 | gzip > reads/ref.unassigned.readids.txt.gz

samtools view reads/query_hap1.bam | cut -f1 | sort | uniq | gzip > reads/query.hap1.readids.txt.gz
samtools view reads/query_hap2.bam | cut -f1 | sort | uniq | gzip > reads/query.hap2.readids.txt.gz
samtools view reads/query_unassigned.bam | cut -f1 | sort | uniq | gzip > reads/query.unassigned.readids.txt.gz

rm reads/*bam
rm reads/*bai

zip -j -r iterpolish/iter0.pilon.zip iter0/pilon/
zip -j -r iterpolish/iter1.pilon.zip iter1/pilon/
zip -j -r iterpolish/iter2.pilon.zip iter2/pilon/


sed -i "1s/.*/>${NAME}.iter0_refpolished.hap1/" iter0/refpolished.hap1.fasta
sed -i "1s/.*/>${NAME}.iter0_refpolished.hap2/" iter0/refpolished.hap2.fasta
sed -i "1s/.*/>${NAME}.iter1.hap1/" iter0/consensus.hap1.fasta
sed -i "1s/.*/>${NAME}.iter1.hap2/" iter0/consensus.hap2.fasta

cat iter0/refpolished.hap1.fasta | gzip > iterpolish/iter0.refpolished.hap1.fa.gz
cat iter0/refpolished.hap2.fasta | gzip > iterpolish/iter0.refpolished.hap2.fa.gz
cat iter0/consensus.hap1.fasta | gzip > iterpolish/iter1.hap1.fa.gz
cat iter0/consensus.hap2.fasta | gzip > iterpolish/iter1.hap2.fa.gz

rm -r iter0

sed -i "1s/.*/>${NAME}.iter1_refpolished.hap1/" iter1/refpolished.hap1.fasta
sed -i "1s/.*/>${NAME}.iter1_refpolished.hap2/" iter1/refpolished.hap2.fasta
sed -i "1s/.*/>${NAME}.iter2.hap1/" iter1/consensus.hap1.fasta
sed -i "1s/.*/>${NAME}.iter2.hap2/" iter1/consensus.hap2.fasta

cat iter1/refpolished.hap1.fasta | gzip > iterpolish/iter1.refpolished.hap1.fa.gz
cat iter1/refpolished.hap2.fasta | gzip > iterpolish/iter1.refpolished.hap2.fa.gz
cat iter1/consensus.hap1.fasta | gzip > iterpolish/iter2.hap1.fa.gz
cat iter1/consensus.hap2.fasta | gzip > iterpolish/iter2.hap2.fa.gz

rm -r iter1

sed -i "1s/.*/>${NAME}.iter2_refpolished.hap1/" iter2/refpolished.hap1.fasta
sed -i "1s/.*/>${NAME}.iter2_refpolished.hap2/" iter2/refpolished.hap2.fasta
sed -i "1s/.*/>${NAME}.hap1/" iter2/consensus.hap1.fasta
sed -i "1s/.*/>${NAME}.hap2/" iter2/consensus.hap2.fasta

cat iter2/refpolished.hap1.fasta | gzip > iterpolish/iter2.refpolished.hap1.fa.gz
cat iter2/refpolished.hap2.fasta | gzip > iterpolish/iter2.refpolished.hap2.fa.gz
cat iter2/consensus.hap1.fasta | gzip > polished.hap1.fa.gz
cat iter2/consensus.hap2.fasta | gzip > polished.hap2.fa.gz

rm -r iter2
