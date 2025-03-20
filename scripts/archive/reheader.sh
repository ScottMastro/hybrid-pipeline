BAM=$1

filename=$(basename -- $BAM)
filename="${filename%.*}"


samtools view -H $BAM | sed "s/SO:coordinate/pb:3.0.4\tSO:coordinate/" | samtools reheader - $BAM > ${filename}.reheader.bam
