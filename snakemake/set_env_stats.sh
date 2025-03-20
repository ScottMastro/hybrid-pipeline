#!/usr/bin/env bash

DIR_=/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/anaconda
__conda_setup="$(${DIR_}'/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "${DIR_}/etc/profile.d/conda.sh" ]; then
    . "${DIR_}/etc/profile.d/conda.sh"
    else
        export PATH="${DIR_}/bin:$PATH"
    fi
fi
unset __conda_setup
echo "Active python version: " `which python`


module load busco/5.2.2
# for offline dataset - https://busco.ezlab.org/busco_userguide.html
BUSCO_DATASETS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/busco_downloads"

module load quast
HG38="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/hg38.fa.gz"
HG38_ANNOTATIONS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/reference/hg38/gencode.v39.annotation.gff3.gz"



module load samtools/1.11
module load bcftools/1.6
module load java/1.8.0_91
module load bwa/0.7.8
module load bedtools/2.29.2
module load canu/1.8
module load minimap2/2.11
module load blast+/2.7.1

#------------------------------------------------

module load perl/5.26.2
ASSEMBLY_STATS=/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/assembly-stats

#gcpp in anaconda
