#!/bin/bash

unset PYTHONPATH
module purge
. /etc/profile.d/modules.sh
module load gcc/6.5.0
module load cmake/3.19.1
module load RepeatMasker/4.1.1
module load samtools/1.11
module load bedtools/2.29.2
#module load pbmm2/1.7.0

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source /hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/SDA/miniconda3/etc/profile.d/conda.sh
conda activate

module load perl/5.14.4


SDA="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/SDA/SDA"
alias bamtofastq="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/bamtofastq-1.2.0"

