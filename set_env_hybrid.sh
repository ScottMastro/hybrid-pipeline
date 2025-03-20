#!/usr/bin/env bash

module load samtools/1.11
module load bcftools/1.6
module load java/1.8.0_91
module load bwa/0.7.8
module load bedtools/2.29.2
module load canu/1.8
module load minimap2/2.11
module load blast+/2.7.1

T2T=/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/t2t_reference/hs1.fa.mmi

#add purge_dups binary files
PATH="$PATH:/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/purge_dups/bin"

HIGH_CONF_HETS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/high_conf_hets.py"
CLEANUP="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/cleanup.sh"
COMPRESS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/compress_scaff.sh"

FIX_10X_HETS="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/longranger_vcf_correct.py"
TRIM_10X_BARCODE="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/trim_10x_barcode.py"
SPLIT_10X_HAP="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/split_10x_hap.py"
CHR_COUNT="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/1_snakemake/chr_count.py"

#------------------------------------------------


#longranger, bamtofastq
PATH="$PATH:/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools"
PICARD="java -jar /hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/picard.jar"
PILON="java -jar /hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/pilon-1.23.jar"

#version of GATK must be 3.3-3.8, or 4 except 3.6
GATK_JAR="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/gatk-package-4.0.0.0-local.jar"

alias bamtofastq="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/bamtofastq-1.2.0"

#gcpp in anaconda
