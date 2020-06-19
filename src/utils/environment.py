#as of now these values only apply to Long Ranger, the most resource intensive step
CORES=8
MEM=12

#==================================================
# Long Ranger requirements
#==================================================

LONGRANGER="/home/scott/bin/longranger-2.2.2/longranger"
#version of GATK must be 3.3-3.8, or 4 except 3.6
GATK_10x_JAR="/home/scott/bin/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar"
BAM2FQ_10X="/home/scott/bin/bamtofastq-1.2.0"
PICARD="java -Xmx8G -jar /home/scott/bin/picard.jar"

#==================================================
# PacBio requirements
#==================================================

PBMM2 = "pbmm2"
LONGSHOT = "longshot"
GENOMIC_CONSENSUS = "gcpp"
MM2 = "/home/scott/bin/minimap2-2.17_x64-linux/minimap2"
PAFTOOLS = "k8 /home/scott/bin/minimap2-2.17_x64-linux/paftools.js"

#==================================================
# Common tools
#==================================================

SAMTOOLS = "samtools"
BGZIP = "bgzip"
TABIX = "tabix"
BWA = "bwa"
PILON="java -Xmx16G -jar /home/scott/bin/pilon-1.23.jar"

#==================================================
# Phasing and graphs
#==================================================

WHATSHAP = "whatshap"
VG="/home/scott/bin/vg"
SEQWISH="/home/scott/bin/seqwish/bin/seqwish"


