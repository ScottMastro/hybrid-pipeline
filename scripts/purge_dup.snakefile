import os

def add_slash(string): 
    if len(string) == 0: return "./"
    return string + ("" if string.endswith("/") else "/")

PURGE_READS_DIR = add_slash(config["purgereads"])
os.environ['PURGE_READS_DIR'] = str(PURGE_READS_DIR)

FASTA = config["fasta"]

localrules: merge_pafs

###############################################
# PREPARE READS IN FQ.GZ FORMAT IF NECESSARY
# (PACBIO BAM -> FASTQ)
###############################################

rule bam_to_fastq:
    input:  "{prefix}.bam"
    output: temp("{prefix}.fq.gz")
    resources: mem=8
    threads:1
    shell: "samtools bam2fq {input} | gzip > {output}"

rule gzip_fastq:
    input:  expand("{prefix}.{fq}", fq=["fastq", "fq"], allow_missing=True)
    output: temp("{prefix}.fq.gz")
    resources: mem=4
    threads:1
    shell: "gzip -c {input} > {output}"
    
###############################################
# PURGE DUPS
###############################################

rule align_for_purge:
    input:
        reads=PURGE_READS_DIR+"{readid}.fq.gz",
        assembly="{assem}"
    output:
        temp("{assem}.align_{readid}.paf.gz")
    resources: mem=64
    threads:8
    shell: "minimap2 -t 8 -xmap-pb {input.assembly} {input.reads} | gzip -c - > {output}"

def get_read_ids():
    fq_ids = glob_wildcards(PURGE_READS_DIR+"{readid}.fq.gz").readid
    if len(fq_ids) > 0: return fq_ids
    bam_ids = glob_wildcards(PURGE_READS_DIR+"{readid}.bam").readid
    return bam_ids

rule merge_pafs:
    input:
        expand("{{assem}}.align_{readid}.paf.gz", readid=get_read_ids())
    output:
        temp("{assem}.all_reads.paf.gz")
    resources: mem=8
    threads:1
    shell:
        "echo {input} ; cat {input} > {output}"

rule purge_dups_cut:
    input:
        "{assem}.all_reads.paf.gz"
    output:
        cutoffs=temp("{assem}.cutoffs"),
        cov=temp("{assem}.PB.base.cov"),
        wig=temp("{assem}.PB.cov.wig"),
        stat=temp("{assem}.PB.stat"),
        log=temp("{assem}.calcults.log")
    resources: mem=32
    threads:1
    shell:
        """
        dir=`dirname {wildcards.assem}` ; pbcstat -O $dir {input}
        mv $dir/PB.stat {output.stat} ; mv $dir/PB.base.cov {output.cov} 
        mv $dir/PB.cov.wig {output.wig}
        calcuts {output.stat} > {output.cutoffs} 2> {output.log}
        """

rule purge_dups_self_align:
    input:
        "{assem}"
    output:
        split=temp("{assem}.split"),
        align=temp("{assem}.split.self.paf.gz")
    resources: mem=32
    threads:1
    shell:
        "split_fa {input} > {output.split} ; minimap2 -xasm5 -DP {output.split} {output.split} | gzip -c - > {output.align}"
    
rule purge_dups:
    input:
        cutoffs="{assem}.cutoffs",
        cov="{assem}.PB.base.cov",
        align="{assem}.split.self.paf.gz"
    output:
        bed="{assem}.dups.bed",
        log="{assem}.purge_dups.log"
    resources: mem=32
    threads:1
    shell:
        "purge_dups -2 -T {input.cutoffs} -c {input.cov} {input.align} > {output.bed} 2> {output.log}"

rule get_purge_fa:
    input:
        assembly=FASTA,
	bed=FASTA+".dups.bed"
    output:
        purgefa=FASTA+".purged.fa",
        hapfa=FASTA+".hap.fa"
    resources: mem=32
    threads:1
    shell:
        """
        cwd=`pwd` ; dir=`dirname {input.assembly}` ; cd $dir ; bed=`basename {input.bed}` ; assembly=`basename {input.assembly}` ; \
        get_seqs $bed $assembly ; mv $dir/hap.fa {output.hapfa} ; mv $dir/purged.fa {output.purgefa} ; \
        cd $cwd
        """

