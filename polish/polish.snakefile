import os

HIGH_CONF_HETS="high_conf_hets.py"
CHR_COUNT="chr_count.py"

ENV_FILE=config["env"]
shell.prefix("source "+ENV_FILE+" ; ")

QREAD_PREFIX = "query_reads"
RREAD_PREFIX = "ref_reads"

include: "polish_ref.snakefile"
include: "polish_query.snakefile"

def add_slash(string): 
    if len(string) == 0: return "./"
    return string + ("" if string.endswith("/") else "/")

OUTDIR = add_slash(config["out"]) if "out" in config else "./polish/"
if not os.path.exists(OUTDIR): os.mkdir(OUTDIR)
OUTDIR = OUTDIR + TIG + "/"
if not os.path.exists(OUTDIR): os.mkdir(OUTDIR)

REF_ALN = config["r"]
QUERY_ALN = config["q"]
HYBRID_FA = config["hybridfa"]

SAMPLE = config["sample"]
ITERS= str(config["iters"])
TIG = config["tig"]

# T2T assembly
T2T = config["t2t"]

localrules: fetch_hybrid_tig, rename_consensus, rename_hap_consensus, first_iter, next_iter, bgzip, done, fai

rule all: 
    input: OUTDIR+".done"

rule bgzip:
    input: "{prefix}.vcf"
    output:
        gz="{prefix}.vcf.gz",
        tbi="{prefix}.vcf.gz.tbi"
    resources: mem=4
    shell: "bgzip {input} ; tabix {output.gz}"


rule fai:
    input: "{prefix}.fasta"
    output: "{prefix}.fasta.fai"
    resources: mem=4
    shell: "samtools faidx {input}"
    
#==================================================
# Fetch contig + generate consensus
#==================================================

rule fetch_hybrid_tig:
    input: HYBRID_FA
    output: 
        fa="{dir}/unpolished.fasta",
        fai="{dir}/unpolished.fasta.fai"
    resources: mem=6
    shell: "samtools faidx {input} "+TIG+" > {output.fa} ; samtools faidx {output.fa}"

rule rename_consensus:
    input: "{dir}/unpolished."+RREAD_PREFIX+".polish_ref."+QREAD_PREFIX+".polish_query.fasta"
    output: "{dir}/consensus.fasta"
    resources: mem=6
    shell: "mv {input} {output}"


rule aln_to_t2t:
    input: "{dir}/unpolished.fasta"
    output: 
        paf="{dir}/unpolished.t2t.paf",
        chr="{dir}/unpolished.t2t.chr.txt"
    resources: mem=16
    shell:
        """ 
        minimap2 -cx asm5 $T2T {input} > {output.paf}
        python $CHR_COUNT {output.paf} {output.chr}
        """

#==================================================
# Phase haplotypes
#==================================================

rule generate_high_conf_hets:
    input:
        rvcf="{dir}/consensus."+RREAD_PREFIX+".variants.vcf.gz",
        qvcf="{dir}/consensus."+QREAD_PREFIX+".variants.vcf.gz",
        chr="{dir}/unpolished.t2t.chr.txt"
    output:
        "{dir}/high_confidence_hets.vcf"
    resources: mem=6
    shell:
        "python $HIGH_CONF_HETS {input.rvcf} {input.qvcf} {output}"

rule whatshap_phase:
    input:
        ref="{dir}/consensus.fasta",
        vcf="{dir}/high_confidence_hets.vcf.gz",
        rbam="{dir}/consensus."+RREAD_PREFIX+".variants.bam",
        rbai="{dir}/consensus."+RREAD_PREFIX+".variants.bam.bai",
        qbam="{dir}/consensus."+QREAD_PREFIX+".variants.bam",
        qbai="{dir}/consensus."+QREAD_PREFIX+".variants.bam.bai"
    output:
        "{dir}/high_confidence_hets.whatshap.vcf"
    resources: mem=16
    shell:
        """
        if [ "$( zcat {input.vcf} | wc -l )" -lt 3 ]
        then
          touch {output}
        else        
          whatshap phase -o {output} --ignore-read-groups --max-coverage 20 --indels --distrust-genotypes \
             --reference {input.ref} {input.vcf} {input.rbam} {input.qbam}
        fi
        """
                
rule whatshap_haplotag:
    input:
        ref="{dir}/consensus.fasta",
        vcf="{dir}/high_confidence_hets.whatshap.vcf.gz",
        rbam="{dir}/consensus."+RREAD_PREFIX+".variants.bam",
        qbam="{dir}/consensus."+QREAD_PREFIX+".variants.bam"
    output:
        rtag="{dir}/"+RREAD_PREFIX+"_haplotag.bam",
        rtagbai="{dir}/"+RREAD_PREFIX+"_haplotag.bam.bai",
        qtag="{dir}/"+QREAD_PREFIX+"_haplotag.bam",
        qtagbai="{dir}/"+QREAD_PREFIX+"_haplotag.bam.bai"
    resources: mem=16
    shell:
        """
        whatshap haplotag -o {output.rtag} --ignore-read-groups --reference {input.ref} {input.vcf} {input.rbam} ; \
        whatshap haplotag -o {output.qtag} --ignore-read-groups --reference {input.ref} {input.vcf} {input.qbam} ; \
        samtools index {output.qtag} ; samtools index {output.rtag}
        """

#==================================================
# Iterative haplotype polishing
# unpolished -> consensus
#==================================================

rule first_iter:
    input: 
        fa="{dir}/consensus.fasta",
        fai="{dir}/consensus.fasta.fai"
    output: 
        fa="{dir}/iter0.hap{x}.unpolished.fasta",
        fai="{dir}/iter0.hap{x}.unpolished.fasta.fai"
    resources: mem=4
    shell: "cp {input.fa} {output.fa} ; cp {input.fai} {output.fai}"

def m_minus_one(wildcards):
	n=str(int(wildcards.m)-1)
	return(wildcards.dir+"/iter"+n+".hap"+wildcards.x+".consensus.fasta")
def m_minus_one_fai(wildcards):
	return (m_minus_one(wildcards) + ".fai")

rule next_iter:
    input: 
        fa=m_minus_one,
        fai=m_minus_one_fai
    output: 
        fa="{dir}/iter{m,[1-9]{1}\d?}.hap{x}.unpolished.fasta",
        fai="{dir}/iter{m,[1-9]{1}\d?}.hap{x}.unpolished.fasta.fai"
    resources: mem=4
    shell: "cp {input.fa} {output.fa} ; cp {input.fai} {output.fai}"

rule rename_hap_consensus:
    input: 
        fasta="{dir}/iter{m}.hap{x}.unpolished."+RREAD_PREFIX+"_haplotype{x}_iter{m}.hap_polish_ref."+QREAD_PREFIX+"_haplotype{x}_iter{m}.hap_polish_query_final.fasta"
    output: 
        fasta="{dir}/iter{m}.hap{x}.consensus.fasta"
    resources: mem=4
    shell: "mv {input.fasta} {output.fasta}"

rule done:
    input:
        a="{outdir}/iter"+ITERS+".hap1.consensus.fasta",
        b="{outdir}/iter"+ITERS+".hap2.consensus.fasta"
    output:
        "{outdir}/.done"
    shell:
        """
        samtools faidx {input.a}
        samtools faidx {input.b}
        touch {output}
        """

