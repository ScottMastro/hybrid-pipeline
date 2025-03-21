QREAD_PREFIX = "query_reads"

GATK_JAR=gatk-package-4.0.0.0-local.jar
PICARD_JAR=picard.jar
PILON_JAR=pilon-1.23.jar

FIX_10X_HETS="longranger_vcf_correct.py"
TRIM_10X_BARCODE="trim_10x_barcode.py"
SPLIT_10X_HAP="split_10x_hap.py"

localrules: get_query_variants, query_polish_rename, hap_query_polish_rename, haplotype_iter_bam

def add_slash(string):
    if len(string) == 0: return "./"
    return string + ("" if string.endswith("/") else "/")

OUTDIR = add_slash(config["out"]) if "out" in config else "./polish/"
TIG = config["tig"]
HYBRID_FA = config["hybridfa"]
QUERY_ALN = config["q"]

'''
#==================================================
# Full alignment to hybrid
#==================================================

def strip_gz(f, ext=".gz"):
    if f.endswith(ext):
        f = f[:-len(ext)]
    return f

rule align_query_to_hybrid:
    input:
        readsdir=QUERY_READS,
        ref=strip_gz(HYBRID_FA)
    output:
        bam="{dir}/hybrid.queryalign.bam",
        bai="{dir}/hybrid.queryalign.bam.bai"
    threads:32
    resources: mem=121
    shell:
        """
        readsdir=`readlink -f {input.readsdir}` ; hybridfa=`readlink -f {input.ref}`
        ln -s $readsdir {wildcards.dir}/reads_query
        ln -s $hybridfa {wildcards.dir}
        cwd=`pwd` ; cd {wildcards.dir}
        filename=`basename {input.ref}` ; prefix=${{filename%.*}} ; idx=refdata-$prefix
        longranger mkref $filename
        $PICARD CreateSequenceDictionary R=$idx/fasta/genome.fa O=$idx/fasta/genome.dict
        longranger align --id=lr --fastqs=reads_query/ --reference=$idx \
           --jobmode=local --localcores=32 --localmem={resources.mem}
        
        mv {wildcards.dir}/lr/outs/phased_possorted_bam.bam {output.bam} ; samtools index {output.bam}
        rm -r {wildcards.dir}/lr ; unlink reads_query ; unlink $filename
        cd $cwd
        """
'''

#==================================================
# Get reads aligned to contig
#==================================================

rule fetch_longranger_reads:    
    input: QUERY_ALN
    output:
        r1="{dir}/"+QREAD_PREFIX+".R1.barcode.fastq.gz",
        r2="{dir}/"+QREAD_PREFIX+".R2.fastq.gz",
        i1="{dir}/"+QREAD_PREFIX+".I1.fastq.gz"
    resources: mem=16
    shell:
        """
        prefix={wildcards.dir}/q 
        bamtofastq --locus="""+TIG+""":0-999999999 {input} $prefix 
        mv $prefix/*/*R1_001.fastq.gz {output.r1}
        mv $prefix/*/*R2_001.fastq.gz {output.r2}
        mv $prefix/*/*I1_001.fastq.gz {output.i1}
        rm -r $prefix
	"""

#==================================================
# Alignment
#==================================================
    
rule align_bwa:
    input:
        ref="{dir}/{fasta}.fasta",
        r1="{dir}/{fq}.R1.fastq.gz",
        r2="{dir}/{fq}.R2.fastq.gz",
    output:
        bam="{dir}/{fasta}.{fq}.bwa.bam",
        bai="{dir}/{fasta}.{fq}.bwa.bam.bai"
    resources: mem=32
    shell:
        """
        bwa index {input.ref} ; \
        bwa mem {input.ref} {input.r1} {input.r2} | samtools sort -o {output.bam} - ; \
        samtools index {output.bam} ; rm {input.ref}.*
        """
        
rule align_longranger:
    input:
        mkref="{dir}/refdata-{fasta}",
        r1="{dir}/{fq}.R1.barcode.fastq.gz",
        r2="{dir}/{fq}.R2.fastq.gz",
        i1="{dir}/{fq}.I1.fastq.gz"
    output:
        bam="{dir}/{fasta}.{fq}.longranger.bam",
        bai="{dir}/{fasta}.{fq}.longranger.bam.bai",
	vcf="{dir}/{fasta}.{fq}.longranger.vcf.gz"
    threads:32
    resources: mem=121
    shell:
        """
        mkdir -p {wildcards.dir}/reads_query
        cp {input.r1} {wildcards.dir}/reads_query/reads_S1_L000_R1_001.fastq.gz 
        cp {input.r2} {wildcards.dir}/reads_query/reads_S1_L000_R2_001.fastq.gz
        cp {input.i1} {wildcards.dir}/reads_query/reads_S1_L000_I1_001.fastq.gz 
        cwd=`pwd` ; cd {wildcards.dir}
        longranger wgs --id=lr --fastqs=reads_query/ --sample=reads --reference=refdata-{wildcards.fasta} \
           --jobmode=local --disable-ui --sex=female --nopreflight --noloupe --vcmode=gatk:$GATK_JAR \
           --localcores=16 --localmem={resources.mem}
        cd $cwd
        rm -r {wildcards.dir}/reads_query
        python $FIX_10X_HETS {wildcards.dir}/lr/outs/phased_variants.vcf.gz
        bgzip {wildcards.dir}/lr/outs/phased_variants.corrected.vcf
        mv {wildcards.dir}/lr/outs/phased_variants.corrected.vcf.gz {output.vcf}
        mv {wildcards.dir}/lr/outs/phased_possorted_bam.bam {output.bam} ; samtools index {output.bam}
        rm -r {wildcards.dir}/lr
        """

#==================================================
# Polish
#==================================================

# {align} allows other rules to decide to use longranger or bwa

rule polish_pilon:
    input:
        fasta="{dir}/{fasta}.fasta",
        bam="{dir}/{fasta}.{reads}.{align}.bam"
    output:
        "{dir}/{fasta}.{reads}.pilon_{align}.fasta"
    resources:
        mem=32
    shell:
        "$PILON_JAR --genome {input.fasta} --frags {input.bam} --changes --tracks \
                --outdir {wildcards.dir} --output {wildcards.fasta}_{wildcards.reads}_pilon ; \
        mv {wildcards.dir}/{wildcards.fasta}_{wildcards.reads}_pilon.fasta {output} ; sed -i '1s/.*/>consensus/' {output} "

rule polish_pilon_final:
    input:
        fasta="{dir}/{fasta}.fasta",
        bam="{dir}/{fasta}.{reads}.{align}.bam"
    output:
        "{dir}/{fasta}.{reads}.pilon_final_{align}.fasta"
    resources:
        mem=32
    shell:
        "$PILON_JAR --genome {input.fasta} --frags {input.bam} --changes --tracks --fix snps,indels \
                --outdir {wildcards.dir} --output {wildcards.fasta}_{wildcards.reads}_pilon; \
        mv {wildcards.dir}/{wildcards.fasta}_{wildcards.reads}_pilon.fasta {output} ; sed -i '1s/.*/>consensus/' {output} "


rule query_polish_rename:
    input: "{dir}/{prefix}.pilon_bwa.fasta"
    output: "{dir}/{prefix}.polish_query.fasta"
    resources: mem=4
    shell: "mv {input} {output}"
    
rule hap_query_polish_rename:
    input: "{dir}/{prefix}.pilon_bwa.fasta"
    output: "{dir}/{prefix}.hap_polish_query.fasta"
    resources: mem=4
    shell: "mv {input} {output}"

rule hap_query_polish_final_rename:
    input: "{dir}/{prefix}.pilon_final_bwa.fasta"
    output: "{dir}/{prefix}.hap_polish_query_final.fasta"
    resources: mem=4
    shell: "mv {input} {output}"

#==================================================
# Variant Caller
#==================================================

# longranger calls variants after alignment

rule get_query_variants:
    input: "{dir}/{prefix}."+QREAD_PREFIX+".longranger.{ext}",
    output: "{dir}/{prefix}."+QREAD_PREFIX+".variants.{ext}"
    resources: mem=4
    shell: "cp {input} {output}"

#==================================================
# Split reads into haplotypes
#==================================================

rule split_10x_hap:    
    input:
        "{dir}/"+QREAD_PREFIX+"_haplotag.bam"
    output:
        u="{dir}/"+QREAD_PREFIX+"_unassigned.bam",
        ubai="{dir}/"+QREAD_PREFIX+"_unassigned.bam.bai",
        h1="{dir}/"+QREAD_PREFIX+"_haplotype1.bam",
        h1bai="{dir}/"+QREAD_PREFIX+"_haplotype1.bam.bai",
        h2="{dir}/"+QREAD_PREFIX+"_haplotype2.bam",
        h2bai="{dir}/"+QREAD_PREFIX+"_haplotype2.bam.bai"
    resources: mem=64
    shell:
        "python $SPLIT_10X_HAP {input} {wildcards.dir} "+QREAD_PREFIX+"_ ; \
        samtools index {output.h1} ; samtools index {output.h2} ; samtools index {output.u}"

rule haplotype_iter_bam:
    input:
        ubam="{dir}/"+QREAD_PREFIX+"_unassigned.bam",
        ubai="{dir}/"+QREAD_PREFIX+"_unassigned.bam.bai",
        bam="{dir}/"+QREAD_PREFIX+"_haplotype{h}.bam",
        bai="{dir}/"+QREAD_PREFIX+"_haplotype{h}.bam.bai"
    output:
        bam="{dir}/"+QREAD_PREFIX+"_haplotype{h}_iter{x}.bam",
        bai="{dir}/"+QREAD_PREFIX+"_haplotype{h}_iter{x}.bam.bai"
    resources: mem=6
    shell:
        "samtools merge {output.bam} {input.bam} {input.ubam} ; samtools index {output.bam}"

#==================================================
# 10x helper functions
#==================================================  

rule trim_10x_barcodes:    
    input: "{prefix}.R1.barcode.fastq.gz",
    output: "{prefix}.R1.fastq.gz",
    resources: mem=6
    shell: "python $TRIM_10X_BARCODE {input} {output}"

rule longranger_mkref:
    input:
        "{dir}/{fasta}.fasta"
    output:
        ref=directory("{dir}/refdata-{fasta}"),
        fa="{dir}/refdata-{fasta}/fasta/genome.fa",
        dict="{dir}/refdata-{fasta}/fasta/genome.dict"
    resources:
        mem=16
    shell:
      	"""
        PWD=`pwd` ; rm -r {output.ref} ; cd {wildcards.dir} 
        longranger mkref {wildcards.fasta}.fasta ; cd $PWD
      	$PICARD_JAR CreateSequenceDictionary R={output.fa} O={output.dict} ; \
      	touch {output.ref} ; touch {output.ref}/* ; touch {output.ref}/*/*
        """

rule longranger_bamtofastq:
    input:
        "{dir}/{prefix}.bam"
    output:
        r1="{dir}/{prefix}.R1.barcode.fastq.gz",
        r2="{dir}/{prefix}.R2.fastq.gz",
        i="{dir}/{prefix}.I1.fastq.gz"
    resources: mem=6
    shell:
        """
        Q={wildcards.dir}/qreads_{wildcards.prefix}_
        rm -rf $Q
        bamtofastq {input} $Q

	sleep 5
	#if no reads exist
	touch {output.r1}
	touch {output.r2}
	touch {output.i}

        mv $Q/*/*_R1_*.fastq.gz {output.r1}
        mv $Q/*/*_R2_*.fastq.gz {output.r2}
        mv $Q/*/*_I1_*.fastq.gz {output.i}
        rm -rf $Q

        """

