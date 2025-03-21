RREAD_PREFIX = "ref_reads"

localrules: get_ref_variants, ref_polish_rename, hap_ref_polish_rename

def add_slash(string): 
    if len(string) == 0: return "./"
    return string + ("" if string.endswith("/") else "/")

OUTDIR = add_slash(config["out"]) if "out" in config else "./polish/"
TIG = config["tig"]
HYBRID_FA = config["hybridfa"]
SDA_DIR = config["sda"]
REF_ALN = config["r"]

'''
#==================================================
# Full alignment to hybrid
#==================================================

def get_read_ids():
    bam_ids = glob_wildcards(REF_READS+"{readid}.bam").readid
    return bam_ids

def strip_gz(f, ext=".gz"):
    if f.endswith(ext):
        f = f[:-len(ext)]
    return f

rule align_ref_to_hybrid:
    input:
        ref=strip_gz(HYBRID_FA),
        reads=REF_READS+"{readid}.bam"
    output:
        bam="{dir}/align_{readid}.refalign.bam",
        bai="{dir}/align_{readid}.refalign.bam.bai"
    resources: mem=64
    threads:16
    shell:
        """
        tempbam={wildcards.dir}/temp.{wildcards.readid}.bam
        pbmm2 align --sort -j 16 {input.ref} {input.reads} $tempbam
        # Reheader BAM to be compliant with PB BAM standards
        samtools view -H $tempbam | sed "s/SO:coordinate/pb:3.0.4\tSO:coordinate/" | \
           samtools reheader - $tempbam > {output.bam}
        rm $tempbam
        samtools index {output.bam}
        """

rule merge_ref_bams:
    input: expand("{{dir}}/align_{readid}.refalign.bam", readid=get_read_ids())
    output: 
        bam="{dir}/hybrid.refalign.bam",
        bai="{dir}/hybrid.refalign.bam.bai"
    resources: mem=8
    shell: "samtools merge {output.bam} {input} ; samtools index {output.bam}"
'''

#==================================================
# Get reads aligned to contig
#==================================================

#extracts reads from hybrid-aligned BAM file and removes reads that were used by SDA
rule fetch_pbmm2_reads:
    input:  
        ref=REF_ALN,
        sda_reads=SDA_DIR+"/sda.phased.readids"
    output: "{dir}/"+RREAD_PREFIX+".bam"
    resources: mem=16
    shell: 
        '''
        samtools view -b {input.ref} '''+TIG+''' > {wildcards.dir}/ref_reads_unfiltered.bam
	samtools view -H {wildcards.dir}/ref_reads_unfiltered.bam > {wildcards.dir}/header.txt
	samtools reheader -c 'grep -v ^@HD' {wildcards.dir}/ref_reads_unfiltered.bam > {wildcards.dir}/ref_reads_reheader.bam
	cat {input.sda_reads} | cut -f2 > {wildcards.dir}/sda.phased.readids.txt
	$PICARD FilterSamReads I={wildcards.dir}/ref_reads_reheader.bam O={wildcards.dir}/ref_reads_filtered.bam READ_LIST_FILE={wildcards.dir}/sda.phased.readids.txt FILTER=excludeReadList
	samtools reheader {wildcards.dir}/header.txt {wildcards.dir}/ref_reads_filtered.bam > {output}
	rm {wildcards.dir}/ref_reads_unfiltered.bam ; rm {wildcards.dir}/ref_reads_reheader.bam ; rm {wildcards.dir}/ref_reads_filtered.bam ; rm {wildcards.dir}/header.txt
	'''
    
#==================================================
# Alignment
#==================================================

rule align_pbmm2:
    input:
        ref="{dir}/{fasta}.fasta",
        reads="{dir}/{bam}.bam"
    output:
        bam="{dir}/{fasta}.{bam}.pbmm2.bam",
        bai="{dir}/{fasta}.{bam}.pbmm2.bam.bai"
    resources: mem=82
    threads:4
    shell:
        "pbmm2 align --sort -j 4 --unmapped {input.ref} {input.reads} {output.bam} ; samtools index {output.bam}"

#==================================================
# Polish
#==================================================

#{align} allows other rules to decide to use pbmm2 or another aligner if needed

rule polish_arrow:
    input:
        fasta="{dir}/{fasta}.fasta",
        bam="{dir}/{fasta}.{reads}.{align}.bam"
    output:
        "{dir}/{fasta}.{reads}.arrow_{align}.fasta"
    resources: mem=82
    threads: 8
    shell:
        "gcpp --log-level INFO -m=0 --algorithm arrow -r {input.fasta} -o {output} {input.bam}"

rule ref_polish_rename:
    input: "{dir}/{prefix}.arrow_pbmm2.fasta"
    output: "{dir}/{prefix}.polish_ref.fasta"
    resources: mem=4
    shell: "mv {input} {output}"
    
rule hap_ref_polish_rename:
    input: "{dir}/{prefix}.arrow_pbmm2.fasta"
    output: "{dir}/{prefix}.hap_polish_ref.fasta"
    resources: mem=4
    shell: "mv {input} {output}"

#==================================================
# Variant Caller
#==================================================

rule run_longshot:
    input:
        ref="{dir}/{fasta}.fasta",
        fai="{dir}/{fasta}.fasta.fai",
        bam="{dir}/{fasta}.{bam}.pbmm2.bam"
    output:
        bam="{dir}/{fasta}.{bam}.longshot.bam",
        bai="{dir}/{fasta}.{bam}.longshot.bam.bai",
        vcf="{dir}/{fasta}.{bam}.longshot.vcf"
    resources: mem=32
    shell:
        """

        if [ `samtools view {input.bam} | head | wc -l` == "0" ]; then
          touch {output.vcf}
          cp {input.bam} {output.bam}
          samtools index {output.bam}

        else

          longshot -F -A --out_bam {output.bam} --bam {input.bam} --ref {input.ref} --out {output.vcf}
        
          if [ ! -f {output.bam} ]; then
            cp {input.bam} {output.bam}
          fi

          samtools index {output.bam}
        fi
        """

rule get_ref_variants:
    input: 
        bam="{dir}/{prefix}."+RREAD_PREFIX+".longshot.bam",
        bai="{dir}/{prefix}."+RREAD_PREFIX+".longshot.bam.bai",
        vcf="{dir}/{prefix}."+RREAD_PREFIX+".longshot.vcf.gz"
    output: 
        bam="{dir}/{prefix}."+RREAD_PREFIX+".variants.bam",
        bai="{dir}/{prefix}."+RREAD_PREFIX+".variants.bam.bai",
        vcf="{dir}/{prefix}."+RREAD_PREFIX+".variants.vcf.gz"
    resources: mem=4
    shell: "mv {input.bam} {output.bam} ; mv {input.bai} {output.bai} ; \
    mv {input.vcf} {output.vcf}"

#==================================================
# Split reads into haplotypes
#==================================================

rule split_reads_hap:
    input:
        "{dir}/"+RREAD_PREFIX+"_haplotag.bam"
    output:
        u="{dir}/"+RREAD_PREFIX+"_unassigned.bam",
        ubai="{dir}/"+RREAD_PREFIX+"_unassigned.bam.bai",
        h1="{dir}/"+RREAD_PREFIX+"_haplotype1.bam",
        h1bai="{dir}/"+RREAD_PREFIX+"_haplotype1.bam.bai",
        h2="{dir}/"+RREAD_PREFIX+"_haplotype2.bam",
        h2bai="{dir}/"+RREAD_PREFIX+"_haplotype2.bam.bai"
    resources: mem=6
    shell:
        """
        prefix={wildcards.dir}/
        samtools view -H {input} > $prefix/header.txt ; \
        samtools view {input} | grep HP:i:1 > $prefix/h1.sam ; \
        samtools view {input} | grep HP:i:2 > $prefix/h2.sam ; \
        samtools view {input} | grep -v HP:i:1 | grep -v HP:i:2 > $prefix/u.sam ; \
        
        cat $prefix/header.txt $prefix/h1.sam | samtools view -b - > {output.h1} ; \
        cat $prefix/header.txt $prefix/h2.sam | samtools view -b - > {output.h2} ; \
        cat $prefix/header.txt $prefix/u.sam  | samtools view -b - > {output.u} ; \
        
        samtools index {output.u} ; samtools index {output.h1} ; samtools index {output.h2} ; \
        rm $prefix/header.txt ; rm $prefix/*.sam         
        """

# to generate random noise for shuffling, but in a consistent way
SEEDER="""
get_seeded_random() \
{{
  seed=$1 ; openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}}

"""

rule split_unphased_reads_random:
    input:
        u="{dir}/"+RREAD_PREFIX+"_unassigned.bam",
        h1="{dir}/"+RREAD_PREFIX+"_haplotype1.bam",
        h2="{dir}/"+RREAD_PREFIX+"_haplotype2.bam"
    output:
        h1=temp("{dir}/"+RREAD_PREFIX+"_haplotype1_iter{m}.bam"),
        h2=temp("{dir}/"+RREAD_PREFIX+"_haplotype2_iter{m}.bam"),
        h1bai=temp("{dir}/"+RREAD_PREFIX+"_haplotype1_iter{m}.bam.bai"),
        h2bai=temp("{dir}/"+RREAD_PREFIX+"_haplotype2_iter{m}.bam.bai"),
    resources: mem=12
    shell:
        SEEDER + """
        prefix={wildcards.dir}/{wildcards.m}_ \

        samtools view {input.u} > ${{prefix}}u.sam ; \
        shuf --random-source=<(get_seeded_random {wildcards.m}) ${{prefix}}u.sam | \
                  split -d -l $(( $(wc -l < ${{prefix}}u.sam) * 50 / 100 + 1 )) - ${{prefix}}usplit ; \

        samtools view -H {input.h1} > ${{prefix}}header.txt ; \
        cat ${{prefix}}header.txt ${{prefix}}usplit*0 | samtools view -b - > ${{prefix}}u1.unsorted.bam
        samtools view -H {input.h2} > ${{prefix}}header.txt ; \
        cat ${{prefix}}header.txt ${{prefix}}usplit*1 | samtools view -b - > ${{prefix}}u2.unsorted.bam
        rm ${{prefix}}header.txt ; rm ${{prefix}}usplit* ; rm ${{prefix}}u.sam ; \

        samtools sort ${{prefix}}u1.unsorted.bam > ${{prefix}}u1.sorted.bam ; \
        samtools sort ${{prefix}}u2.unsorted.bam > ${{prefix}}u2.sorted.bam ; \

        samtools merge {output.h1} {input.h1} ${{prefix}}u1.sorted.bam ; samtools index {output.h1} ; \
        samtools merge {output.h2} {input.h2} ${{prefix}}u2.sorted.bam ; samtools index {output.h2} ; \

        rm ${{prefix}}u*sorted.bam
        """

