import os

REF_FASTA = config["r"]
QUERY_FASTA = config["q"]
OUTDIR = config["out"] +"/"

localrules: all, split_query, combine_blast

###############################################
# CONSTRUCT BLOCKS
###############################################
scattergather:
    split=50

rule all:
    input: 
        blocks=OUTDIR+"blastn.summary.txt"

rule make_blast_db:
    resources: mem=32
    input:  REF_FASTA
    output: 
        nsq=temp("{dir}ref_db.nsq"),
        nhr=temp("{dir}ref_db.nhr"),
        nin=temp("{dir}ref_db.nin")
    shell:
        """
        # handle gzip
        DECOMP_CMD="cat {input}"
        [[ "{input}" == *.gz ]] && DECOMP_CMD="zcat {input}"
        $DECOMP_CMD | makeblastdb -in - -title ref_db -dbtype nucl -out {wildcards.dir}ref_db
        touch {wildcards.dir}ref_db
        """

rule split_query:
    input: QUERY_FASTA
    output:
        chunks=temp("{dir}query.1000bp.fa"),
        chunkcount=temp("{dir}count.txt")
    resources: mem=4
    shell:
        """
        # Define a command to handle decompression
        DECOMP_CMD="cat {input}"
        [[ "{input}" == *.gz ]] && DECOMP_CMD="zcat {input}"
        $DECOMP_CMD | perl -nae 'if($_=~/^>/){{ if($.!=1){{ print "\\n"; }} print $_; }} else{{ chomp; print $_; }}' | \
        perl -nae 'BEGIN{{ $split='1000'; }} chomp; if($_=~/^>/){{ @id=@F; $f=$_; $i=1; $n=0; }} else{{ if($n==0){{ $n=length($_)/1000; @num=split(/\./, $n); if($num[1] ne "00" || $num ne ""){{ $n=$num[0]+1; }} }} for($l=0; $l<=length($_); $l+=$split){{ print ">predicted:".substr($id[0],1).":$n:part$i\\n".substr($_, $l, $split)."\\n"; $i++; }} }}' > {output.chunks} ; \
        n=`cat {output.chunks} | wc -l` ; n=$((n/2)) ; echo $n > {output.chunkcount}
        """

rule split_chunks:
    input:
        chunks="{dir}query.1000bp.fa",
        chunkcount="{dir}count.txt"
    output:
        temp("{dir}query.{scatteritem}.fa")
    resources: mem=4
    shell:
        """
        n=`cat {input.chunkcount}` ; step=$((n/50+1)) ; m=`basename {output} | grep -o -E '[0-9]+' | head -1` ; m=$((m-1)) ; \
        from=$((m * step * 2 + 1)) ; to=$((step * 2 + from)) ; sed -n "$((from)),$((to-1))p;$((to))q" {input.chunks} > {output}
        """

rule blast_chunks:
    input:
        part="{dir}query.{scatteritem}.fa",
        nsq="{dir}ref_db.nsq",
        nhr="{dir}ref_db.nhr",
        nin="{dir}ref_db.nin"
    output:
        temp("{dir}blastn.{scatteritem}.out")
    resources: mem=64
    threads:8
    shell:
       """
       db={input.nsq} ; dbbase=${{db%.*}} ; \
       blastn -query {input.part} -db $dbbase -out {output} -outfmt 6' std qlen slen gaps qseq sseq sframe' -evalue 1e-35 -num_threads 8 -max_target_seqs 10 -max_hsps 5
       """

rule combine_blast:
    input:  gather.split("{{dir}}blastn.{scatteritem}.out")
    output: temp("{dir}blastn.out")
    resources: mem=4
    shell: "cat {input} > {output}"

rule create_tabular_output_filter:
    input:  "{dir}blastn.out"
    output: temp("{dir}blastn.filter.out")
    resources: mem=6
    shell:
        """
        cat {input} | perl -naF'\\t' -e 'BEGIN{{ $id='95'; $cov='40' }} $p=(($F[7]-$F[6]+1)/$F[12])*100; if($F[2]>=$id && $p>=$cov){{ chomp($F[1]); @s=split(":",$F[0]); $n=$s[3]; $n=~s/(exon|part)//g; printf("$s[1]\\t$s[2]\\t$n\\t$F[1]\\t$F[8]\\t$F[9]\\t$F[12]\\t$F[13]\\t%.2f\\t$_",$p); }}' > {output}
	"""

rule create_tabular_output_sort:
    input:  "{dir}blastn.filter.out"
    output: temp("{dir}blastn.filter.sort.out")
    resources: mem=6
    shell: "cat {input} | sort -k1,1d -k3,3n -k4,4n -k5,5n -k6,6nr > {output}"

rule create_tabular_output_summary:
    input:
        "{dir}block/blastn.filter.sort.out"
    output:
        "{dir}blastn.summary.txt"
    resources: mem=6
    shell:
        """
        mkdir -p {wildcards.dir}block ; \
        cat {input} | perl -nae '$k="$F[0]\\t$F[1]\\t$F[3]\\t$F[7]"; if(! exists($h{{$k}})){{ $h{{$k}}=$F[2]; $o[$i]=$k; $hs{{$k}}=$F[4]; $he{{$k}}=$F[5]; $hqs{{$k}}=$F[15]; $hqe{{$k}}=$F[16]; $hql{{$k}}=$F[21]; $hid{{$k}}=$F[11]; $i++; }} else{{ $h{{$k}}.=",".$F[2]; $hs{{$k}}.=",".$F[4]; $he{{$k}}.=",".$F[5]; $hqs{{$k}}.=",".$F[15]; $hqe{{$k}}.=",".$F[16]; $hql{{$k}}.=",".$F[21]; $hid{{$k}}.=",".$F[11]; }} END{{ foreach(@o){{ print "$_\\t$h{{$_}}\\t$hs{{$_}}\\t$he{{$_}}\\t$hqs{{$_}}\\t$hqe{{$_}}\\t$hql{{$_}}\\t$hid{{$_}}\n" }} }}' | perl -naF'\\t' -e 'chomp($F[10]); @s=split(",",$F[4]); chomp; print join("\\t",@F[0..8])."\\t$s[0]\\t".@s."\\t$F[9]\\t$F[10]\n";' | sort -k1,1d -k10,10n -k11,11nr -k4,4nr > {output}
        """
