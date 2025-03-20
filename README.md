![hybrid image](https://github.com/ScottMastro/hybrid-pipeline/blob/master/hybrid.svg)
# hybrid-pipeline 

Prerequirements:

- Query assembly: A supernova assembly constructed from 10X Genomics linked-reads
- Reference assembly: A canu assembly constructed from PacBio CLR reads

In theory, other assemblies could be used. But the pipeline was created under the assumption that the query is accurate and highly fragmented, while the reference assembly is capable of resolving the low-complexity regions.


https://zenodo.org/records/15059067



## Step 0: purge_dups (optional):

A step to remove haplotigs from the assembly prior to the hybrid process. Here, the PacBio CLR reads will be aligned to a reference and duplicated content will be identified and removed.

Ensure `minimap2`, `samtools` and the `purge_dups` bin directory (with all the subtools) are in the PATH.

`snakemake -s Snakefile_purge_dup --config purgereads=/path/to/reads/ fasta=/path/to/assembly.fa`

The files provided in the `/example` directory have already been purged.

## Step 1: BLASTn alignments






snakemake -s ../snakemake/Snakefile_hybrid --config env=set_env_hybrid.sh r=HG002.contigs.fasta.purged.fa.gz q=NA24385_supernova.pseudohap.fasta.gz.purged.fa.gz \
 out=${DIR}/hybrid purgereads=${DIR}/pacbio/reads unitigs=$PBUNITIGS  \
        scripts=${BASEDIR}/../hybrid-pipeline/  
