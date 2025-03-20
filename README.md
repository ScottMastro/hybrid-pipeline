![logo](https://github.com/ScottMastro/hybrid-pipeline/blob/master/images/hybrid.svg)

## Required data:

- Query assembly: A **supernova assembly** constructed from 10X Genomics linked-reads
- Reference assembly: A **canu assembly** constructed from PacBio CLR reads

In theory, other assemblies could be used. But the pipeline was created under the assumption that the query is accurate and highly fragmented, while the reference assembly is capable of resolving the low-complexity regions.

![pipline details](https://github.com/ScottMastro/hybrid-pipeline/blob/master/images/pipeline.svg)

# Example data

Supernova and canu assemblies for sample HG002 can be found [at this zenodo record](https://zenodo.org/records/15059067). These are already purged.

Intermediate files and final hybrid assembly are also provided.

# Full Pipeline Details 

## Step 0: purge_dups (optional):

A step to remove haplotigs from the assembly prior to the hybrid process. This is to produce a strictly haploid assembly. Here, the PacBio CLR reads are aligned to both assemblies. Duplicated content will be identified and removed.

Ensure `minimap2`, `samtools` and the [`purge_dups`](https://github.com/dfguan/purge_dups) bin directory (with all the subtools) are in the PATH.

This Snakemake file defines all the steps:

```
snakemake -s hybrid-pipeline/scripts/purge_dup.snakefile --config purgereads=/path/to/reads fasta=/path/to/assembly.fa
```

## Step 1: BLASTn alignments

This step splits the query assembly into 1 kp chunks, creates a BLAST database for the reference assembly, parallelizes BLAST searches, and aggregates results in a single TSV.

Ensure `blastn` and `makeblastdb` are available.

```
snakemake -s hybrid-pipeline/scripts/blast_chunks.snakefile --config out=output_dir q=HG002_supernova.pseudohap.purged.fa.gz r=HG002_canu.contigs.purged.fa.gz
```

## Step 2: Hybrid assembly steps (main script)

> Note: supplying canu unitig data is optional but helps during scaffolding
A script to properly format the BED file is provided:
`python hybrid-pipeline/scripts/clean_bed.py canu.unitigs.bed canu.unitigs.clean.bed`

```
python hybrid-pipeline/src/hybrid.py --confident HG002_canu.unitigs.clean.bed -o output_dir {input.blocks} HG002_supernova.pseudohap.purged.fa.gz HG002_canu.contigs.purged.fa.gz
```
The output is a haploid hybrid assembly that is a mix of both haplotypes. Subsequent steps can be used to phase contigs.
