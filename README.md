![logo](https://github.com/ScottMastro/hybrid-pipeline/blob/master/images/hybrid.svg)

## Required data:

- Query assembly: A **supernova assembly** constructed from 10X Genomics linked-reads
- Reference assembly: A **canu assembly** constructed from PacBio CLR reads

In theory, other assemblies could be used. But the pipeline was created under the assumption that the query is accurate and highly fragmented, while the reference assembly is capable of resolving the low-complexity regions.

# Example data

Supernova and canu assemblies for sample HG002 can be found [at this zenodo record](https://zenodo.org/records/15058963). These are already purged.

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

This step splits the query assembly into 1 kb chunks, creates a BLAST database for the reference assembly, parallelizes BLAST searches, and aggregates results in a single TSV.

Ensure NCBI BLAST+ tools (`blastn` and `makeblastdb`) are installed and available.

```
snakemake -s hybrid-pipeline/scripts/blast_chunks.snakefile --config out=output_dir q=HG002_supernova.pseudohap.purged.fa.gz r=HG002_canu.contigs.purged.fa.gz
```

## Step 2: Hybrid assembly 

The hybrid assembly process consists of three key steps: stitching alignments, welding blocks, and scaffolding paths.

### Stitching alignments
Using BLAST alignments, the pipeline identifies continuous chunk alignments to construct blocks—segments of sequence shared between the two assemblies. These blocks are further grouped into megablocks, which capture larger structural similarities across assemblies.

### Welding blocks
To precisely align block ends, pairwise alignment is used to establish shared positions between assemblies. This process identifies "forks", points where one assembly can transition into the other at consistent positions. These forks define a navigable path between assemblies, helping to fill gaps and correct local misassemblies.

### Scaffolding paths
In the final step, the identified paths are connected to form scaffolds, effectively linking contigs into longer scaffolds.

> Note: supplying canu unitig data is optional but helps during scaffolding. A script to properly format the BED file is provided: `python hybrid-pipeline/scripts/clean_bed.py canu.unitigs.bed canu.unitigs.clean.bed`

```
python hybrid-pipeline/src/hybrid.py --confident HG002_canu.unitigs.clean.bed -o output_dir HG002.blastn.summary.txt HG002_supernova.pseudohap.purged.fa.gz HG002_canu.contigs.purged.fa.gz
```
The output is a haploid hybrid assembly that is a mix of both haplotypes. 

> Note: We can assign each scaffold to a human chromosome by finding the best alignment against the T2T CHM13 assembly.

## Step 3: Phasing

The phasing pipeline is a complex multi-step process that integrates multiple tools, tailored specifically to the datasets used in this project. The Snakemake workflows and supporting scripts can be found in `hybrid-pipeline/polish/*.snakefile` and serve as documentation of the methodology, but likely require modification to run.

First, `pbmm2` is used to align the original PacBio reads back to the hybrid assembly and [`SDA`](https://github.com/mrvollger/SDA) is used to identify and adjust for segmental duplications.

The process below is done on a per-scaffold basis:

The following steps are performed iteratively to refine the assembly:

- Align PacBio reads using `pbmm2`.
- Polish with `gcpp` (Arrow algorithm).
- Align 10X Genomics reads using `bwa`.
- Polish with `pilon`.

From this process, a consensus assembly is created. Reads are again aligned against this consensus, `longranger` (10XG) and `longshot` (PacBio) are used to call variants. High-confidence heterozygous variants (shared between the two callsets) are identified, `whatshap` is used to phase all the reads.

After read phasing, the iterative polishing process is repeated using haplotype-specific reads. This step refines each haplotype separately and can be repeated multiple times until convergence. The final output consists of two haplotype-resolved scaffolds.

The complete phasing pipeline is visually represented below:

![pipeline details](https://github.com/ScottMastro/hybrid-pipeline/blob/master/images/pipeline.svg)
