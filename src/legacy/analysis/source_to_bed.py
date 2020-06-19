import gzip
from Bio import SeqIO

def read_fasta(fasta):
    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rb") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fasta, "fasta"))
        
    return dict(zip([str(r.id) for r in records], [str(r.seq) for r in records]))

SOURCE_FA = "/media/scott/Rotom/assembly_data/hg002/hybrid_source.fa"
source = read_fasta(SOURCE_FA)

OUT_BED = "/media/scott/Rotom/assembly_data/hg002/hybrid_source.bed"
bed = open(OUT_BED, "w")
tab = "\t"
sourceMap = {'r' : 'canu', 'q' : 'supernova', "N" : "NNN"}

scoreMap = {'r' : '255', 'q' : '600', "N" : "999"}

for tigId, seq in source.items():
    print(tigId)
    null = "_"
    base = null
    start = 0
    
    for i,x in enumerate(seq):
        if base == x: continue
        if base == null: 
            base = x
            continue                

        bed.write(tab.join([tigId, str(start), str(i-1), sourceMap[base], scoreMap[base], "."]) + "\n")
        start, base = i, null
        
    bed.write(tab.join([tigId, str(start), str(i), "+", sourceMap[base], "."]) + "\n")
    