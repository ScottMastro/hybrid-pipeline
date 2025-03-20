import sys
import gzip
i=sys.argv[1]
o=sys.argv[2]

#the 16-base 10x barcode plus 7 additional bases
def trim_10x_barcode(fastq, trimmedFq, trim=23):
    reader = gzip.open(fastq, "rt" if fastq.endswith("gz") else "r")
    writer = gzip.open(trimmedFq, 'wt')
    
    line = reader.readline()
    while(line):        
        if not line.startswith("@") and not line.startswith("+"):
            line = line[trim:]            
        writer.write(line)
        line = reader.readline()
    reader.close()
    writer.close()

trim_10x_barcode(i,o)
