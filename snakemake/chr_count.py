import sys
import pandas as pd

paf_file = sys.argv[1]
out_file = sys.argv[2]

chr = dict()

with open(paf_file) as f:
   for line in f:
      split = line.split("\t")
      if split[5] not in chr:
         chr[split[5]] = 0
      chr[split[5]] += int(split[8])-int(split[7]) 

chrlist = [(k, str(chr[k])) for k in sorted(chr, key=chr.get, reverse=True)]


with open(out_file, "w") as outfile:
   for c in chrlist:
      "\t".join(c)
      outfile.write("\t".join(c) + "\n")
