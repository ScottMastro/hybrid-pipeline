from pybedtools import BedTool
import argparse
import re

def main():
    
    parser = argparse.ArgumentParser(description="Fixes Canu unitig BED File")
    #Positional
    parser.add_argument("input", metavar="IN_BED", nargs="?", 
                        help="Canu unitig BED file" )
    parser.add_argument("output", metavar="OUT_BED", nargs="?",
                        help="Output cleaned BED file")

    args = parser.parse_args()
    bedFile = args.input
    outFile = args.output

    writer = open(outFile, "w")
    reader = open(bedFile, "r")

    for line in reader: 
        columns = line.split('\t')
        s,e = int(columns[1]), int(columns[2])
        if not s < e: s,e = e,s
        name = re.sub("ctg", "tig", columns[0])
        
        cleaned = '\t'.join([name, str(s), str(e)] + columns[3:])
        writer.write(cleaned)
        
    reader.close()
    writer.close()
    
    #try to open...
    bed = BedTool(outFile)

        
if __name__== "__main__":
  main()
  exit()
