import argparse
import gzip

HEADERINFO="##vcf_correction=LongRanger 2.2.2 GT calls have been corrected\n" + \
            "##INFO=<ID=HetCorrected,Number=1,Type=Integer,Description=\"1 for variants that were corrected to het (LongRanger 2.2.2 bug)\">\n"
              
CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE=0,1,2,3,4,5,6,7,8,9

#==================================================
# This script attempts to correct a bug in Long Ranger 2.2.2
#==================================================

def _fix_line(line):
    cols = line.split("\t")
    
    try:
        
        sampDict = { f:s for f,s in zip(cols[FORMAT].split(":"), cols[SAMPLE].split(":")) }
        #infos = cols[INFO].split(";")
        #infoDict = { v.split("=")[0]:v for v in infos }

        gtCount = int(sampDict["GT"][0]) + int(sampDict["GT"][-1])
        pl = [ int(x) for x in sampDict["PL"].split(",") ]
                
        if gtCount > 1 and pl[1] < pl[2]:
            print("Correcting " + cols[CHROM] + ":" + cols[POS], 
                  "\n\tGT=" + sampDict["GT"], "PL=" + str(pl))
            
            fixedLine = [ cols[CHROM], cols[POS], cols[ID], cols[REF], cols[ALT], cols[QUAL], cols[FILTER] ]
            fixedLine.append(cols[INFO] + ";" + "HetCorrected=1")
            sampDict["GT"] = "0/1"
            fixedLine.append(cols[FORMAT])
            samp = [ sampDict[f] for f in cols[FORMAT].split(":") ]
            fixedLine.append(":".join(samp))
        else:
            return line
        fixedLine = "\t".join(fixedLine)
        return fixedLine
    except:
        return line  
    
def run_correction(inputVCF, outputVCF):
    
    #construct_barcode_dict(inputVCF)
    
    headerWritten = False
    
    if inputVCF.endswith(".gz"):
        reader = gzip.open(inputVCF, "rt")
    else:
        reader = open(inputVCF, "r")

    if outputVCF.endswith(".gz"):
        writer = gzip.open(outputVCF, 'wt')
    else:
        writer = open(outputVCF, 'w')


    for line in reader:        

        if len(line) > 1:
            
            if not headerWritten:
                if line[0:2] == "##":
                    writer.write(line)
                elif line[0] == "#":
                    writer.write(HEADERINFO)
                    writer.write(line)
                    headerWritten = True
                continue
                    
        correctedLine = _fix_line(line)
        writer.write(correctedLine)
        
    reader.close()
    writer.close()
    
    return outputVCF

if __name__== "__main__":
    parser = argparse.ArgumentParser(description="Longranger 2.2.2 VCF Correction Tool (single sample only)")
    parser.add_argument("vcf", metavar="vcf", default="phased_variants.vcf.gz", nargs="?", 
                    help="Longranger 2.2.2 \"phased_variants\" VCF" )
    args = parser.parse_args()

    inputVCF = args.vcf
    outputVCF = "/".join(inputVCF.split("/")[:-1]) + "/phased_variants.corrected.vcf"

    run_correction(inputVCF, outputVCF)
