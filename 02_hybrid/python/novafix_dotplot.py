import flexidot as dot
import glob, os, re
from Bio.Seq import Seq
from Bio import pairwise2

fixedDir = "../out/nfixed.fa"
canuDir = "../out/nfixed_canu.fa"
novaDir = "../out/nfixed_nova.fa"

f = open(fixedDir, "r")
fixedLines = f.readlines()[::-1]
f.close()

f = open(canuDir, "r")
canuLines = f.readlines()[::-1]
f.close()

f = open(novaDir, "r")
novaLines = f.readlines()[::-1]
f.close() 

tempfile = "./temp.fa"
counter=0
wordsize=11
while(len(fixedLines) >= 0):
    fixedLines.pop()
    canuLines.pop()
    novaLines.pop()

    stitch=fixedLines.pop()
    nova=novaLines.pop().replace("N", "")
    canu=canuLines.pop() 
    temp_lines = [">canu " + canu[0:15] + "\n", canu,
                  ">nova " + nova[0:15] + "\n", nova,
                  ">stitch " + stitch[0:15] + "\n", stitch]
    
    f = open(tempfile, "w")
    f.writelines(temp_lines)
    f.close()
    counter=counter+1
    prefix="../out/dotplot/" + str(counter) + "_"
    dot.main(tempfile, wordsize, prefix=prefix)
    file1=prefix + "-Polydotplot_LCS_Shading_Legend_*.png"
    file2=prefix + "-Polydotplot_lcs_data_file_6shades_ref0_ori0.txt"
    file3=prefix + "-Polydotplot_wordsize" + str(wordsize) + "_6shades_ref0_ori0.png"
    for x in glob.glob(file1): os.remove(x)
    os.rename(file2, prefix+"data.txt")
    os.rename(file3, prefix+"plot.png")
    
    
    
    
    
 
    
    
    
def trim(aln, base_buffer):
    canu_aln = aln[0][0].rstrip()
    nova_aln = aln[0][1].rstrip()
    
    while(nova_aln[0] == "-" or canu_aln[0] == "-"):
        canu_aln=canu_aln[1:]
        nova_aln=nova_aln[1:]

    while(nova_aln[-1] == "-" or canu_aln[-1] == "-"):
        canu_aln=canu_aln[:-1]
        nova_aln=nova_aln[:-1]
        
    nidx = nova_aln.index("N")
    counter = base_buffer
    left = nova_aln[:nidx]
    left_buffer = ""
    index = nidx-1
    while(index >=0 and counter > 0):
        if not canu_aln[index] == "-":
            left_buffer = canu_aln[index] + left_buffer
            counter=counter-1
            
        left=left[:-1]
        index=index-1
    
    left = left + left_buffer
    
    counter = base_buffer
    right = nova_aln[nidx:]
    right_buffer = ""
    index = nidx-1
    while(index < len(canu_aln) and counter > 0):
        if not canu_aln[index] == "-":
            right_buffer = right_buffer +  canu_aln[index]
            counter=counter-1
            
        right=right[1:]
        index=index+1
        
    right = right_buffer + right
    trimmed = (left + right).replace("N", "").replace("-", "")
    return trimmed
        
    
def stitch_middle(aln, pre, post):
    base_buffer=7
    
    pre = pre[:-base_buffer]
    post = post[base_buffer:]

    canu_aln = aln[0][0].rstrip()
    nova_aln = aln[0][1].rstrip()

    while(len(pre) > 0):
        if pre[0] == nova_aln[0]:
            pre = pre[1:]
    
        canu_aln = canu_aln[1:]
        nova_aln = nova_aln[1:]
        
    while(len(post) > 0):
        if post[-1] == nova_aln[-1]:
            post = post[:-1]
    
        canu_aln = canu_aln[:-1]
        nova_aln = nova_aln[:-1]

    canu_aln = canu_aln.replace("-", "")
    
    if(len(canu_aln) <= 0 ):
        return trim(aln, base_buffer)
    return canu_aln


def do_alignment(num):
    f = open(fixedDir, "r")
    fixedLines = f.readlines()[::-1]
    f.close()
    
    f = open(canuDir, "r")
    canuLines = f.readlines()[::-1]
    f.close()
    
    f = open(novaDir, "r")
    novaLines = f.readlines()[::-1]
    f.close() 

    for i in range(num):
        fixedLines.pop()
        canuLines.pop()
        novaLines.pop()
    
        stitch = fixedLines.pop().rstrip()
        novaSeq=novaLines.pop().rstrip()
        canuSeq=canuLines.pop().rstrip()
        
    pre = novaSeq[0:novaSeq.index("N")]
    post = novaSeq[novaSeq.rindex("N")+1:len(novaSeq)]
    

    aln = pairwise2.align.globalms(canuSeq, novaSeq, 2, -3, -5, -2, one_alignment_only=True)
    canuRev=str(Seq(canuSeq).reverse_complement())
    alnRev = pairwise2.align.globalms(canuRev, novaSeq, 2, -3, -5, -2, one_alignment_only=True) 


    if(aln[0][2] > alnRev[0][2]):
        stitched = stitch_middle(aln, pre, post)
    else:
        stitched = stitch_middle(alnRev, pre, post)
        
    

                