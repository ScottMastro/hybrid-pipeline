#!/bin/env python3

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from random import randint
import os
import argparse

# Parse the input arguements:
parser = argparse.ArgumentParser(description = "BLAT intervals of contig against a reference")
parser.add_argument("-contig", help = "fasta file of the contig (one contig per file)")
parser.add_argument("-n", help = " number of nucleotides per query to BLAT")
parser.add_argument("-outdir", help = "the path to the directory for writing output files")
args = parser.parse_args()

n = int(args.n)
contig_base = os.path.basename(args.contig)

# Select sections of each contig and saves it in a separate .fasta file


seq_list = []

with open(args.contig, "r") as contigs:
    for seq in SeqIO.parse(contigs, "fasta"):
        header = seq.id
        nt = seq.seq
        nt_len = len(nt)
        
        if nt_len < n:
            print("Contig length shorter than n")
            break
    
        # number of stretches of nt to test (10% coverage if no overlap)
        l = int((len(nt)/n)/10)
        ni = 0
        
        # Minimun sections per contig = 2
        if l == 0:
            print(" contig: {} short, only 2 sections to BLAT".format(header))
            ind = randint(0,(nt_len-n))
            q_nt = nt[ind:ind+n]
            record = SeqRecord(q_nt, 
                             id = header + "|" + str(nt_len) + "_" + str(ind) + "\n",
                             name = header,
                             description = "Contig location {} {}".format(str(ni), str(ni+n)))
            
            seq_list.append(record)            
        
       
        # for 10% coverage (assuming no overlap)
        for x in range(0, l+1):
            x_str = str(x)
            ind = randint(0,(nt_len-n))
                
            q_nt = nt[ind:(ind + n)]
            record = SeqRecord(q_nt, 
                             id = header + "|" + str(nt_len) + "_" + str(ind) + "\n",
                             name = header,
                             description = "Contig location {} {}".format(str(ind), str(ind+n)))
            
            seq_list.append(record)            
                
                
# File where the chopped up contigs are saved
# Separate sections ("chops") into different fasta files
# fa_per_file = 1/(n^2)*(10^10)
# if n = 1000, fa_per_file = 10 000
# if n = 250, fa_per_file = 160 000
fa_per_file = int((1/n**2)*(10**10))
i = 0
for fasta in range(0, int(len(seq_list)/fa_per_file)+1):
    chop_dir = args.outdir + "/" + str(contig_base[:contig_base.index(".")]) + "_chops_{}.fa".format(fasta)
    chop = open(chop_dir, "w+")    
    print("Writing to {}".format(chop_dir))
    SeqIO.write(seq_list[i:i+fa_per_file], chop_dir, "fasta")
    i += fa_per_file
    chop.close()
    
    
        



            

            
            



