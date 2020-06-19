#!/usr/bin/env python3

import matplotlib.pyplot as plt 
import matplotlib.patches as patches
import pandas as pd
import math
import argparse

# Parse the input arguements:
parser = argparse.ArgumentParser(description = "counts the number of top BLAT matches to each chromosome")
parser.add_argument("-hits", help = "file containing the top hits for the contigs")
args = parser.parse_args()


# open the text file containing condensed information for each top hit
hit = open(args.hits, "r")
hits = pd.read_csv(hit, sep="\t", lineterminator="\n", names=None)

# 25 total chromosomes 
chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 
               19, 20, 21, 22, 'X', 'Y', 'Other']
# according to GRCh38.p12
chr_size = [248956422,242193529,198295559,190214555,181538259,170805979,
            159345973,145138636,138394717,133797422,135086622,133275309,
            114364328,107043718,101991189,90338345,83257441,80373285,
            58617616,64444167,46709983,50818468,156040895,57227415,
            10**7]

# dict of each chromosome and a corresponding colour
chr_col = { 1:"deeppink", 2:"magenta", 3:"darkviolet", 4:"mediumblue", 
            5:"steelblue", 6:"cornflowerblue", 7:"turquoise", 8:"green", 9:"darkseagreen", 
            10:"lime", 11:"gold", 12:"goldenrod", 13:"lightsalmon", 
            14:"tomato", 15:"red", 16:"brown", 17:"indianred", 
            18:"palevioletred", 19:"plum", 20:"mediumorchid", 
            21:"mediumvioletred", 22:"deeppink", 'X':"darkviolet",
            'Y':"mediumblue", 'Other':'black'}


genome_size = sum(chr_size) #3088286401
bil = 1*(10**9)
# height of rectangle:
height = 10**8 
gap_width = bil/50


# create a dictionary of the start position of each chromosome on the plot
x_start = 0
chr_start_pos = {}
for x_span, name in zip(chr_size, chromosomes):
    chr_start_pos[str(name)] = x_start
    x_start += (gap_width + x_span)

# midpoint coordinate of the contig, to have it centered 
contig_mid = x_start/2 
# min size of contig on plot, if 10nt (log10(10) = 1)
contig_base = 5*10**5 


#Function that creates a new plot and adds the human chromosomes
# parameters:
# contig = name of the contig, str
# contig_size = size of the contig, int
# target = name of the matched target, str or int
# pos_list = a list of list of the position on the contig, target, and name of target. 
#     eg. [[contig, target, target_name]...]
def genome_pos(contig, contig_size, contig_span, pos_list):
    # create new plot
    fig = plt.figure(figsize=(10,5), dpi=500)
    ax = fig.add_subplot(111, aspect='equal') 
    
    x_start = 0
    
    # add the rectangle for the genome
    for x_span, name in zip(chr_size, chromosomes): 
        if name == "Other":
            ax.add_patch(patches.Rectangle((3.6*10**9, bil), x_span, height, 
                                           color = chr_col[name]))
            plt.text(3.57*10**9, (bil+height*1.3), name, fontsize=5)
        else:
            ax.add_patch(patches.Rectangle((x_start, bil), x_span, height, 
                                           color = chr_col[name]))
            if (type(name) == int and name < 10) or name == 'X':
                mid_chr = x_start + x_span/2.6
                plt.text(mid_chr, (bil+height*1.3), name, fontsize=5)
                x_start += (gap_width + x_span)
            elif (type(name) == int and name < 17) or name == 'X':
                mid_chr = x_start + x_span/3
                plt.text(mid_chr, (bil+height*1.3), name, fontsize=5)
                x_start += (gap_width + x_span)
            else:
                plt.text(x_start, (bil+height*1.3), name, fontsize=5)
                x_start += (gap_width + x_span)
                
    contig_start = contig_mid - (contig_span/2)
    
    #contig rectangle (starts =  897071600.25, spans = 1794143200.5)
    ax.add_patch(patches.Rectangle((contig_start, 0), contig_span, height, 
                                   color = "grey"))
    plt.text((contig_mid - contig_span/2), -(height/1.5), 
             'Contig: {} ({}nt)'.format(contig, contig_size), fontsize=5)

    plt.axis([0, genome_size*1.3, 0, genome_size/2.5])
    plt.axis('off')
    
    # thicker line width for plots with less lines for improved visibility 
    if len(pos_list) < 10:
        linew = 0.3
    elif len(pos_list) < 100:
        linew = 0.2
    else:
        linew = 0.1
    
    
    # Draw the lines and determine the dominant chr
    chr_list = {}
    for pos in pos_list:
        chr_line = pos[2]
        plt.plot(pos[:2], [height, bil], color = chr_col[chr_line], linewidth=linew) 
        
        if str(chr_line).isalpha():
            if chr_line == "X":
                chr_line = 23
            elif chr_line == "Y":
                chr_line = 24
            elif chr_line == "Other":
                chr_line = 25
            
        if chr_line not in chr_list:
            chr_list[chr_line] = 1
        else:
            chr_list[chr_line] += 1
                    
    sorted_chr_list = sorted( ((v,k) for k,v in chr_list.items()), reverse=True)
    
    dom_chr = sorted_chr_list[0][1]
    if dom_chr > 22:
        if dom_chr == 23:
            dom_chr = 'X'
        elif dom_chr == 24:
            dom_chr = 'Y'
        elif dom_chr == 25:
            dom_chr = 'Other'

    plt.savefig('chr{} {} - {}nt.png'.format(str(dom_chr), contig, 
                                                    contig_size), dpi = 500)
    plt.close()
    
    
        
# loop through all the rows of the data frame to create a list of positions on the contig and target
# then apply the information using the function genome_pos

prev_contig = ""
for index, row in hits.iterrows():
    contig_frag = row[9]
    contig = contig_frag[:contig_frag.index("|")]
    target_name = row[13]
    try:
        target = int(target_name[3:])
    except:
        target = target_name[3:]    
        if target not in chromosomes:
            target = 'Other'
        
    cname_end_ind = len(contig_frag) - contig_frag[::-1].index("_") -1 
    c_len = int(contig_frag[contig_frag.index("|")+1 : cname_end_ind])
    c_span = math.log10(c_len)**4 * contig_base
    if c_span > genome_size:
        c_span = genome_size
    c_start = int(contig_frag[cname_end_ind + 1:])
    c_spt = (c_start/c_len)*c_span + (contig_mid - c_span/2)
    
    g_spt = row[15] + chr_start_pos[str(target)]
    if target == "Other" :
        g_spt = 3.61*10**9

    
    if index == 0:       #first row
        pos_list = []
        
    elif contig != prev_contig: #next contig
        # create plot using the function genome_pos
        genome_pos(prev_contig, prev_contig_size, prev_contig_span, pos_list)
        pos_list = []
    
    elif index == len(hits)-1:  #last row
        pos_list.append([c_spt, g_spt, target])
        prev_contig = contig
        genome_pos(prev_contig, prev_contig_size, prev_contig_span, pos_list)
    
    pos_list.append([c_spt, g_spt, target])
    prev_contig = contig
    prev_contig_size = c_len
    prev_contig_span = c_span
    
    
hit.close()
