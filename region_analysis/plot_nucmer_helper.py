import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
#from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import math
import re

#100 units for plotting space, 20 units for gene track
trackheight=100.0
geneheight=30.0
leftbuffer=0.05
gap=trackheight/50
tigheight=trackheight/8   #20

assemblyFontSize=24   #18
assemblyFontAlpha=0.7
contigFontSize=12 #8
geneFontSize=16 #10
mainFontSize=16 #10


snpscmap=plt.cm.gist_ncar
snpcolours = dict()
snpcolours["T"]= int(74*255/100) #ff3730
snpcolours["G"]= int(67*255/100) #d07004
snpcolours["C"]= int(11*255/100) #0000fd
snpcolours["A"]= int(33*255/100) #01ff02
snpcolours["."]= int(3*255/100) #8676d3 deletion
snpcolours["-"] = int(95*255/100) #ee00ff #insertion

def get_snp_colour(snp):

    #insertion
    if snp["refbase"] == "." :
        return snpcolours["-"]
    if snp["altbase"] in snpcolours:
        return snpcolours[snp["altbase"]]
    return 0

def snp_legend(ax):
    labs=["T", "C", "A", "G", "Insertion", "Deletion"]
    custom_lines = [Line2D([0], [0], color=snpscmap(snpcolours[labs[0]]), lw=4),
                    Line2D([0], [0], color=snpscmap(snpcolours[labs[1]]), lw=4),
                    Line2D([0], [0], color=snpscmap(snpcolours[labs[2]]), lw=4),
                    Line2D([0], [0], color=snpscmap(snpcolours[labs[3]]), lw=4),
                    Line2D([0], [0], color=snpscmap(snpcolours["-"]), lw=4),
                    Line2D([0], [0], color=snpscmap(snpcolours["."]), lw=4)]
    
    ax.legend(custom_lines, labs, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=6, fontsize=mainFontSize)

def get_totalwidth(files):
    a = pd.read_csv(files[0], sep='\t', header=None)
    return int(a[18].iloc[0])


def initialize_plot(fig, ax, tracknames, width, linewidth=2, colour='k', alpha=1):
    plt.ylim(0, trackheight + geneheight)
    plt.xlim(-width*leftbuffer, width)
    plt.box(False)
    ax.axes.get_yaxis().set_visible(False)
    #plt.title(title)
    fig.set_size_inches(12, 6)
    height=trackheight/len(tracknames)
    def makeTrackRectangle(index):
        return Rectangle( (0, (index)*height), width, height)

    tracks = [makeTrackRectangle(i) for i, name in enumerate(tracknames)]
    tk = PatchCollection(tracks, facecolor='none', alpha=alpha, edgecolor=colour, linewidth=linewidth)
    ax.add_collection(tk)
    
    for i, name in enumerate(tracknames):
        plt.annotate(str(name), xy=(width*0.01, i*height + height/4), size=assemblyFontSize, alpha=assemblyFontAlpha, clip_on=True)


def add_genes(ax, gff_file, totalwidth, source="ensembl", cmap="nipy_spectral", alpha=0.65):
    genes = pd.read_csv(gff_file, sep='\t', header=None)
    genes.columns = ["id", "source", "type", "start", "end", "x1", "direction", "x2", "info"]
    genes = genes[genes["type"].isin(["gene", "exon"])]
    genes = genes[genes["source"].str.contains(source)]    

    if(len(genes) < 1):
        print("WARNING: no genes or exons in gff file")
        
    def extract_value(info, value):
        info = re.sub(';\s', ';', info)
        if info[-1] is ";":
            info = info[:-1]
        info=re.sub('[\"]', '', info).split(";")
        infodict=dict(zip([x.split(" ")[0] for x in info], [x.split(" ")[1] for x in info]))
        if value in infodict:
            return str(infodict[value])
        else:
            return "NA"
            
    genes["gene_name"] = [extract_value(x["info"], "gene_name") for i, x in genes.iterrows()]

    #annotate gene names on plot and compute gene colour
    colourdict=dict()
    for index, gene in genes[genes["type"] == "gene"].iterrows():
    
        if(gene["direction"] == "+"): 
            ypos = trackheight + geneheight/4
            halign="left"
            xpos = gene["start"]

        else: 
            ypos = trackheight + (3*geneheight/4)
            halign="right"
            xpos = gene["end"]
            
        plt.annotate(gene["gene_name"], xy=(xpos, ypos), horizontalalignment=halign,
           verticalalignment="bottom", size=geneFontSize,  clip_on=False)
        
        colourdict[gene["gene_name"]] = math.sqrt((gene["start"] + gene["end"])/2)
        
    def make_rectangle(row):
        left=min(row["start"], row["end"])
        right=max(row["start"], row["end"])
        if(row["type"] == "exon"): height=geneheight/6
        else: height=(geneheight/4)/6
        if(row["direction"] == "+"): ypos = trackheight + geneheight/8 - height/2
        else: ypos = trackheight + (5*geneheight/8) - height/2
        rect = Rectangle( (left, ypos), right - left, height)
        colour = colourdict[row["gene_name"]]
        return (rect, colour)
            
    #annotate gene rectangles on plot
    generects, colours=zip(*[make_rectangle(row) for index, row in genes.iterrows()])
    if(len(generects) > 0):
        pc = PatchCollection(generects, cmap=cmap, alpha=alpha)
        pc.set_array(np.array(colours))
        pc.set_clim([0, math.sqrt(totalwidth)])
        ax.add_collection(pc)

def read_alignments(align_file, lenFilter=-1, idFilter=-1):
    
    alignments = pd.read_csv(align_file, sep='\t', header=None)
    alignments.columns = ["contig", "date", "tiglength", "nucmer", "dir", "region", "tigstart", "tigend", "refstart", "refend", 
                                     "identity", "similarity", "length", "x1", "x2", "x3", "x4", "direction", "reflength", "x5", "x6"]

    if lenFilter > 0:
        alignments = alignments[alignments["length"] >= lenFilter]
    if idFilter > 0:
        alignments = alignments[alignments["identity"] >= idFilter]
        
    lengths = {tig:np.sum(alignments[alignments["contig"]==tig]["length"]) for tig in np.unique(alignments["contig"])}
    ranks = dict()
    r=0
    for key in sorted(lengths, key=lambda k: lengths[k]):
        ranks[key] = r
        r=r+1
    
    alignments["rank"] = [ranks[row["contig"]] for idx, row in alignments.iterrows()]
    return alignments
    
def add_alignments(ax, alignments, index, ntracks, width, snp_file="", 
                   colour="grey", alpha=0.6):
       
    for tig in np.unique(alignments["contig"]):
        tigname=str(tig).split("_")[0].split(" ")[0]
        rank = alignments[alignments["contig"] == tig]["rank"].iloc[0]
        ypos=(index+1)*trackheight/ntracks - tigheight*(rank+1) - gap*(rank+1)
        plt.annotate(tigname, xy=(-width*0.001, ypos), size=contigFontSize,
             verticalalignment="bottom", horizontalalignment="right")

    def make_rectangle(row):
        rank = row["rank"]
        ypos=(index+1)*trackheight/ntracks - tigheight*(rank+1) - gap*(rank+1)

        left=min(row["refstart"], row["refend"])
        right=max(row["refstart"], row["refend"])
        
        return Rectangle( (left, ypos), right - left, tigheight)
    
    rects = [make_rectangle(row) for idx, row in alignments.iterrows()]
    pc = PatchCollection(rects, facecolor=colour, alpha=alpha, edgecolor="k")
    ax.add_collection(pc)
                
    
    #snps
    if len(snp_file) > 0:
        
        def make_snp_rectangle(snp, rank):
            
            ypos=(index+1)*trackheight/ntracks - tigheight*(rank+1) - gap*(rank+1)
            left=snp["refpos"] - 0.5
            right=snp["refpos"] + 0.5
            rect = Rectangle( (left, ypos), right - left, tigheight)
            colour = get_snp_colour(snp)
            return (rect, colour)

        
        snplist = pd.read_csv(snp_file, sep='\t', header=None)
        snplist.columns = ["refpos", "refbase", "altbase", "tigpos", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "contig"] 

        for tig in np.unique(alignments["contig"]):
            snps = snplist[snplist["contig"] == tig]
            rank = alignments[alignments["contig"] == tig]["rank"].iloc[0]

            rects, colours=zip(*[make_snp_rectangle(row, rank) for idx, row in snps.iterrows()])
            pc_ = PatchCollection(rects, cmap=snpscmap, alpha=1, edgecolor="face")
            pc_.set_array(np.array(colours))
            pc_.set_clim([0, 255])

            ax.add_collection(pc_)

            snp_legend(ax)

def adjust_xaxis(ax, startPos):
    
    ax.xaxis.tick_top()  
    ax.axhline(y=trackheight+geneheight, xmin=leftbuffer, linewidth=4, color='grey')
    ax.tick_params(axis='x', colors='grey')
    ax.spines['top'].set_visible(True)
    labels = plt.xticks()[1]
    ticLabs = [str(int(x.get_text().replace("−", "-")) + startPos) for x in labels]
    ax.set_xticklabels(ticLabs, fontsize=mainFontSize)     

