#!/usr/bin/env python2
"""Plot concordance of read alignments to a genome"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
import gzip
from Bio import SeqIO
from collections import Counter
import concordance_helper as ch

def read_fasta(fasta):
    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rb") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fasta, "fasta"))
        
    return dict(zip([str(r.id) for r in records], [str(r.seq) for r in records]))

def read_data(path):
    return pd.read_csv(path, sep=",", header=0)


def plot_confusion(df):
    
    covDf = df[df[ch.TOTAL_] <= 50]
    
    sb.pairplot(covDf[[ch.DIFF_, ch.DEL_, ch.INS_]],
                plot_kws={'alpha':0.1})    

def plot_discorant_hist(df, minDepth=10, maxDepth=100, minDiscordance=0, bins=50, singletonFilter=None):
    
    a, b = 0.4, bins
    percentFilter=0.9  
    subDf = df[(df[ch.TOTAL_] > minDepth) & (df[ch.TOTAL_] < maxDepth) &  \
               (df[ch.PC_DISCORD_] > minDiscordance)] 
    if singletonFilter is not None:
        if singletonFilter:
            subDf = subDf[(subDf[ch.PC_SINGLE_] >= 0.5)] 
        if not singletonFilter:
            subDf = subDf[(subDf[ch.PC_SINGLE_] < 0.5)] 

    fig, axs = plt.subplots(3, sharey=True, gridspec_kw={"hspace": 0})
    
    def filterGroup(g, pc=0): 
        return subDf.loc[(subDf[ch.DOM_] == g) & (subDf[ch.PC_DOM_] >= pc), ch.PC_DISCORD_]

    axs[0].set(title="Discordance Histogram")
    axs[-1].set(xlabel="Percent Discordance")

    axs[0].hist(filterGroup(ch.DIFF), color=ch.colourMap[ch.DIFF], label="DIFF", alpha=a, bins=b)
    axs[0].hist(filterGroup(ch.DIFF, percentFilter), color=ch.colourMap[ch.DIFF], label="DIFF_", alpha=a, bins=b)
    axs[0].set(ylabel="SUBSTUTION Frequency")

    axs[1].hist(filterGroup(ch.DEL), color=ch.colourMap[ch.DEL], label="DEL", alpha=a, bins=b)
    axs[1].hist(filterGroup(ch.DEL, percentFilter), color=ch.colourMap[ch.DEL], label="DEL_", alpha=a, bins=b)
    axs[1].set(ylabel="DELETION Frequency")

    axs[2].hist(filterGroup(ch.INS), color="b", label=ch.colourMap[ch.INS], alpha=a, bins=b)
    axs[2].hist(filterGroup(ch.INS, percentFilter), color=ch.colourMap[ch.INS], label="INS_", alpha=a, bins=b)
    axs[2].set(ylabel="INSERTION Frequency")

    fig.set_size_inches(8,12 , forward=True)
    fig.tight_layout()
    plt.show()
    
def plot_coverage(df, maxDepth=80, minDiscordance=0, bins=50):    
    ceilDf = df
    ceilDf[ch.TOTAL_] = df[ch.TOTAL_].apply(lambda x : min(x, maxDepth))
    plt.hist(ceilDf[ch.TOTAL_], bins=50)
    

    
    
def position_distribution(df, nchrom=1, minDepth=10, maxDepth=100, minDiscordance=0):
    
    a = 0.1
    fig, axs = plt.subplots(1, nchrom, sharey=True, gridspec_kw={'wspace': 0})
    if nchrom == 1: axs = [axs]
    
    subDf = df[(df[ch.TOTAL_] > minDepth) & (df[ch.TOTAL_] < maxDepth) &  \
               (df[ch.PC_DISCORD_] > minDiscordance)]

    
    def filterChrom(chromId): 
        return subDf[(subDf[ch.CHROM_] == chromId)]
    
    for i,ax in enumerate(axs):
        ax.set(title="Contig " + str(i))

        chromDf = filterChrom(i)
        ax.scatter(chromDf[ch.POS_], chromDf[ch.PC_DISCORD_], alpha=a, \
                   c=chromDf[ch.DOM_].apply(lambda x: ch.colourMap[x]))  #colour by type
                   #c=chromDf[TOTAL_]) #colour by depth


    fig.set_size_inches(26,16 , forward=True)
    fig.tight_layout()
    plt.show()

def position_distribution2(df, chromId=0, start=0, end=6e6, minDiscordance=0.9):
    
    a = 0.3
    fig, axs = plt.subplots(3, sharex=True, gridspec_kw={"hspace": 0})    
    subDf = df[(df[ch.CHROM_] == chromId) & 
             (df[ch.POS_] >= start) & (df[ch.POS_] <= end) &
                     (df[ch.PC_DISCORD_] > minDiscordance)]
    
    def filterGroup(g): 
        return subDf[(subDf[ch.DOM_] == g)]
    
    axs[0].set(title="Contig " + str(chromId))
    
    diffDf = filterGroup(ch.DIFF)
    delDf = filterGroup(ch.DEL)
    insDf = filterGroup(ch.INS)

    axs[0].scatter(diffDf[ch.POS_], diffDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.DIFF], label="DIFF", alpha=a)
    axs[1].scatter(delDf[ch.POS_], delDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.DEL], label="DEL", alpha=a)
    axs[2].scatter(insDf[ch.POS_], insDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.INS], label="INS", alpha=a)


    fig.set_size_inches(36,4 , forward=True)
    fig.tight_layout()
    plt.show()


def position_distribution3(df, chromId=0, window=100):
    
    a = 0.3
    subDf = df[(df[ch.CHROM_] == chromId)]
    end = subDf[ch.POS_].max()
    def count(start): 
        return subDf[(subDf[ch.POS_] >= start) & (subDf[ch.POS_] < start+window)].shape[0]
    
    y = [count(start) for start in range(0, 10000000, window)]
    x = [start for start in range(0, 10000000, window)]

    plt.scatter(x, y, alpha=a)

    fig, axs = plt.subplots(3, sharex=True, gridspec_kw={"hspace": 0})    
    
    def filterGroup(g): 
        return subDf[(subDf[ch.DOM_] == g)]
    
    axs[0].set(title="Contig " + str(chromId))
    
    diffDf = filterGroup(ch.DIFF)
    delDf = filterGroup(ch.DEL)
    insDf = filterGroup(ch.INS)

    axs[0].scatter(diffDf[ch.POS_], diffDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.DIFF], label="DIFF", alpha=a)
    axs[1].scatter(delDf[ch.POS_], delDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.DEL], label="DEL", alpha=a)
    axs[2].scatter(insDf[ch.POS_], insDf[ch.PC_DISCORD_],  color=ch.colourMap[ch.INS], label="INS", alpha=a)


    fig.set_size_inches(36,4 , forward=True)
    fig.tight_layout()
    plt.show()
    

def source_analysis(df, source):
    
    chrom="hybrid_0"
    
    def get_source(chrom, pos):
        return source[chrom][pos-1].upper()

    positions = np.array(df[ch.POS_])
    sources = [get_source(chrom, pos) for pos in positions]
    letter_counts = Counter(sources)
    df = pd.DataFrame.from_dict(letter_counts, orient='index')
    df.plot(kind='bar')



DIR = "/media/scott/Rotom/assembly_data/hg002/"
GENOME_CSV = DIR + "genome_concordance_ccs.csv"
SOURCE_FA = DIR + "hybrid_source.fa"
df = read_data(GENOME_CSV)
df = ch.add_columns(df)

plot_discorant_hist(df, singletonFilter=None, bins=200)
plot_discorant_hist(df, singletonFilter=True, bins=200)
plot_discorant_hist(df, singletonFilter=False, bins=200)

#source = read_fasta(SOURCE_FA)
