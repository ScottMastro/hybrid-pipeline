#!/usr/bin/env python2
"""Compute concordance statistics of read alignments to a genome"""

from __future__ import division
import pandas as pd
import concordance_helper as ch

def read_data(path):
    return pd.read_csv(path, sep=",", header=0)

def add_umns(df):
    
    df[ch.DOM_] = None
    df[ch.PC_DISCORD_] = 1.0*df[ch.DISCORD_] / df[ch.TOTAL_]
    df[ch.PC_SINGLE_] = 1.0*df[ch.SINGLE_] / df[ch.DISCORD_]

    df[ch.PC_DOM_] = df[[ch.DIFF_, ch.DEL_, ch.INS_]].max(axis = 1) / (df[ch.DISCORD_])

    def update(g, g2, g3, value): 
        df.loc[(df[g] > df[g2]) & (df[g] > df[g3]), ch.DOM_] = value

    update(ch.DIFF_,  ch.DEL_, ch.INS_, ch.DIFF)
    update( ch.INS_, ch.DIFF_, ch.DEL_, ch.DEL)
    update( ch.DEL_, ch.DIFF_, ch.INS_, ch.INS)
    
    return df

def filter_( df, column, fn ):
    return df[df[column].apply(fn)]

def process(df, minDepth=10, maxDepth=80):

    df = ch.add_columns(df)
    subDf = df
    
    subDf = filter_( subDf, ch.PC_DISCORD_, (lambda x: x > 0.9) )
    subDf = filter_( subDf, ch.TOTAL_, (lambda x: x >= minDepth) )
    subDf = filter_( subDf, ch.TOTAL_, (lambda x: x <= maxDepth) )

    #source = read_fasta(SOURCE_FA)


    

    
    return df

#lambda x: x + 5

DIR = "/media/scott/HDD/sickkids/NA24385/"
GENOME_CSV = DIR + "genome_concordance_ccs.csv"
SOURCE_FA = DIR + "hybrid_source.fa"
df = read_data(GENOME_CSV)


