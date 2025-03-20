#!/usr/bin/env python2

# Enum for CIGAR operations (pysam convention)
MATCH, INS, DEL = 0, 1, 2
SHOW_MATCH_WARNING = False
REF_SKIP, SOFT_CLIP, HARD_CLIP = 3, 4, 5
PAD, EQUAL, DIFF = 6, 7, 8

# Groups of related CIGAR operations
QUERYOPS = (MATCH, EQUAL, DIFF, INS) # operations that include query bases
REFOPS = (MATCH, EQUAL, DIFF, DEL)   # operations that include reference bases
CLIPPINGOPS = (HARD_CLIP, SOFT_CLIP) # operations that clip query bases

CHROM_ = "#chrom"
POS_ = "pos"
TOTAL_ = "alignmentCount"
DISCORD_ = "diffCount"
DIFF_ = "substitutionCount"
DEL_ = "deletionCount"
INS_ = "insertionCount"
DOM_ = "dominantType"
SINGLE_ = "singletonCount"
PC_DISCORD_ = "differenceProportion"
PC_SINGLE_ = "singletonProportion"
PC_DOM_ = "dominantProportion"

colourMap = {DIFF:'green', INS:'blue', DEL:'red', None:'black'}

def add_columns(df):
    
    df[DOM_] = None
    df[PC_DISCORD_] = 1.0*df[DISCORD_] / df[TOTAL_]
    df[PC_SINGLE_] = 1.0*df[SINGLE_] / df[DISCORD_]

    df[PC_DOM_] = df[[DIFF_, DEL_, INS_]].max(axis = 1) / (df[DISCORD_])

    def update(g, g2, g3, value): 
        df.loc[(df[g] > df[g2]) & (df[g] > df[g3]), DOM_] = value

    update(DIFF_,  DEL_, INS_, DIFF)
    update( INS_, DIFF_, DEL_, DEL)
    update( DEL_, DIFF_, INS_, INS)
    
    return df
