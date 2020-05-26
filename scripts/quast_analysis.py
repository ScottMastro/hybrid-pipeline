import pandas as pd
import numpy as np

import argparse

def main():

    DEFAULT="/media/scott/Zapdos/out/quast/quast_stats.txt"
    
    parser = argparse.ArgumentParser(description="QUAST Summary")

    #Positional
    parser.add_argument("quast_stats", metavar="stats", default=DEFAULT, nargs="?",
                        help="CSV with QUAST statistics")

    args = parser.parse_args()
    #set parameters from user input
    stats = args.quast_stats
    quastDf = pd.read_csv(stats)
    ids = set(quastDf["id"])
    assemblies = set(quastDf["assembly"])

    def average(column, name=None, d=2):
        print("---------- " + column + " ----------")

        if name is None: name = column
        for assembly in assemblies:
            x = np.mean(quastDf[quastDf["assembly"]==assembly][column].astype(float)) 
            print(column + " mean:\t" + str(round(x,d)) + "\t" + assembly )
        print("")

        '''
        if name is None: name = column
        z=""
        for assembly in assemblies:
            x = np.mean(quastDf[quastDf["assembly"]==assembly][column].astype(float)) 
            z = z + "\t" + str(round(x,d))
        print(column +  z )
        '''
    average("total_contigs")
    average("total_span")
    average("largest")
    average("gc")

    average("n50")
    #average("na50")

    average("misassembly")
    average("scaffold_gap")
    average("unaligned_span")
    average("genome_fraction")
    average("dup_ratio", d=3)

    average("n_100kb", d=3)
    average("mismatches_100kb", d=3)
    average("indel_100kb", d=3)

if __name__== "__main__":
  main()
  #exit()
    
        
        
        
