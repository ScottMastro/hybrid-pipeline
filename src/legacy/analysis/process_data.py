import kmer_tools as kt
import os
import pandas as pd
import matplotlib.pyplot as plt

Q,R,H,HG = "query", "ref", "hybrid", "hg38"

def get_data_frame(prefix, baseDir):
    
    counter = 0
    info = dict()
    
    for file in os.listdir(baseDir):
        if file.startswith(prefix):
            try:
                txt = baseDir + "/" + file + "/info.txt"
                reader = open(txt, "r")
                header = reader.readline().strip().split("\t")
                data = reader.readline().split("\t")
                reader.close()
                
                for key,value in zip(header,data):
                    if key not in info:
                        info[key] = [None]*counter
                        
                    try:
                        info[key].append(float(value))
                    except:
                        info[key].append(value)

                counter += 1
            except:
                pass

    return pd.DataFrame.from_dict(info)
    
def analyze_blockpaths(param):

    print("Analyzing block data...")       
    df = get_data_frame("block_", param.OUTPUT_DIR)
    
    
    def len_weighted(name, x, fullName=None):
        if fullName is None : fullName=name + x
        return (df["length_" + x] * df[fullName]).sum() / df["length_" + x].sum() 
                
    totalLen = { x: df["length_" + x].sum() for x in [Q,R,H,HG] }
    totalGaps = { x: int(df["gaps_" + x].sum()) for x in [Q,R,H,HG] }
    totalN = { x: int(df["ambiguous_bases_" + x].sum()) for x in [Q,R,H,HG] }

    weightedGC = { x: len_weighted("gc_percent_", x) for x in [Q,R,H,HG] }

    
    def draw_kmer_pie(compressed=False):
        p = "compressed_" if compressed else ""

        sharedKmersTotal = df[p + "kmer_overlap_query_ref_hg38"].sum()
        qrKmersTotal = df[p + "kmer_overlap_query_ref"].sum()
        qhgKmersTotal = df[p + "kmer_overlap_query_hg38"].sum()
        rhgKmersTotal = df[p + "kmer_overlap_ref_hg38"].sum()
        qKmersTotal = df[p + "kmer_overlap_query"].sum()
        rKmersTotal = df[p + "kmer_overlap_ref"].sum()
        hgKmersTotal = df[p + "kmer_overlap_hg38"].sum()
        
        kmerTotal = sharedKmersTotal + qrKmersTotal + qhgKmersTotal + rhgKmersTotal + qKmersTotal + rKmersTotal + hgKmersTotal
        print(kmerTotal)
        kt.draw_pie3_counts(Q, R, H, sharedKmersTotal, qrKmersTotal, qhgKmersTotal, rhgKmersTotal,
                            qKmersTotal, rKmersTotal, hgKmersTotal, param.OUTPUT_DIR + p + "/kmer_total")
    
    draw_kmer_pie(compressed=False)
    
        
    qrWeightedKmerSimilarity = round(100*len_weighted(None, HG, fullName="kmer_similarity_query_ref"),4)
    rhgWeightedKmerSimilarity = round(100*len_weighted("kmer_similarity_ref_", HG),4)
    qhgWeightedKmerSimilarity = round(100*len_weighted("kmer_similarity_query_", HG),4)

    qrWeightedKmerSimilarity = round(100*len_weighted(None, HG, fullName="compressed_kmer_similarity_query_ref"),4)
    rhgWeightedKmerSimilarityCompressed = round(100*len_weighted("compressed_kmer_similarity_ref_", HG),4)
    qhgWeightedKmerSimilarityCompressed = round(100*len_weighted("compressed_kmer_similarity_query_", HG),4)

    
    weightedPBConcordance = { x: round(len_weighted("pacbio_concordance_", x),4) for x in [Q,R,H,HG] }


    rHybridComposition = len_weighted(None, H, fullName="hybrid_r_composition")
    qHybridComposition = 1 - rHybridComposition
    
    pacbioMatch = { x: int(df["pacbio_match_" + x].sum()) for x in [Q,R,H,HG] }
    pacbioMismatch = { x: int(df["pacbio_mismatch_" + x].sum()) for x in [Q,R,H,HG] }
    pacbioIndel = { x: int(df["pacbio_indel_" + x].sum()) for x in [Q,R,H,HG] }
    pacbioReadCount = { x: int(df["pacbio_read_count_" + x].sum()) for x in [Q,R,H,HG] }

    def pacbio_pie(x):
        labs = ["match", "mismatch", "indel"]
        values = [pacbioMatch[x], pacbioMismatch[x], pacbioIndel[x]]

        plt.pie(values, labels=labs, autopct='%1.1f%%')
        plt.title(x)
        plt.axis("equal")
        plt.show()
    
    for x in [Q,R,H,HG]: pacbio_pie(x)