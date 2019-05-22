import pandas as pd
import csv
csv.field_size_limit(999999999)

def parse_alignments(csv_file):
    
    colnames = ["qid", "qsize", "rid", "rlen",
                "chunks",  "rstart", "rend", 
                "qstart", "qend", 
                "First part aligning to subject",
                "nchunks", "alen", "pcid"]
    
    df = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    return df


def line_to_df(line):
    contig=str(line.iloc[0,0])
    n_chunks=int(line.iloc[0,1])
    chrom=str(line.iloc[0,2])
    chunks=[int(x) for x in line.iloc[0,4].split(',')]
    starts=[int(x) for x in line.iloc[0,5].split(',')]
    ends=[int(x) for x in line.iloc[0,6].split(',')]
    qends=[int(x) for x in line.iloc[0,7].split(',')]
    qstarts=[int(x) for x in line.iloc[0,8].split(',')]

    row_df=pd.DataFrame(pd.to_numeric(chunks))
    row_df.columns=["chunk"]
    row_df["contig"]=contig
    row_df["chrom"]=chrom
    row_df["nchunks"]=n_chunks
    row_df["start"]=pd.to_numeric(starts)
    row_df["end"]=pd.to_numeric(ends)
    row_df["qstart"]=pd.to_numeric(qstarts)
    row_df["qend"]=pd.to_numeric(qends)
    row_df["left"]=row_df[['start','end']].min(axis=1)
    row_df["right"]=row_df[['start','end']].max(axis=1)
    row_df=row_df.sort_values(by=["left"])
    #row_df["size"]=row_df["right"] - row_df["left"] +1
    #row_df["qsize"]=row_df["qstart"] - row_df["qend"]+1

    return row_df
