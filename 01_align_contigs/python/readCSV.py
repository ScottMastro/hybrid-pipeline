import pandas as pd
from pandas import ExcelWriter

def read(csv_file):
    
    colnames = ["Query sequence ID", 
                "Total number of parts",
                "Subject sequence ID",
                "Subject sequence length",
                "Parts aligning to subject sequence", 
                "Start of alignment in subject", 
                "End of alignment in subject", 
                "Start of alignment in query", 
                "End of alignment in query", 
                "First part aligning to subject",
                "Total number of aligning parts",
                "Length of query alignment per part"]
    
    df = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    return df


def line2df(line):
    contig=line.iloc[0,0]
    n_chunks=line.iloc[0,1]
    chrom=line.iloc[0,2]
    chunks=[x for x in line.iloc[0,4].split(',')]
    starts=[x for x in line.iloc[0,5].split(',')]
    ends=[x for x in line.iloc[0,6].split(',')]
    qends=[x for x in line.iloc[0,7].split(',')]
    qstarts=[x for x in line.iloc[0,8].split(',')]

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