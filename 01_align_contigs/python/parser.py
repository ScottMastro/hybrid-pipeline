import pandas as pd
import csv
csv.field_size_limit(999999999999)

def parse_alignments(csv_file):
    
    colnames = ["Query sequence ID", 
                "Total number of chunks",
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
