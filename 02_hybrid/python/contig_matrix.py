#import libraries
#----------------
import pandas as pd

#functions
#----------------

def dfToDict(df):
    ''' Convert matrix df to a dictionary consist of two parts: 
    d1={columns:[connected rows]} and d2={rows:[connected columns]}
    the connected rows/columns have non-NaN entries
    '''
    d={};
    for columns in df:
        d[columns]= list(df[df[columns].notna()].index)     
    
    for rows in df.index:
        d[rows]=list(df.loc[rows][df.loc[rows].notna().values].index)
        
    return d

def construct_matrix(novaTigs, canuTigs, blockdf, threshold=35, asDict=False):
    blockdf = blockdf.sort_values(by=["contig", "chrom", "left"])
    blockdf["cov"] = round(100.0*(blockdf["right"]-blockdf["left"])/blockdf["span"], 2)
    blockdf["below_nova"] = blockdf["contig"].shift(-1)
    blockdf["below_canu"] = blockdf["chrom"].shift(-1)
    
    matrix = pd.DataFrame(columns=novaTigs, index=canuTigs)
    covs = blockdf.groupby(['chrom', 'contig'])["cov"].sum()

    for pair in covs.index:
        cov = covs[pair]
        if(cov >= threshold):
            matrix.at[pair[0], pair[1]] = covs[pair]

    if(asDict):
        return dfToDict(matrix)
    
    return matrix

