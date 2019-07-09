from Bio import SeqIO
import os
import gzip
import pickle as pkl
import pandas as pd

import csv
csv.field_size_limit(999999999)

def pickle(obj, path, overwrite=False):
    """Pickles an object.
    Optionally overwrite if file already exists at path."""

    exists = os.path.isfile(path)
    if overwrite or not exists:
         with open(path, 'wb') as handle:
            pkl.dump(obj, handle)
    
def unpickle(path):
    """Unpickles an object and return contents.
    Returns None if path does not exist"""

    exists = os.path.isfile(path)
    if exists:
        with open(path, 'rb') as handle:
            return pkl.load(handle)
    return None

def read_fasta(fasta):
    """Loads fasta (or gzipped fasta) file into memory with SeqIO.
    Returns dictionary of fasta ID to sequence characters.
    """

    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rb") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fasta, "fasta"))
        
    return dict(zip([str(r.id) for r in records], [str(r.seq) for r in records]))


def parse_alignments(csv_file):
    """Loads file containing BLAST alignments."""

    colnames = ["qid", "qsize", "rid", "rlen",
                "chunks",  "rstart", "rend", 
                "qstart", "qend", "firstchunk", #todo: remove firstchunk col
                "nchunks", "alen", "pcid"]
    
    df = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    return df