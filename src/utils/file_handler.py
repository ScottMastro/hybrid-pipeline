from pybedtools import BedTool

import os
import re
import pickle as pkl
import pandas as pd

import csv
csv.field_size_limit(999999999)

def parse_alignments(csv_file):
    """Loads file containing BLAST alignments."""

    colnames = ["qid", "qsize", "rid", "rlen",
                "chunks",  "rstart", "rend", 
                "qstart", "qend", "firstchunk", #todo: remove firstchunk col
                "nchunks", "alen", "pcid"]
    
    aligndf = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    alignDict = {str(tigId):rows for (tigId,rows) in tuple(aligndf.groupby('qid'))}

    return alignDict

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

def validate_ids(rids, qids):
    """Verifies that IDs between query and reference FASTA are unique.
    Returns True if IDs are ok, False if duplicate is found"""

    intersect = set(rids).intersection(qids) 
    if len(intersect) > 0:
        print("ERROR: IDs for input fasta files overlap.")
        print("\"" + str(intersect.pop()) + "\" is found in both fasta files.")
        print("Please ensure each ID is unique.")
        return False
    return True

    
def parse_confident_regions(bedFile, rids, param):
    """Loads BED file containing unitigs."""
    
    bed = BedTool(bedFile)
    
    #convert chrom name in bedFile from "ctg" to "tig" (Canu)
    isTig = True
    for rid in rids:
        if not rid.startswith("tig"):
            isTig=False
            break
    isCtg = True
    for inv in bed:
        if not inv.chrom.startswith("ctg"):
            isCtg=False
            break

    if isTig and isCtg:
        reader = open(bedFile, "r")
        name = os.path.basename(bedFile)
        renamed = param.OUTPUT_DIR + "/" + re.sub(".bed", ".renamed.bed", name)
        writer = open(renamed, "w+")
        for line in reader: 
            writer.write(re.sub("ctg", "tig", line))
        reader.close()
        writer.close()
        bed = BedTool(renamed)
        
    return bed
    
 