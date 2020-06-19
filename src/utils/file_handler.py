import shutil
import os
import glob
import pickle as pkl

import pandas as pd
import csv
from pybedtools import BedTool

#==================================================
# Basic convenience functions for file manipulation
#==================================================

def copy_file(file, newFile):
    shutil.copyfile(file, newFile)
    return newFile

def move_file(file, newFile):
    shutil.move(file, newFile)
    return newFile

def rename_file(file, newFile):
    os.rename(file, newFile)
    return newFile

def make_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

def file_exists(file):
    return os.path.isfile(file)

def dir_exists(directory):
    return os.path.exists(directory)

def delete_file(file, deleteDirTree=False):
    if os.path.isfile(file):
        os.remove(file)
    if os.path.isfile(file + ".bai"):
        os.remove(file + ".bai")
    if os.path.isfile(file + ".fai"):
        os.remove(file + ".fai")
    if os.path.isfile(file + ".tbi"):
        os.remove(file + ".tbi")

    if os.path.isdir(file):
        try:
            os.rmdir(file)
        except:
            if deleteDirTree: shutil.rmtree(file)

def reset_directory(directory):
    try:
        shutil.rmtree(directory)
    except:
        pass
    
    try:
        os.mkdir(directory)
    except:
        pass
    
def cat(directory, extension, outfile, others=[]):
    with open(outfile, 'wb') as out:
        for other in others:
            with open(other, 'rb') as readfile:
                shutil.copyfileobj(readfile, out)
                
        for filename in glob.glob(directory + "*" + extension):
            if filename == outfile:
                # don't want to copy the output into the output
                continue
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, out)
                
    return outfile

#==================================================
# Pickling functions
#==================================================

def pickle(obj, path, overwrite=False):
    exists = os.path.isfile(path)
    if overwrite or not exists:
         with open(path, 'wb') as handle:
            pkl.dump(obj, handle)
    
def unpickle(path):
    exists = os.path.isfile(path)
    if exists:
        with open(path, 'rb') as handle:
            return pkl.load(handle)
    return None

#==================================================
# Hybrid-specific parsers
#==================================================

csv.field_size_limit(999999999) #needed to properly parse this file
def parse_alignments(csv_file):
    """Loads file containing BLAST alignments."""

    colnames = ["qid", "qsize", "rid", "rlen",
                "chunks",  "rstart", "rend", 
                "qstart", "qend", "firstchunk", #todo: remove firstchunk col
                "nchunks", "alen", "pcid"]
    
    aligndf = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    alignDict = {str(tigId):rows for (tigId,rows) in tuple(aligndf.groupby('qid'))}

    return alignDict

def parse_bed(bedFile):
    """Loads BED file containing unitigs."""
    
    if bedFile is None:
        return None
    
    bed = BedTool(bedFile)        
    return bed