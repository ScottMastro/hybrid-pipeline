from Bio import SeqIO
import os
import re
import gzip
import pickle as pkl
import pandas as pd
from pybedtools import BedTool

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

def read_fasta(fasta):
    """Loads fasta (or gzipped fasta) file into memory with SeqIO.
    Returns dictionary of fasta ID to sequence characters.
    """

    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rt") as handle:
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
    
def parse_confident_regions(bedFile):
    """Loads BED file containing unitigs."""
    
    bed = BedTool(bedFile)
    return bed
    
 
def path_to_sequence(path, seqData):
    """
    Takes in a path, returns a list of sequences and a list of sources that represent
    the hybrid sequence when concatenated together.
    """

    sequence = []
    source = []
    
    rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
    def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])    

    def add_seq(startFork, endFork):
        
        if startFork.is_Nfork() or endFork.is_Nfork():
            sequence.append("N"*32)
            source.append("N"*32)
            return
        
        tigId = startFork.after_id()
        start = startFork.after_pos()
        end = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        
        
        '''
        if(startFork.after_strand() != endFork.before_strand()):
            print(startFork)
            print(endFork)
            print("strand issue")
            input()
        '''
        
        seq = seqData[str(tigId)]
        
        if strand == -1:
            start = len(seq) - start
            end = len(seq) - end
            start, end = end, start
            
        segment = seq[start:end]
        if strand == -1:
            segment = reverse_complement(segment)

        sequence.append(segment)
        source.append(src*len(segment))
    
    startFork, endFork = None, path[0]
    
    for fork in path[1:]:
        startFork = endFork
        endFork = fork
        add_seq(startFork, endFork)
    
    return (sequence, source)

def write_fasta(fid, sequence, filePath):   
    writer = open(filePath, "w+")
    writer.write(">" + fid + "\n" + \
             re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")


def fasta2dict(faPath, toUpper=False):
    fastaDict = dict()
    seqs = SeqIO.parse(open(faPath),'fasta')
    for f in seqs:
        fastaDict[f.id] = str(f.seq)
        if toUpper: fastaDict[f.id] = fastaDict[f.id].upper()
    return fastaDict

def dict2fasta(faPath, fastaDict, toUpper=False):
    
    writer = open(faPath, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")

    writer.close()
    return faPath


def write_hybrid(scaffolds, seqData, param):
    """
    Writes hybrid and source as FASTA.
    Also writes out source as a BED file.
    """
    
    directory = os.path.dirname(param.OUTPUT_DIR)
    if directory == "": directory = "."
    
    hybridFasta = directory + "/hybrid_assembly.fasta"
    sourceFasta = directory + "/hybrid_source.fasta"
    sourceBed   = directory + "/hybrid_source.bed"
    
    sourceMap = {'r' : 'reference', 'q' : 'query', "N" : "NNN"}
    scoreMap = {'r' : '255', 'q' : '600', "N" : "999"}
    
    hf = open(hybridFasta, "w+")
    sf = open(sourceFasta, "w+")
    sb = open(sourceBed, "w+")

    for i, scaffold in enumerate(scaffolds):
        sequence, source = path_to_sequence(scaffold, seqData)
        tigId = "hybrid" + "_" + str(i)

        hf.write(">" + tigId + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")

        sf.write(">" + tigId + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")

        pos = 0
        for component in source:
            if len(component) > 0:
                x = component[0]
                sb.write("\t".join([tigId, str(pos), str(len(component)-1 + pos), \
                                    sourceMap[x], scoreMap[x], "."]) + "\n")
            pos = pos + len(component)
               
    hf.close()
    sf.close()
    sb.close()