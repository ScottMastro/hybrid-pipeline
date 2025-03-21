import shutil
import re
import gzip

from Bio import SeqIO
import pyfaidx

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

#==================================================
# Basic FASTA IO functionality 
#==================================================

def read_fasta(fasta, toUpper=False):
    """
    Loads fasta (or gzipped fasta) file into memory with SeqIO.
    Returns dictionary of id->sequence.
    """

    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fasta, "fasta"))

    def get_str(seq):
        return str(seq).upper() if toUpper else str(seq)
    
    faDict = dict(zip([str(r.id) for r in records], [get_str(r.seq) for r in records]))
    return faDict

def write_fasta(faPath, fastaDict, toUpper=False, index=False, charPerLine=64):
    """
    Takes a directory path and dictionary of id->sequence.
    Writes a FASTA file from dictionary sequence.
    """

    if not faPath.endswith(".fa") and not faPath.endswith(".fasta"):
        faPath = faPath + ".fasta"
        
    def get_seq_string(fid):
        seq = re.sub("(.{" + str(charPerLine) + "})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL)
        if toUpper: return seq.upper()
        return seq

    writer = open(faPath, "w+")    
    for fid in fastaDict:
        writer.write(">" + fid + "\n" +  get_seq_string(fid) + "\n")
    writer.close()

    if index: index_fasta(faPath)
    return faPath

def fasta_fetch(faFile, fid, toUpper=False):
    '''Grabs a single FASTA sequence with ID fid.'''
    fa = pyfaidx.Fasta(faFile)
    seq = str(fa[fid])
    if toUpper : return seq.upper()
    return seq

def faidx_fetch(faFile, region, toUpper=False):
    '''Grabs the FASTA sequence corresponding to the Region object.'''
    fa = pyfaidx.Faidx(faFile)
    seq = str(fa.fetch(region.chrom, region.start, region.end))
    if toUpper : return seq.upper()
    return seq

def index_fasta(faFile):
    '''Index a FASTA file.'''
    faidx = pyfaidx.Faidx(faFile)
    return faidx

#==================================================
# FASTA IO convenience functionality 
#==================================================

def get_fasta_len(faFile, fid=None):
    """
    Gets the sequence length of a FASTA sequence with ID fid.
    Returns the first sequence length of fid is None.
    """
    faDict = read_fasta(faFile)
    
    if fid is None:
        fid = list(faDict.keys())[0]
    
    return len(faDict[fid])

def get_first_id(faFile):
    """Returns a string of the first id in the FASTA file."""
    faDict = read_fasta(faFile)
    fid = list(faDict.keys())[0]
    return fid

def rename_single_fasta(faFile, name, toUpper=False):
    """ Rename the ID of a single (not multisequence!) fasta to name. """
    
    faDict = read_fasta(faFile, toUpper=toUpper)
    fid = list(faDict.keys())[0]
    newFaDict = {name : faDict[fid]}
    
    temp = faFile + "__temporaryfa__.fasta"
    write_fasta(temp, newFaDict)
    shutil.move(temp, faFile)
    
    index_fasta(faFile)

    return faFile

def get_fasta_seq(faFile, region=None, toUpper=False):
    """
    Returns a string of the sequence given by region.
    If region is None, returns the first FASTA sequence.
    """
    faDict = read_fasta(faFile, toUpper=toUpper)

    if region is None:
        fid = list(faDict.keys())[0]
        return faDict[fid]

    return faDict[region.chrom][region.start:region.end]

#==================================================
# Hybrid-specific functionality 
#==================================================

def validate_ids(rids, qids):
    """
    Verifies that IDs between query and reference FASTA are unique.
    Returns True if IDs are ok, False if duplicate is found
    """

    intersect = set(rids).intersection(qids) 
    if len(intersect) > 0:
        print("ERROR: IDs for input fasta files overlap.")
        print("\"" + str(intersect.pop()) + "\" is found in both fasta files.")
        print("Please ensure each ID is unique.")
        return False
    return True

def path_to_sequence(path, seqData, numberOfNsInGap=32):
    """
    Takes in a path, returns a list of sequences and a list of sources that represent
    the hybrid sequence when concatenated together.
    """

    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        if startFork.is_Nfork() or endFork.is_Nfork():
            sequence.append("N"*numberOfNsInGap)
            source.append("N"*numberOfNsInGap)
            return
        
        tigId = startFork.after_id()
        start = startFork.after_pos()
        end = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        
                
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

def write_hybrid(scaffolds, seqData, param):
    """
    Writes hybrid and source as FASTA.
    Also writes out source as a BED file.
    """
        
    hybridFasta = param.OUTPUT_DIR + "/hybrid_assembly.fasta"
    sourceFasta = param.OUTPUT_DIR + "/hybrid_source.fasta"
    sourceBed   = param.OUTPUT_DIR + "/hybrid_source.bed"
    
    sourceMap = {'r' : 'reference', 'q' : 'query', "N" : "NNN"}
    scoreMap = {'r' : '255', 'q' : '600', "N" : "999"}
    
    hf = open(hybridFasta, "w+")
    sf = open(sourceFasta, "w+")
    sb = open(sourceBed, "w+")

    for i, scaffold in enumerate(scaffolds):
        sequence, source = path_to_sequence(scaffold, seqData)
        tigId = scaffold.pid

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
    
def write_leftover(regions, seqData, param, minSize=10000, compressNs=10, spacerNs=500):
        
    fa = param.OUTPUT_DIR + "/hybrid_assembly_leftover.fasta"
    bed   = param.OUTPUT_DIR + "/hybrid_source_leftover.bed"
        
    f = open(fa, "w+")
    b = open(bed, "w+")

    fastaId="leftovers"
    f.write(">" + fastaId + "\n")

    pos=0
    sequences=[]
    for i, region in enumerate(regions):
        seq = seqData[region.chrom][region.start:region.end]
        split = seq.split('N')
        split = [x for x in split if len(x) > 0]
        if sum([len(x) for x in split]) < minSize:
            continue            

        seq = ("N"*compressNs).join(split)
        sequences.append(seq)
        
        b.write("\t".join([fastaId, str(pos), str(len(seq)-1 + pos), str(region)]) + "\n")
        pos += len(seq) + spacerNs
        
    sequence = ("N"*spacerNs).join(sequences)
    f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL))
        
    f.close()
    b.close()

