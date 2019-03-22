from path import Path
from path import clean_path

from reversecomp import reverse_complement

import aligner
import re

def validate_forks(path, seqData, nbases=25):

    for fork in path:
        
        print( "=========================" )
        
        source = " NOVA" if fork.switch == 'r' else " CANU"

        id = fork.before_id()
        pos = fork.before_pos()
        seq = str(seqData[str(id)])
        if fork.before_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.before_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq[:nbases] + nbases*"-" + source)
        print(forkSeq + source)
        
        #---------------------
        
        source = " CANU" if fork.switch == 'r' else " NOVA"

        id = fork.after_id()
        pos = fork.after_pos()
        seq = str(seqData[str(id)])
        if fork.after_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.after_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq + source)       
        print(nbases*"-" + forkSeq[nbases:] + source)


def find_Ns(sequence, minIdx=-1, maxIdx=-1):
    """
    Find the (start,end) of Ns in sequence
    """
    idx = minIdx
    startIdx = -1
    lastIdx = -1

    while True:
        idx = sequence.find('N', idx+1)
        
        if maxIdx > 0 and idx > maxIdx:
            if startIdx < lastIdx:
                yield (startIdx, lastIdx)
            break

        if(idx != lastIdx + 1):
            if startIdx == -1:
                startIdx = idx
            else:           
                yield (startIdx, lastIdx)
                startIdx = idx
            
        if idx == -1:
            break
        
        lastIdx = idx
'''
def sequence_validation(sequences):
    f = open("validate.fasta", "w+")
    i = 0
    for sequence in sequences:
        
        if len(sequence) < 3001:
            f.write(">segment" + "_" + str(i) + "\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence, 0, re.DOTALL) + "\n")
        else:
            f.write(">segment" + "_" + str(i) + "start\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[:1000], 0, re.DOTALL) + "\n")
            f.write(">segment" + "_" + str(i) + "end\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[-1000:], 0, re.DOTALL) + "\n")

        i = i + 1
    f.close()      
'''

def get_edges(sequence, nbases=2000, toPrint=False):
    edges = []
    for i in range(len(sequence)-1):
        edge = sequence[i][-int(nbases/2):] +\
            sequence[i+1][:int(nbases/2)]
        edges.append(edge)
        if toPrint:
            print(">" + str(i))
            print(edge)
        
    return edges


def fasten_megablock(megablock, seqData, param):
    rid = megablock.rid
    qid = megablock.qid
    qSeq = str(seqData[str(qid)])
    rSeq = str(seqData[str(rid)])
    path = Path()
    
    def fill_Ns(block, start, end):

        alignBuffer=1000
        minPos = min(start, end)
        maxPos = max(start, end)
        pathN = Path()

        for Nstart, Nend in find_Ns(qSeq, minPos, maxPos):
            qstart = int(Nstart - alignBuffer/2)
            qend = int(Nend + alignBuffer/2)

            rstart = block.closest_corresponding_position(qstart, q=True)
            rend = block.closest_corresponding_position(qend, q=True)

            (startFork, endFork) = aligner.create_bubble(rid, qid, rSeq, qSeq, \
                rstart, rend, qstart, qend, alignBuffer)
            
            pathN.add_fork(startFork)
            pathN.add_fork(endFork)

        return pathN

    #prevBlock = megablock[0]
    #block = megablock[1]
    prevBlock = None

    for block in megablock:
        if prevBlock is not None:
            rstart = prevBlock.end(q=False)
            rend = block.start(q=False)
            qstart = prevBlock.end(q=True)
            qend = block.start(q=True)
            
            #aligns the boundary points between blocks
            #creates a ref/query bubble for the space between blocks
            (startFork, endFork) = aligner.create_bubble(rid, qid, \
                rSeq, qSeq, rstart, rend, qstart, qend, alignBuffer=2000)
            
            # find and fill Ns
            Nstart = 0 if path.last_fork() is None else path.last_fork().get_pos(q=True)
            Nend = startFork.get_pos(q=True)

            pathN = fill_Ns(prevBlock, Nstart, Nend)
            path.add_path(pathN)
            path.add_fork(startFork)
            path.add_fork(endFork)

        prevBlock = block

    return path


def path_to_sequence(path, seqData):
    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        tigId = startFork.after_id()
        s = startFork.after_pos()
        e = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        

        seq = str(seqData[str(tigId)])
                
        if strand == 1:
            start, end = s, e
        elif strand == -1:
            start = None if e is None else (len(seq) - e)
            end = None if s is None else (len(seq) - s)
            
        segment = seq[start:end]
        if strand == -1:
            segment = reverse_complement(segment)
    
        sequence.append(segment)
        source.append(src*len(segment))
    
    startFork = None
    endFork = path.head()       
    
    for fork in path:
        startFork = endFork
        endFork = fork
        add_seq(startFork, endFork)

    add_seq(endFork, path.tail())
    
    return (sequence, source)


def fasten_contig(contig, seqData, lengthData, param):

    f = open("hybrid6.fasta", "w+")
    g = open("source6.fasta", "w+")

    #megablock = contig.mblocks[0]
    for megablock in contig.mblocks:
        path = fasten_megablock(megablock, seqData, param)
        path = clean_path(path, lengthData, verbose=True)
        (sequence, source) = path_to_sequence(path, seqData)
        
        f.write(">" + str(megablock.qid) + "_" + str(megablock.rid) + "\n")
        f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
        g.write(">" + str(megablock.qid) + "_" + str(megablock.rid) + "\n")
        g.write(re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")
    f.close()
    
    validate_forks(path, seqData, 25)
    qforks = path.get_fork_sequence(seqData, 3000, q=True)
    rforks = path.get_fork_sequence(seqData, 3000, q=False)
