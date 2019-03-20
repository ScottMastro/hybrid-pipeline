from segment import Segment
from path import Path

import aligner
import re

def segments_to_sequence(segments, refData, queryData):
    sequence = []
    source = []
    for segment in segments:
        
        segId = str(segment.id)
        
        if segId in queryData:
            qSeq = str(queryData[segId])
            sequence.append(segment.get_sequence(qSeq))
            source.append("S"*len(sequence[-1]))
        if segId in refData:
            rSeq = str(refData[segId])
            sequence.append(segment.get_sequence(rSeq))
            source.append("P"*len(sequence[-1]))

        if 'N' in sequence[-1]:
            print(segment)
            
    return (sequence, source)

def validate_edges(segments, refData, queryData):

    prevSegment = None
    for segment in segments:
        if prevSegment == None:
            prevSegment = segment
            continue
        
        segId = str(segment.id)
        prevSegId = str(prevSegment.id)
        nbases=25

        print( "=========================" )
        if prevSegId in queryData:
            print("NOVA")
            qSeq = str(queryData[prevSegId])
            print(prevSegment.get_sequence(qSeq)[-nbases:] + nbases*"-")
            print(prevSegment.get_validation_sequence(qSeq,'r', nbases*2))

        if prevSegId in refData:
            print("CANU")
            rSeq = str(refData[prevSegId])
            print(prevSegment.get_sequence(rSeq)[-nbases:] + nbases*"-")
            print(prevSegment.get_validation_sequence(rSeq,'r', nbases*2))

        if segId in queryData:
            print("NOVA")
            qSeq = str(queryData[segId])
            print(segment.get_validation_sequence(qSeq,'l', nbases*2))
            print(nbases*"-" + segment.get_sequence(qSeq)[:nbases] )

        if segId in refData:
            print("CANU")
            rSeq = str(refData[segId])
            print(segment.get_validation_sequence(rSeq,'l', nbases*2) )
            print(nbases*"-" +segment.get_sequence(rSeq)[:nbases])

        prevSegment = segment
        
        
def get_edges(sequence, nbases=1000):

    edges = []
    prevSeq = None
    for seq in sequence:
        if prevSeq is not None:
            edges.append(prevSeq[-nbases:] + seq[:nbases])
        
        prevSeq = seq

    return edges
        

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

def clean_segments(segments, param):
    cleaned = []
    prevSegment = None
    flipStrand = False
    
    for segment in segments:
        
        if flipStrand:
            segment.make_reverse_complement()
            
        if not segment.valid_strands():
            segment.flip_end()            
            if segment.is_valid():
                flipStrand = not flipStrand
        
        if segment.is_valid():
    
            if prevSegment is not None:
                mergedSegment = prevSegment.try_merge(segment)
                if mergedSegment is None:
                    cleaned.append(prevSegment)
                    prevSegment = segment
                else:
                    segment = mergedSegment
                    
            prevSegment = segment
    
    if prevSegment is not None:
        cleaned.append(prevSegment)
        
    return cleaned

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
    
def fasten_megablock(megablock, refData, queryData, param):
    rid = megablock.rid
    qid = megablock.qid
    qSeq = str(queryData[str(qid)])
    rSeq = str(refData[str(rid)])
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


def fasten_contig(contig, refData, queryData, param):

    f = open("hybrid6.fasta", "w+")
    g = open("source6.fasta", "w+")

    #megablock = contig.mblocks[0]
    for megablock in contig.mblocks:
        segments = fasten_megablock(megablock, refData, queryData, param)
        (sequence, source) = segments_to_sequence(segments, refData, queryData)
        
        f.write(">" + str(megablock.qid) + "_" + str(megablock.rid) + "\n")
        f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
        g.write(">" + str(megablock.qid) + "_" + str(megablock.rid) + "\n")
        g.write(re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")
    f.close()
    
    validate_edges(segments, refData, queryData)
    edges = get_edges(sequence, 500)
