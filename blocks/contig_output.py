from reversecomp import reverse_complement

def validate_forks(path, seqData, nbases=25):

    for fork in path:
        
        print( "=========================" )
        
        source = " NOVA" if fork.switch == 'r' else " CANU"

        id = fork.before_id()
        pos = fork.before_pos()
        seq = seqData[str(id)]
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
        seq = seqData[str(id)]
        if fork.after_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.after_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq + source)       
        print(nbases*"-" + forkSeq[nbases:] + source)

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

def get_edges(sequence, nbases=1000, toPrint=False):
    edges = []
    for i in range(len(sequence)-1):
        edge = sequence[i][-int(nbases/2):] +\
            sequence[i+1][:int(nbases/2)]
        edges.append(edge)
        if toPrint:
            print(">" + str(i))
            print(edge)
        
    return edges


def path_to_sequence(path, seqData, invert=False):
    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        tigId = startFork.after_id()
        s = startFork.after_pos()
        e = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        

        if invert:
            tigId = startFork.before_id()
            s = startFork.before_pos()
            e = endFork.after_pos()
            strand = startFork.before_strand()
            src = startFork.before_switch()  
            
            if tigId is None or e is None:
                sequence.append("-")
                source.append("-")
                return
            
            if(startFork.before_strand() != endFork.after_strand()):
                e = len(seqData[str(endFork.after_id())]) - endFork.after_pos()
                
            if strand == 1 and s > e: (e,s) = (s,e)
            if strand == -1 and e < s: (e,s) = (s,e)

        seq = seqData[str(tigId)]
        
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
