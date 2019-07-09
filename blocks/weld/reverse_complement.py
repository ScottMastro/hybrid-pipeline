dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq):
    '''
    Returns the reverse complement of a sequence.
    '''
    return ''.join([dRc[x] for x in seq[::-1]])