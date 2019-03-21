reverseMap = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N', '-':'-'}

def reverse_complement(seq):
    return ''.join([reverseMap[x] for x in seq[::-1]])
