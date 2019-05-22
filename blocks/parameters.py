class Parameters:
    def __init__(self):
        self.CHUNK_SIZE         =   1000    #size of each chunk in base pairs
        self.CHUNK_PER_BLOCK    =   3       #minimum number of chunks per block
        self.MAX_DIST_CHUNK     =   5000    #maximum distance between chunks
        self.CHUNK_SKIP         =   2       #maximum # of chunk skip allowed
        self.CHUNK_SKIP_BP      =   2000    #maximum skip allowed in bp
        self.CHUNK_OVERLAP      =   0.25    #max overlap for a repeated chunk
        self.BLOCK_DIST_THRESH  =   1.5     #threshold for block insert dist
        self.MEGABLOCK_OVERLAP  =   0.01    #max overlap for megablocks
        self.VERBOSE            =   1       #verbosity
        self.STATISTICS         =   True    #generate and output statistics
        
        self.reports = []
        
    def add_report(self, report):
        if self.STATISTICS:
            self.reports.append(report)