class Parameters:
    def __init__(self):
        
        prefix1="/home/scott/Dropbox/hybrid-pipeline/blocks"
        #prefix2="/media/scott/HDD/sickkids"
        prefix2="/media/scott/Rotom/assembly_data"
        
        '''
        self.summaryFile = prefix1 + "/summary2.txt"
        self.queryFasta = prefix2 + "/CF062/OSK7121_03C_supernova.pseudohap2.1.fasta"
        self.referenceFasta = prefix2 + "/CF062/CF062B2D.contigs.fasta.PILON2.fasta"
        '''
        self.summaryFile = prefix1 + "/summary_giab.txt"
        self.queryFasta = prefix2 + "/NA24385/NA24385_supernova.pseudohap2.1.fasta"
        self.referenceFasta = prefix2 + "/NA24385/HG002_NA24385_son_57X.contigs.fasta.PILON2"
        
        self.CHUNK_SIZE         =   1000    #size of each chunk in base pairs
        self.CHUNK_PER_BLOCK    =   3       #minimum number of chunks per block
        self.CHUNK_MAX_DIST     =   5000    #maximum distance between chunks
        self.CHUNK_MIN_IDENT    =   97      #minimum percent identity for chunk
        self.CHUNK_SKIP         =   2       #maximum # of chunk skip allowed
        self.CHUNK_SKIP_BP      =   2000    #maximum skip allowed in bp
        self.CHUNK_OVERLAP      =   0.25    #max overlap for chunk
        
        self.MBLOCK_RDIST       =   20000   #max ref distance between blocks
        self.MBLOCK_QDIST       =   200000  #max query distance between blocks
        self.MBLOCK_QLEN_FACTOR =   0.05    #tolerate distance based on qlen

        self.ALIGN_BUFFER       =   200     #bp buffer used for aligning

        self.BLOCK_DIST_THRESH  =   1.5     #threshold for block insert dist
        self.MEGABLOCK_OVERLAP  =   0.2     #max overlap for megablocks
        self.VERBOSE            =   1       #verbosity
        self.STATISTICS         =   True    #generate and output statistics

        self.reports = []
        
    def add_report(self, report):
        if self.STATISTICS:
            self.reports.append(report)