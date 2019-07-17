class Parameters:
    def __init__(self):
        
        prefix1="/Users/allen bao/Documents/hybrid-pipeline/blocks"
        prefix2="/Users/allen bao/Documents/assembly_data"
        prefix3="/Users/allen bao/Documents/assembly_data/summary_files"
        
        
        self.summaryFile = prefix3 + "/CF062.summary.txt"
        self.queryFasta = prefix2 + "/CF062/OSK7121_03C_supernova.pseudohap.fasta"
        self.referenceFasta = prefix2 + "/CF062/CF062B2D.contigs.fasta.PILON2.fasta"
        self.arks_output = prefix1 + "/CF062B2D.contigs.fasta.PILON2.fasta.tigpair_checkpoint.tsv"
        
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
        
        self.arks               =   True    #run arks/linked reads step of pipeline for improved scaffolding
        self.LR_THRESHOLD       =   2       #minimum number of acceptable barcodes in common between canu contigs for scaffolding to proceed
        
    def add_report(self, report):
        if self.STATISTICS:
            self.reports.append(report)