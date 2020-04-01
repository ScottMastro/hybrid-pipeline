import argparse

class Parameters:
    def __init__(self):
        
        self.REF_REGION = None

        
        self.REF_FA = None
        self.REF_ALIGNED_READS = None

        self.CURATED_ALIGN = None
        self.CURATED_QUERY_FA = None
        self.CURATED_REF_FA = None

        self.ALIGN = None
        self.QUERY_FA = None
        self.QUERY_ALIGNED_READS = None
        self.CANU_DIR = None
        
        self.ARKS = None
        self.UNITIGS = None

        self.OUTPUT_DIR         =   "."     #output directory
        
        self.PLOT_BLOCKS        =   None
        self.DOT_PLOT           =   None

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
        self.MBLOCK_OVERLAP     =   100     #max overlap for megablocks in bp
        self.MIN_MBLOCK_SIZE    =   5000    #min size for megablocks in bp
        self.MAX_MBLOCK_DIST    =   1e7     #max distance between two megablocks

        self.VERBOSE            =   1       #verbosity
        self.WAIT               =   True    #wait for input at high verbosity 
        self.STATISTICS         =   True    #generate and output statistics      
        
        self.BLOCK_DIST_THRESH  =   1.5     #threshold for block insert dist

        self.arks               =   False    #run arks/linked reads step of pipeline for improved scaffolding
        self.LR_THRESHOLD       =   2      #minimum number of acceptable barcodes in common between canu contigs for scaffolding to proceed
        
        self.TRUST_REF_CHUNKS   =   True    #trust the reference at the chunk level

        
    def add_report(self, report):
        if self.STATISTICS:
            self.reports.append(report)


def get_parameters_reference_polish(CFID="CF002"):
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()
    '''

    p.REF_FA =  "/media/scott/HDD/sickkids/slc9a3/hg38.fa"
    p.REF_ALIGNED_READS = "/media/scott/HDD/sickkids/slc9a3/CF001.slc9a3.bam"
    p.OUTPUT_DIR = "/media/scott/HDD/sickkids/slc9a3/out"
    
    '''

    p.REF_FA =  "/media/scott/Rotom/hybrid2/slc9a3/hg38.fa"
    p.REF_ALIGNED_READS = "/media/scott/Rotom/hybrid2/slc9a3/" + CFID + ".slc9a3.bam"
    p.OUTPUT_DIR = "/media/scott/Rotom/hybrid2/slc9a3/out/" + CFID

    p.REF_REGION = "chr5:393462-696129"

    parser = argparse.ArgumentParser(description="Reference sequence region polish")

    # Positional
    parser.add_argument("refFa", metavar="reference", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    parser.add_argument("refBam", metavar="reads", default=p.REF_ALIGNED_READS, nargs="?",
                    help="Read BAM aligned to reference")
   
    parser.add_argument("refRegion", metavar="region", default=p.REF_REGION, nargs="?",
                help="Region to polish, formatted as chrom:start-end")

    # Optional
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                    help="Directory where output will be written." )
    
    # Parse
    args = parser.parse_args()

    p.REF_FA = args.refFa
    p.REF_REGION = args.refRegion
    p.REF_ALIGNED_READS = args.refBam

    p.OUTPUT_DIR = args.outdir
    
    return p
   
           
def get_parameters_polish():
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()

    p.ALIGN = "/media/scott/HDD/sickkids/CF062/CF062_supernova_to_canu.paf"
    p.QUERY_FA = "/media/scott/HDD/sickkids/CF062/OSK7121_03C_supernova.megabubbles.fasta.gz"
    p.REF_FA =  "/media/scott/HDD/sickkids/CF062/CF062B2D.contigs.fasta.PILON2.fasta"
    p.OUTPUT_DIR = "/media/scott/HDD/sickkids/CF062/out"
    p.REF_ALIGNED_READS = "/media/scott/HDD/sickkids/CF062/pacbio_to_canu.pbmm2.reheader.bam"
    p.CANU_DIR = "/media/scott/HDD/sickkids/CF062/canu_data"

    #defaultUNITIGS = "/media/scott/HDD/sickkids/CF062/HG002_NA24385_son_57X.unitigs.bed"
    
    '''
    p.CURATED_ALIGN = "/media/scott/Rotom/hybrid2/curated/supernova_to_canu.curated.paf"
    p.CURATED_QUERY_FA = "/media/scott/Rotom/hybrid2/curated/OSK7121_supernova.psuedohap.10k.fasta.gz"
    p.CURATED_REF_FA = "/media/scott/Rotom/hybrid2/curated/CF062.canu.curated.fasta"

    p.ALIGN = "/media/scott/Rotom/hybrid2/supernova_to_canu.paf"
    p.UNITIG_ALIGN = "/media/scott/Rotom/hybrid2/unitigs_to_canu.paf"

    p.QUERY_FA = "/media/scott/Rotom/hybrid2/OSK7121_03C_supernova.megabubbles.fasta.gz"
    p.REF_FA = "/media/scott/Rotom/hybrid2/CF062.canu.fasta"
    p.REF_ALIGNED_READS = "/media/scott/Rotom/hybrid2/pacbio_to_canu.pbmm2.reheader.bam"
    p.OUTPUT_DIR = "/media/scott/Rotom/hybrid2/out"
    p.QUERY_ALIGNED_READS = "/media/scott/Rotom/hybrid2/pacbio_to_supernova.bam"
    p.CANU_DIR = "/media/scott/Rotom/hybrid2/canu_data"

    '''



    parser = argparse.ArgumentParser(description="Hybrid assembly tool")


    #Positional
    parser.add_argument("alignments", metavar="ALN", default=p.ALIGN, nargs="?", 
                        help="Tab-separated output of BLAST aligner" )
    parser.add_argument("qfasta", metavar="Q_FA", default=p.QUERY_FA, nargs="?",
                        help="Query FASTA file")
    parser.add_argument("rfasta", metavar="R_FA", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    #Optional
    parser.add_argument("--unitigs", type=str, default=p.UNITIGS, 
                    help="Unitig BED file")

    parser.add_argument("--chunk_size", type=int, default=p.CHUNK_SIZE, 
                        help="Size of each chunk in base pairs")
    parser.add_argument("--chunk_per_block", type=int, default=p.CHUNK_PER_BLOCK, 
                        help="Minimum number of chunks per block")
    parser.add_argument("--chunk_max_dist", type=int, default=p.CHUNK_MAX_DIST,
                        help="Maximum distance between chunks")
    parser.add_argument("--chunk_min_ident", type=float, default=p.CHUNK_MIN_IDENT,
                        help="Minimum percent identity for chunk")
    parser.add_argument("--chunk_skip", type=int, default=p.CHUNK_SKIP,
                        help="Maximum # of chunk skip allowed")
    parser.add_argument("--chunk_skip_bp", type=int, default=p.CHUNK_SKIP_BP, 
                        help="Maximum skip allowed in bp")
    parser.add_argument("--chunk_overlap", type=float, default=p.CHUNK_OVERLAP,
                        help="Max overlap for chunk")
    parser.add_argument("--mblock_rdist", type=int, default=p.MBLOCK_RDIST, 
                        help="Max ref distance between blocks")
    parser.add_argument("--mblock_qdist", type=int, default=p.MBLOCK_QDIST, 
                        help="Max query distance between blocks")
    parser.add_argument("--mblock_qlen_factor", type=float, default=p.MBLOCK_QLEN_FACTOR, 
                        help="Tolerate distance based on qlen")
    parser.add_argument("--align_buffer", type=int, default=p.ALIGN_BUFFER,
                        help="Bp buffer used when aligning")
    parser.add_argument("--block_dist_thresh", type=float, default=p.BLOCK_DIST_THRESH,
                        help="Threshold for block insert dist")
    parser.add_argument("--mblock_overlap", type=int, default=p.MBLOCK_OVERLAP,
                        help="Max overlap for megablocks in bp")
    parser.add_argument("--min_mblock_size", type=int, default=p.MIN_MBLOCK_SIZE,
                    help="Min size for megablocks in bp")
    parser.add_argument("--max_mblock_dist", type=int, default=p.MAX_MBLOCK_DIST,
                help="Max distance between two megablocks in bp")
 
    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                        help="Verbosity")
    parser.add_argument("--wait", type=bool, default=p.WAIT,
                    help="Wait for user input at high verbosity levels")

    parser.add_argument("--plot_blocks", default=p.PLOT_BLOCKS, nargs="?",
                        help="Directory to save plots of the blocks formed")
    parser.add_argument("--dot_plot", default=p.DOT_PLOT, nargs="?",
                        help="Directory to save dot plots of gaps filled")

    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                    help="Directory where output will be written." )
    
    parser.add_argument("--refBam", type=str, default=p.REF_ALIGNED_READS,
                    help="Read BAM aligned to reference." )

    args = parser.parse_args()

    p.ALIGN = args.alignments
    p.QUERY_FA = args.qfasta
    p.REF_FA = args.rfasta
    p.UNITIGS = args.unitigs
    p.REF_REGION = args.region
    p.REF_ALIGNED_READS = args.refBam

    p.CHUNK_SIZE = args.chunk_size
    p.CHUNK_PER_BLOCK = args.chunk_per_block
    p.CHUNK_MAX_DIST = args.chunk_max_dist
    p.CHUNK_MIN_IDENT = args.chunk_min_ident
    p.CHUNK_SKIP = args.chunk_skip
    p.CHUNK_SKIP_BP = args.chunk_skip_bp
    p.CHUNK_OVERLAP = args.chunk_overlap
    p.MBLOCK_RDIST = args.mblock_rdist
    p.MBLOCK_QDIST = args.mblock_qdist
    p.MBLOCK_QLEN_FACTOR = args.mblock_qlen_factor
    p.ALIGN_BUFFER = args.align_buffer
    p.BLOCK_DIST_THRESH = args.block_dist_thresh
    p.MBLOCK_OVERLAP = args.mblock_overlap
    p.MIN_MBLOCK_SIZE = args.min_mblock_size
    p.MAX_MBLOCK_DIST = args.max_mblock_dist
    
    p.PLOT_BLOCKS = args.plot_blocks
    p.DOT_PLOT = args.dot_plot

    p.VERBOSE = args.verbose
    p.WAIT = args.wait
    
    p.OUTPUT_DIR = args.outdir
    
    return p
   

def get_parameters():
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()
    
    p.ALIGN = "/media/scott/HDD/sickkids/CF062/CF062_supernova_to_canu.paf"
    p.QUERY_FA = "/media/scott/HDD/sickkids/CF062/OSK7121_03C_supernova.megabubbles.fasta.gz"
    p.REF_FA =  "/media/scott/HDD/sickkids/CF062/CF062B2D.contigs.fasta.PILON2.fasta"
    p.OUTPUT_DIR = "/media/scott/HDD/sickkids/CF062/out"
    p.CANU_DIR = "/media/scott/HDD/sickkids/CF062/canu_data"

    #defaultUNITIGS = "/media/scott/HDD/sickkids/CF062/HG002_NA24385_son_57X.unitigs.bed"
    
    '''

    p.ALIGN = "/media/scott/Rotom/hybrid2/supernova_to_canu.paf"
    p.UNITIG_ALIGN = "/media/scott/Rotom/hybrid2/unitigs_to_canu.paf"

    p.QUERY_FA = "/media/scott/Rotom/hybrid2/OSK7121_03C_supernova.megabubbles.fasta.gz"
    p.REF_FA = "/media/scott/Rotom/hybrid2/CF062.canu.fasta"
    p.REF_ALIGNED_READS = "/media/scott/Rotom/hybrid2/pacbio_to_canu.pbmm2.reheader.bam"
    p.OUTPUT_DIR = "/media/scott/Rotom/hybrid2/out"
    p.QUERY_ALIGNED_READS = "/media/scott/Rotom/hybrid2/pacbio_to_supernova.bam"
    p.CANU_DIR = "/media/scott/Rotom/hybrid2/canu_data"

    '''

    parser = argparse.ArgumentParser(description="Hybrid assembly tool")


    #Positional
    parser.add_argument("alignments", metavar="ALN", default=p.ALIGN, nargs="?", 
                        help="Tab-separated output of BLAST aligner" )
    parser.add_argument("qfasta", metavar="Q_FA", default=p.QUERY_FA, nargs="?",
                        help="Query FASTA file")
    parser.add_argument("rfasta", metavar="R_FA", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    #Optional
    parser.add_argument("--unitigs", type=str, default=p.UNITIGS, 
                    help="Unitig BED file")

    parser.add_argument("--chunk_size", type=int, default=p.CHUNK_SIZE, 
                        help="Size of each chunk in base pairs")
    parser.add_argument("--chunk_per_block", type=int, default=p.CHUNK_PER_BLOCK, 
                        help="Minimum number of chunks per block")
    parser.add_argument("--chunk_max_dist", type=int, default=p.CHUNK_MAX_DIST,
                        help="Maximum distance between chunks")
    parser.add_argument("--chunk_min_ident", type=float, default=p.CHUNK_MIN_IDENT,
                        help="Minimum percent identity for chunk")
    parser.add_argument("--chunk_skip", type=int, default=p.CHUNK_SKIP,
                        help="Maximum # of chunk skip allowed")
    parser.add_argument("--chunk_skip_bp", type=int, default=p.CHUNK_SKIP_BP, 
                        help="Maximum skip allowed in bp")
    parser.add_argument("--chunk_overlap", type=float, default=p.CHUNK_OVERLAP,
                        help="Max overlap for chunk")
    parser.add_argument("--mblock_rdist", type=int, default=p.MBLOCK_RDIST, 
                        help="Max ref distance between blocks")
    parser.add_argument("--mblock_qdist", type=int, default=p.MBLOCK_QDIST, 
                        help="Max query distance between blocks")
    parser.add_argument("--mblock_qlen_factor", type=float, default=p.MBLOCK_QLEN_FACTOR, 
                        help="Tolerate distance based on qlen")
    parser.add_argument("--align_buffer", type=int, default=p.ALIGN_BUFFER,
                        help="Bp buffer used when aligning")
    parser.add_argument("--block_dist_thresh", type=float, default=p.BLOCK_DIST_THRESH,
                        help="Threshold for block insert dist")
    parser.add_argument("--mblock_overlap", type=int, default=p.MBLOCK_OVERLAP,
                        help="Max overlap for megablocks in bp")
    parser.add_argument("--min_mblock_size", type=int, default=p.MIN_MBLOCK_SIZE,
                    help="Min size for megablocks in bp")
    parser.add_argument("--max_mblock_dist", type=int, default=p.MAX_MBLOCK_DIST,
                help="Max distance between two megablocks in bp")
 
    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                        help="Verbosity")
    parser.add_argument("--wait", type=bool, default=p.WAIT,
                    help="Wait for user input at high verbosity levels")

    parser.add_argument("--plot_blocks", default=p.PLOT_BLOCKS, nargs="?",
                        help="Directory to save plots of the blocks formed")
    parser.add_argument("--dot_plot", default=p.DOT_PLOT, nargs="?",
                        help="Directory to save dot plots of gaps filled")

    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                    help="Directory where output will be written." )
    
    args = parser.parse_args()

    p.ALIGN = args.alignments
    p.QUERY_FA = args.qfasta
    p.REF_FA = args.rfasta
    p.UNITIGS = args.unitigs

    p.CHUNK_SIZE = args.chunk_size
    p.CHUNK_PER_BLOCK = args.chunk_per_block
    p.CHUNK_MAX_DIST = args.chunk_max_dist
    p.CHUNK_MIN_IDENT = args.chunk_min_ident
    p.CHUNK_SKIP = args.chunk_skip
    p.CHUNK_SKIP_BP = args.chunk_skip_bp
    p.CHUNK_OVERLAP = args.chunk_overlap
    p.MBLOCK_RDIST = args.mblock_rdist
    p.MBLOCK_QDIST = args.mblock_qdist
    p.MBLOCK_QLEN_FACTOR = args.mblock_qlen_factor
    p.ALIGN_BUFFER = args.align_buffer
    p.BLOCK_DIST_THRESH = args.block_dist_thresh
    p.MBLOCK_OVERLAP = args.mblock_overlap
    p.MIN_MBLOCK_SIZE = args.min_mblock_size
    p.MAX_MBLOCK_DIST = args.max_mblock_dist
    
    p.PLOT_BLOCKS = args.plot_blocks
    p.DOT_PLOT = args.dot_plot

    p.VERBOSE = args.verbose
    p.WAIT = args.wait
    
    p.OUTPUT_DIR = args.outdir
    
    return p
   
    