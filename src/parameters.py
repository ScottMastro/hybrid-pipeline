#==================================================
# Parameter class 
#==================================================

class Parameters:

    #default parameters established here
    def __init__(self):
        
        self.OUTPUT_DIR   = "./out"     #output directory
        self.VERBOSE      = 1           #verbosity
        self.WAIT         = False       #wait for input at high verbosity 

        self.HG38         = None
        self.HG38_INDEX   = None

        #hybrid specific parameters
        self.SUMMARY      = None
        self.QUERY_FA     = None
        self.REF_FA       = None
        self.REF_BED      = None

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
        self.MIN_MBLOCK_SIZE    =   2500    #min size for megablocks in bp
        self.MAX_MBLOCK_DIST    =   1e7     #max distance between two megablocks
        self.BLOCK_DIST_THRESH  =   1.5     #threshold for block insert dist
        self.TRUST_REF_CHUNKS   =   False   #trust the reference at the chunk level

        # polish specific parameters
        self.FASTA                  = None
        self.REF_ALIGNED_READS      = None
        self.QUERY_ALIGNED_READS    = None
        self.TARGET_CONTIG          = None
        self.KEEP_INTERMEDIATE      = False
        self.ANALYSIS               = False

        # graph specific parameters
        self.TSV_INPUT       = None
        self.REF_GENOME      = None
        self.SKIP_EXTRACTION = False
        self.TARGET_REGION   = None
        self.FILTER_X        = 200     #filter out paf alignments less than this size in bp
        self.COMBINE_Y       = 20000   #combine paf alignments if they are closer than this size in bp

#==================================================
# CLI help 
#==================================================

def set_hybrid_parameters(parser, defaultParameters=None):

    p = Parameters() if defaultParameters is None else defaultParameters
    
    #Positional
    parser.add_argument("alignments", metavar="ALN", default=p.SUMMARY, nargs="?", 
                        help="Tab-separated output of BLAST aligner" )
    parser.add_argument("qfasta", metavar="Q_FA", default=p.QUERY_FA, nargs="?",
                        help="Query FASTA file")
    parser.add_argument("rfasta", metavar="R_FA", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    #Optional
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                help="Directory where output will be written." )

    parser.add_argument("--confident", type=str, default=p.REF_BED, 
                    help="BED file of high-confidence reference regions")

    parser.add_argument("--chunk_size", type=int, default=p.CHUNK_SIZE, 
                        help="Size of each chunk in base pairs. Default=" + str(p.CHUNK_SIZE))
    parser.add_argument("--chunk_per_block", type=int, default=p.CHUNK_PER_BLOCK, 
                        help="Minimum number of chunks per block. Default=" + str(p.CHUNK_PER_BLOCK))
    parser.add_argument("--chunk_max_dist", type=int, default=p.CHUNK_MAX_DIST,
                        help="Maximum distance between chunks. Default=" + str(p.CHUNK_MAX_DIST))
    parser.add_argument("--chunk_min_ident", type=float, default=p.CHUNK_MIN_IDENT,
                        help="Minimum percent identity for chunk. Default=" + str(p.CHUNK_MIN_IDENT))
    parser.add_argument("--chunk_skip", type=int, default=p.CHUNK_SKIP,
                        help="Maximum # of chunk skip allowed. Default=" + str(p.CHUNK_SKIP))
    parser.add_argument("--chunk_skip_bp", type=int, default=p.CHUNK_SKIP_BP, 
                        help="Maximum skip allowed in bp. Default=" + str(p.CHUNK_SKIP_BP))
    parser.add_argument("--chunk_overlap", type=float, default=p.CHUNK_OVERLAP,
                        help="Max overlap for chunk. Default=" + str(p.CHUNK_OVERLAP))
    parser.add_argument("--mblock_rdist", type=int, default=p.MBLOCK_RDIST, 
                        help="Max ref distance between blocks. Default=" + str(p.MBLOCK_RDIST))
    parser.add_argument("--mblock_qdist", type=int, default=p.MBLOCK_QDIST, 
                        help="Max query distance between blocks. Default=" + str(p.MBLOCK_QDIST))
    parser.add_argument("--mblock_qlen_factor", type=float, default=p.MBLOCK_QLEN_FACTOR, 
                        help="Tolerate distance based on qlen. Default=" + str(p.MBLOCK_QLEN_FACTOR))
    parser.add_argument("--align_buffer", type=int, default=p.ALIGN_BUFFER,
                        help="Bp buffer used when aligning. Default=" + str(p.ALIGN_BUFFER))
    parser.add_argument("--block_dist_thresh", type=float, default=p.BLOCK_DIST_THRESH,
                        help="Threshold for block insert dist. Default=" + str(p.BLOCK_DIST_THRESH))
    parser.add_argument("--mblock_overlap", type=int, default=p.MBLOCK_OVERLAP,
                        help="Max overlap for megablocks in bp. Default=" + str(p.MBLOCK_OVERLAP))
    parser.add_argument("--min_mblock_size", type=int, default=p.MIN_MBLOCK_SIZE,
                    help="Min size for megablocks in bp. Default=" + str(p.MIN_MBLOCK_SIZE))
    parser.add_argument("--max_mblock_dist", type=int, default=p.MAX_MBLOCK_DIST,
                help="Max distance between two megablocks in bp. Default=" + str(p.MAX_MBLOCK_DIST))
 
    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                        help="Verbosity. Default=" + str(p.VERBOSE))
    parser.add_argument("--wait", type=bool, default=p.WAIT,
                    help="Wait for user input at high verbosity levels. Default=False")

    return parser

def set_polish_parameters(parser, defaultParameters=None):

    p = Parameters() if defaultParameters is None else defaultParameters

    #Positional
    parser.add_argument("fasta", metavar="FA", default=p.FASTA, nargs="?",
                        help="Assembly FASTA file")
    parser.add_argument("qbam", metavar="Q_BAM", default=p.QUERY_ALIGNED_READS, nargs="?",
                    help="Query reads aligned to FA (BAM file)")
    parser.add_argument("rbam", metavar="R_BAM", default=p.REF_ALIGNED_READS, nargs="?",
                    help="Reference reads aligned to FA (BAM file)")
    #Optional
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                help="Directory where output will be written." ) 
    parser.add_argument("-k", "--keep", type=bool, default=p.KEEP_INTERMEDIATE,
            help="Keep all intermediate files. Default=False" )
    parser.add_argument("-t", "--tig", type=str, default=p.TARGET_CONTIG,
                help="Target contig." )
    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                    help="Verbosity. Default=" + str(p.VERBOSE))


    return parser

def set_graph_parameters(parser, defaultParameters=None):
    p = Parameters() if defaultParameters is None else defaultParameters

    parser.add_argument("-i", "--input", type=str, default=p.TSV_INPUT,
            help="Tab-separated list of input FASTA. First column is sample ID, second column is FASTA directory")
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                help="Directory where output will be written." )
    parser.add_argument("-r", "--ref", type=str, default=p.REF_GENOME,
                help="FASTA file with sequence to target." )
    parser.add_argument("-s", "--subset", type=str, default=p.TARGET_REGION,
            help="Extract a subsequence from ref FASTA as target (optional). Format as chr:start-end." )
    parser.add_argument("-x", "--filter", type=int, default=p.FILTER_X,
            help="Filter alignments smaller than this distance in bp." )
    parser.add_argument("-y", "--collapse", type=int, default=p.COMBINE_Y,
            help="Combine alignments that are separated by a distance smaller than this threshold in bp." )
    parser.add_argument("-k", "--skip", type=bool, default=p.SKIP_EXTRACTION,
            help="Skip over FASTA extraction step and build graph from existing PAF files." )
    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                        help="Verbosity. Default=" + str(p.VERBOSE))

    return parser

#==================================================
# CLI parsing 
#==================================================

def get_hybrid_parameters(parser):

    p = Parameters()
    args = parser.parse_args()

    p.QUERY_FA = args.qfasta
    p.REF_FA = args.rfasta
    p.SUMMARY = args.alignments
    p.OUTPUT_DIR = args.outdir
    p.REF_BED = args.confident
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
    p.VERBOSE = args.verbose
    p.WAIT = args.wait
    
    return p

def get_polish_parameters(parser):

    p = Parameters()
    args = parser.parse_args()
    
    #set parameters from user input
    p.FASTA = args.fasta
    p.QUERY_ALIGNED_READS = args.qbam
    p.REF_ALIGNED_READS = args.rbam
    p.TARGET_CONTIG = args.tig
    p.OUTPUT_DIR = args.outdir
    p.VERBOSE = args.verbose
    p.WAIT = args.wait
    p.KEEP_INTERMEDIATE = args.keep

    return p

def get_graph_parameters(parser):
   
    p = Parameters()
    args = parser.parse_args()

    #set parameters from user input
    p.TSV_INPUT = args.input
    p.OUTPUT_DIR = args.outdir
    p.REF_GENOME = args.ref
    p.SKIP_EXTRACTION = args.skip
    p.FILTER_X = args.filter
    p.COMBINE_Y = args.collapse

    try:
        p.TARGET_REGION = args.subset
    except:
        p.TARGET_REGION = None

    return p

#==================================================
# Parameter validation 
#==================================================

def validate_hybrid_parameters(p):

    if p.QUERY_FA is None: return False
    if p.REF_FA   is None: return False
    if p.SUMMARY  is None: return False
    
    return True

def validate_polish_parameters(p):
    
    if p.FASTA                is None: return False
    if p.QUERY_ALIGNED_READS  is None: return False
    if p.REF_ALIGNED_READS    is None: return False
    
    return True


def validate_graph_parameters(p):
    
    return True


#==================================================
# Explicitly specify parameters (for debugging, running from IDE) 
#==================================================

def explicit_hybrid_params():
    
    p = Parameters()

    #override default parameters
    prefix = "/media/scott/Zapdos/CF062_19/"
    p.SUMMARY = prefix + "summary.txt"
    p.QUERY_FA = prefix + "OSK7121_supernova.psuedohap.10k.fasta.gz"
    p.REF_FA = prefix + "CF062B2D.contigs.fasta" 
    p.REF_BED =  prefix + "canu_data/CF062B2D.unitigs.corrected.bed"
    p.OUTPUT_DIR = prefix + "out"
    p.HG38_INDEX = "/media/scott/HDD/sickkids/hg38.mmi"
    p.HG38= "/media/scott/HDD/sickkids/hg38.fa"
    return p

def example_hybrid_params():
    
    p = Parameters()

    #override default parameters
    prefix = "/media/scott/Zapdos/CF062_19/"
    p.SUMMARY = prefix + "summary.txt"
    p.QUERY_FA = prefix + "OSK7121_supernova.psuedohap.10k.fasta.gz"
    p.REF_FA = prefix + "CF062B2D.contigs.fasta" 
    p.REF_BED =  prefix + "canu_data/CF062B2D.unitigs.corrected.bed"
    p.OUTPUT_DIR = prefix + "out2"
    p.HG38_INDEX = "/media/scott/HDD/sickkids/hg38.mmi"
    p.HG38= "/media/scott/HDD/sickkids/hg38.fa"
    p.VERBOSE=3

    return p


def explicit_polish_params():
    
    p = Parameters()

    #override default parameters
    prefix = "/media/scott/Zapdos/hybrid_project/polish_test/"
    p.FASTA = prefix + "scaff134.fasta"
    p.REF_ALIGNED_READS = prefix + "scaff134.pacbio.bam"
    p.QUERY_ALIGNED_READS = prefix + "scaff134.10x.bam"
    p.OUTPUT_DIR = prefix + "out"
    p.HG38_INDEX = "/media/scott/HDD/sickkids/hg38.mmi"
    p.HG38 = "/media/scott/HDD/sickkids/hg38.fa"

    return p

def explicit_graph_params():
    
    p = Parameters()
    
    p.TSV_INPUT = "/home/scott/Dropbox/hybrid-pipeline/example/input.txt"
    p.OUTPUT_DIR = "/media/scott/Zapdos/out/all/out"
    p.REF_GENOME = "/media/scott/Zapdos/out/all/chr5_polish/chr5_393462_677667.fasta"
    p.TARGET_REGION = "chr5_393462_677667:1-100"
    p.FILTER_X = 200
    p.COMBINE_Y = 20000
    p.SKIP_EXTRACTION = False

    return p