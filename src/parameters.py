import sys
import argparse

SOFTWARE_VERSION="v1.0.0"
SOFTWARE_NAME="Hybrid Assembly Pipeline"
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

#==================================================
# CLI help 
#==================================================

def get_parameters(defaultParameters=None):

    p = Parameters()

    parser = argparse.ArgumentParser(
        prog=SOFTWARE_NAME,
        description=f"A hybrid assembly tool for integrating supernova and canu assemblies. {SOFTWARE_VERSION}"
    )
    
    parser.add_argument("--version", 
        action="store_true",
        help="Show program's version number and exit."
    )

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

    return set_parameters(parser, defaultParameters)

#==================================================
# CLI parsing 
#==================================================

def set_parameters(parser, defaultParameters=None):

    p = Parameters() if defaultParameters is None else defaultParameters

    args = parser.parse_args(sys.argv[1:])
    
    if args.version:
        print(f"{SOFTWARE_NAME} {SOFTWARE_VERSION}")
        sys.exit(0)

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
    
    if not validate_parameters(p):
        parser.print_help()
        sys.exit(1)
    
    return p

#==================================================
# Parameter validation 
#==================================================

def validate_parameters(p):

    if p.QUERY_FA is None: return False
    if p.REF_FA   is None: return False
    if p.SUMMARY  is None: return False
    
    return True

