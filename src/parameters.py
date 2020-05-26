import argparse

class Parameters:
    #default parameters established here
    def __init__(self):
        
        self.OUTPUT_DIR   =   "./out"     #output directory
        self.VERBOSE            =   1       #verbosity
        self.WAIT               =   False    #wait for input at high verbosity 

        self.HG38 = None
        self.HG38_INDEX = None

        #hybrid specific parameters
        self.SUMMARY = None
        self.QUERY_FA     =   None
        self.REF_FA       =   None

        self.REF_BED = None

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
        self.TRUST_REF_CHUNKS   =   False    #trust the reference at the chunk level

        # polish specific parameters
        self.FASTA = None
        self.REF_ALIGNED_READS = None
        self.QUERY_ALIGNED_READS = None
        self.TARGET_CONTIG = None
        self.KEEP_INTERMEDIATE = False
        self.ANALYSIS = False
        
        # ref polish specific parameters
        self.ID = "_id_missing_"
        self.DOWNSAMPLE = 1

def get_parameters_hybrid():
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()

    #override default parameters
    prefix = "/media/scott/Zapdos/CF062_19/"
    p.SUMMARY = prefix + "summary.txt"
    p.QUERY_FA = prefix + "OSK7121_supernova.psuedohap.10k.fasta.gz"
    p.REF_FA = prefix + "CF062B2D.contigs.fasta" 
    p.REF_BED =  prefix + "canu_data/CF062B2D.unitigs.bed"
    p.OUTPUT_DIR = prefix + "out"

    p.HG38_INDEX = "/media/scott/HDD/sickkids/hg38.mmi"
    p.HG38= "/media/scott/HDD/sickkids/hg38.fa"


    parser = argparse.ArgumentParser(description="Hybrid assembly tool")

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


    args = parser.parse_args()
    #set parameters from user input
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

def get_parameters_polish():
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()

    #override default parameters
    prefix = "/media/scott/Zapdos/"
    p.FASTA = prefix + "scaff42.fasta"
    
    p.REF_ALIGNED_READS = prefix + "pb_scaff42.bam"
    p.QUERY_ALIGNED_READS = prefix + "10x_scaff42.bam"
    p.OUTPUT_DIR = prefix + "out"

    p.HG38_INDEX = "/media/scott/HDD/sickkids/hg38.mmi"
    p.HG38= "/media/scott/HDD/sickkids/hg38.fa"

    parser = argparse.ArgumentParser(description="Hybrid polish tool")

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

    parser.add_argument("-v", "--verbose", type=int, default=p.VERBOSE,
                        help="Verbosity. Default=" + str(p.VERBOSE))
    parser.add_argument("--wait", type=bool, default=p.WAIT,
                    help="Wait for user input at high verbosity levels. Default=False")
    parser.add_argument("-t", "--tig", type=str, default=p.TARGET_CONTIG,
                help="Target contig." )
    parser.add_argument("-a", "--analysis", type=bool, default=p.ANALYSIS,
                        help="Verbosity. Default=" + str(p.VERBOSE))

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


def get_parameters_reference_polish(CFID="CF002"):
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()
    '''

    p.REF_FA =  "/media/scott/HDD/sickkids/slc9a3/hg38.fa"
    p.REF_ALIGNED_READS = "/media/scott/HDD/sickkids/slc9a3/CF001.slc9a3.bam"
    p.OUTPUT_DIR = "/media/scott/HDD/sickkids/slc9a3/out"
    
    '''
    
    p.ID = CFID
    p.OUTPUT_DIR = "/media/scott/Zapdos/slc9a3_polish/out/" + CFID

    p.REF_FA =  "/media/scott/Zapdos/reference/hg38.fa"
    p.READS = "/media/scott/Zapdos/slc9a3_polish/" + CFID + "_slc9a3.bam"
    p.OUTPUT_DIR = "/media/scott/Zapdos/slc9a3_polish/out/" + CFID

    #p.REF_REGION = "chr5:393462-696129"
    p.REF_REGION = "chr5:393462-677667"

    parser = argparse.ArgumentParser(description="Reference sequence region polish")

    # Positional
    parser.add_argument("refFa", metavar="reference", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    parser.add_argument("reads", metavar="reads", default=p.READS, nargs="?",
                    help="Pacbio reads")
    
    parser.add_argument("sampleId", metavar="id", default=p.ID, nargs="?",
             help="ID for sample")

    parser.add_argument("refRegion", metavar="region", default=p.REF_REGION, nargs="?",
                help="Region to polish, formatted as chrom:start-end")

    # Optional
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                    help="Directory where output will be written." )
    
    parser.add_argument("-x", "--downsample", type=float, default=p.DOWNSAMPLE,
                    help="Randomly sample this propotion of the reads." )
    

    # Parse
    args = parser.parse_args()

    p.REF_FA = args.refFa
    p.REF_REGION = args.refRegion
    p.READS = args.reads
    p.ID = args.sampleId
    p.DOWNSAMPLE = args.downsample

    p.OUTPUT_DIR = args.outdir
    
    return p
