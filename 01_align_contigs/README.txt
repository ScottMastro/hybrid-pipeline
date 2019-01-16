1. align_contigs.sh
 -constructs blast database of canu assembly
 -splits supernova sequence into x bp length segments (1000)
 -blasts each supernova segment against canu database, triggers multiple PBS jobs

2. summarize_hits.sh
 -combines all blast result files
 -filters blast results (cov=40, identity=95)
 -summarizes blast results
 -constructs blocks of contiguous hits

