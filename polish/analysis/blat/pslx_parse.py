#!/usr/bin/env python3

print("Parsing psl file")

from Bio import SearchIO
import argparse


# Parse the input arguements:
parser = argparse.ArgumentParser(description = "counts the number of top BLAT matches to each chromosome")
parser.add_argument("-pslx", help = "pslx file of the BLAT search")
parser.add_argument("-n", help = " number of nucleotides per query to BLAT")
args = parser.parse_args() 


# a file containing only the top hits, in plsx format
toppslx_dir = args.pslx + "_tophits.pslx"
pslx_top = open(toppslx_dir, "w+")

pslx_top.write("matches" + "\tmisMatches" + "\trepMatches" + "\tnCount" + 
               "\tqNumInsert" + "\tqBaseInsert" + "\ttNumInsert" + 
               "\ttBaseInsert" + "\tstrand" + "\tqName" + "\tqSize" + 
               "\tqStart" + "\tqEnd" + "\ttName" + "\ttSize" + "\ttStart" + 
               "\ttEnd" + "\tblockCount" + "\tblockSizes" + "\tqStarts" + 
               "\ttStarts" + "\tqSeq" + "\ttSeq")


# a file with SNP information to input into VEP
vep_dir = "vep_input.txt"
vep = open(vep_dir, "w+")
#vep columns ("chr \tstart \tend \tallele \tstrand \tquery_id")


for results in SearchIO.parse(args.pslx, 'blat-psl', pslx=True):
        
    top_score = 0
    for hits in results.hsps:
        # if hsp has higher score, replace top_hit placeholder
        if hits.score > top_score:
            top_score = hits.score
            top_hit = hits
        # if hsp has same score and less fragments, use new hsp instead
        elif hits.score == top_score:
            if len(hits.fragments) < len(top_hit.fragments):
                top_score = hits.score
                top_hit = hits                
            
            
    # Rewrite the top hit into another pslx file
    tsize = results[top_hit.hit_id].seq_len
    q_frag_len = int(args.n)
    block_num = len(top_hit.fragments)
    t_starts = top_hit.hit_start_all    
    mismatch = top_hit.mismatch_num
    target = top_hit.hit_id
    chr_num = target[3:]
    
    if block_num == 1:
        strand = top_hit.query_strand
        query_seq = str(top_hit.query._seq) + ","
        hit_seq = str(top_hit.hit._seq) + ","
        qstarts = str(top_hit.query_start_all)[1:-1] + ","
                                               
    else:
        strand = 0
        query_seq = ""
        hit_seq = ""
        
        for frag in top_hit:
            strand += frag.query_strand
            query_seq += str(frag.query._seq) + ","
            hit_seq += str(frag.hit._seq) + ","
        strand = int(strand / block_num)
        
        if strand == -1:
            qstarts = []
            for q in top_hit.query_end_all:
                qstarts.append(q_frag_len - q)
            qstarts = str(qstarts)[1:-1] 
        else:
            qstarts = top_hit.query_start_all
            qstarts = str(qstarts)[1:-1] + ","  
        
    if strand == 1:
        strand = '+'
    else:
        strand = '-'
        
    # if the alignment contains a mismatch 
    if mismatch != 0:
        qseq_list = query_seq.split(',')[:block_num]
        tseq_list = hit_seq.split(',')[:block_num]
        
        # for every fragment with a mismatch
        for q, t, pos in zip(qseq_list, tseq_list, top_hit.hit_start_all):
            
            # for each nt in fragment alignment
            for q_nt, t_nt in zip(q, t):
                
                # +1 to pos before SNP call because BLAT was 0-based and VEP uses 1-based
                pos += 1 
                if q_nt != t_nt:
                    snp_pos = pos
                    allele = t_nt + "/" + q_nt
                    # write out SNP position for VEP input
                    contig_name = top_hit.query_id[:top_hit.query_id.index("|")]
                    contig_pos = top_hit.query_id[len(top_hit.query_id) - 
                                                  top_hit.query_id.index("_")-1:]
                    vep.write("\n{}\t{}\t{}\t{}\t{}\t{}".format(
                        chr_num, snp_pos, snp_pos, allele, strand, 
                        contig_name + ":" + contig_pos + "_" + allele))  
               
        
    pslx_top.write("\n{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}"
                    "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}"
                     "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}" "\t{}".format(
                    top_hit.match_num, mismatch, 
                    top_hit.match_rep_num, top_hit.n_num, 
                    top_hit.query_gapopen_num, top_hit.query_gap_num, 
                    top_hit.hit_gapopen_num, top_hit.hit_gap_num, strand, 
                    top_hit.query_id, q_frag_len , top_hit.query_start, 
                    top_hit.query_end, target, tsize, top_hit.hit_start, 
                    top_hit.hit_end, len(top_hit.query_span_all), 
                    str(top_hit.query_span_all)[1:-1]+",", qstarts, 
                    str(t_starts)[1:-1]+",", query_seq, hit_seq))

    
pslx_top.close()
vep.close()