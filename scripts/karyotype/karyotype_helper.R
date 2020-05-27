#---------------------------------------------
library("karyoploteR")
library("regioneR")
library("GenomicRanges")
library("stringr")

genome <- getGenomeAndMask("hg38")$genome
genome <- filterChromosomes(genome, chr.type="canonical")
chrs = as.vector(genome@seqnames@values)
width = genome@ranges@width


empty_data <- function(w=1e5){
  
  split_chr <- function(i){
    chr = chrs[i]
    len = width[i]
    windows = seq(1, len, by=w)
    return(data.frame(chr=chr, s=windows,e=windows+w-1, value=0))
  }
  df = bind_rows(lapply(X=1:length(chrs), FUN=split_chr))
  return(df)
}

top_contigs <- function(aln, n=50){
  
  sorted <- aln[order(-aln$E2),]
  tigs = unique(as.vector(sorted$Contig))
  top = tigs[1:n]
  filtered = aln[aln$Contig %in% top,]
  
  df = data.frame(chr=filtered$chr, s=filtered$S1,e=filtered$E1, 
                  value=match(filtered$Contig,top))
  return(df)
}

alignments_to_data <- function(aln, w=1e5){
  df = empty_data(w)
  grw = with(df, GRanges(chr, IRanges(start=s, end=e)))
  gra = with(aln, GRanges(chr, IRanges(start=S1, end=E1)))
  hits = findOverlaps(query=gra, subject=grw)
  overlaps <- pintersect(gra[queryHits(hits)], grw[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(grw[subjectHits(hits)])
  agg = aggregate(percentOverlap, by=list(subjectHits(hits)), FUN=sum)
  df[agg$Group.1,"value"] = agg$x
  return(df)
}

phaseblock_to_data <- function(phasedir, aln){
  phase = read.csv(phasedir)
  
  map_pos2 <- function(x, aln){
    
    start = min(aln$S1)
    end = max(aln$E1)
    tiglen=max(max(aln$S2), max(aln$E2))
    
    return( (x/tiglen)*(end-start) + start )
  }
  
  
  map_pos <- function(x, aln){
    
    idx1 = which.min(abs(aln$S2 - x))
    idx2 = which.min(abs(aln$E2 - x))
    dist1 = aln$S2[idx1] - x
    dist2 = aln$E2[idx2] - x
    
    a=aln$S2[idx1]
    r=aln$S1[idx1]
    strand=1
    if(aln$S2[idx1] > aln$E2[idx1]){
      #r=aln$E1[idx1]
      strand=-1
    }
    
    if(abs(dist2) < abs(dist1)){ 
      a=aln$E2[idx2]
      r=aln$E1[idx2]
      strand=1
      if(aln$S2[idx2] > aln$E2[idx2]){
        #r=aln$S1[idx2]
        strand=-1
      }
    }
    
    return( r + (x-a)*strand )
    
  }
  
  phase$x1=0
  phase$x2=0
  for (i in 1:nrow(phase)){
    
    phase[i,"x1"]=  map_pos2(phase[i,"from"], aln)
    phase[i,"x2"]= map_pos2(phase[i,"to"], aln)
    
  }
  
  return(phase)
}