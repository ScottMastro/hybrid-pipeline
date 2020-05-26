#---------------------------------------------
library("karyoploteR")
library("regioneR")
library("GenomicRanges")
library("stringr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("GenomicRanges")


#source("http://bioconductor.org/biocLite.R")
#biocLite("IRanges")
#library("BiocManager")
#install("BSgenome.Hsapiens.UCSC.hg38")
genome <- getGenomeAndMask("hg38")$genome
genome <- filterChromosomes(genome, chr.type="canonical")
chrs = as.vector(genome@seqnames@values)
width = genome@ranges@width

get_alignments <- function(dir){
  alignments = read.csv(sep="\t", dir)
  alignments = alignments[! is.na(alignments$IDY), ]
  x = str_split(alignments$Reference, "_", n=2)
  chroms = unlist(lapply(x, `[[`, 1))
  alignments$chr = chroms
  alignments = alignments[! grepl("rand", alignments$Reference, fixed=T),]
  
  alignments$S1 = strtoi(alignments$S1)
  alignments$E1 = strtoi(alignments$E1)
  alignments$S2 = strtoi(alignments$S2)
  alignments$E2 = strtoi(alignments$E2)
  alignments <- alignments[order(alignments$chr, alignments$S1),]
  return (alignments)
}

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


s=205882176
e=205912588
buffer= 100000 #100kb
#buffer= 10000000 # 10MB
size=(s-e + buffer*2)

detail.region <- toGRanges(data.frame("chr1", s-buffer, e+buffer))
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes=c("chr1"), zoom=detail.region)
kpAddBaseNumbers(kp, tick.dist=size/5, minor.tick.dist=size/10, add.units=TRUE)
kpDataBackground(kp, data.panel = 1)

kpRect(kp, chr="chr1", 
       x0=s, x1=e,
       y0=0,y1=1, r0=0, r1=0.1,
       col="#ff6699")


startfrom=0.1

expected=17
step= (1-startfrom)/(expected)
gap=4
index=0
files <- list.files(path="C:/Users/scott/Dropbox/hybrid-pipeline/to_plot/quast_alignments", full.names=TRUE, recursive=FALSE)
for (file in files) {

  print(file)
  if (file == "C:/Users/scott/Dropbox/hybrid-pipeline/to_plot/quast_alignments/desktop.ini"){
    break
  }
  
  aln = get_alignments(file)
  aln=aln[aln$chr == "chr1",]
  
  phasedir = gsub("quast", "whatshap", file)
  
  
  kpRect(kp, chr="chr1", 
            x0=min(aln$S1), x1=max(aln$E1),
            y0=0,y1=1, r0=startfrom + index*step + (step/gap)/2, r1=startfrom + (index+1)*step - (step/gap)/2,
            col="#D06Ec4", border=NA)
  
  if(file.exists(phasedir)){
    phase=phaseblock_to_data(phasedir, aln)
    
    kpRect(kp, chr="chr1", 
           x0=phase$x1, x1=phase$x2,
           y0=0,y1=1, r0=startfrom + index*step + (step/gap)/2, r1=startfrom + (index+1)*step - (step/gap)/2,
           col="#D06Ec4")
  }

  index = index+1
}










plot_top <- function(x, c0, c1, c2){
  
  kp <- plotKaryotype(genome="hg38")
  kpDataBackground(kp, data.panel = 1)
  
  expected=20
  last=0
  
  files <- list.files(path="C:/Users/scott/Desktop/quast/quast_alignments", full.names=TRUE, recursive=FALSE)
  for (file in files) {
    print(file)
    aln = get_alignments(paste(sep="/", file, x))
    hdf = top_contigs(aln, n=50)
    
    kpHeatmap(kp, chr=hdf$chr, 
              x0=hdf$s, x1=hdf$e,
              y=hdf$value, r0=0.05 + last, r1=0.05 + last+0.75/expected,
              colors = c(c0,c1,c2))
    last = last+1.0/expected
  }
}

plot_top("hybrid_alignments.tsv", "#341752", "#D06Ec4", "#F3C6FF")
plot_top("supernova_alignments.tsv", "#052B52", "#99CCFF", "#9DD2FF")
plot_top("canu_alignments.tsv", "#764F13", "#FF8840", "#FFDDAA")


w=1e5
files <- list.files(path="C:/Users/scott/Desktop/quast/quast_alignments", full.names=TRUE, recursive=FALSE)
hdf = empty_data(w)
cdf = empty_data(w)
sdf = empty_data(w)
count = 0

for (x in files) {
  print(x)
  aln = get_alignments(paste(sep="/", x, "hybrid_alignments.tsv"))
  data = alignments_to_data(aln, w)
  hdf$value = hdf$value + data$value
  
  aln = get_alignments(paste(sep="/", x, "canu_alignments.tsv"))
  data = alignments_to_data(aln, w)
  cdf$value = cdf$value + data$value
  
  aln = get_alignments(paste(sep="/", x, "supernova_alignments.tsv"))
  data = alignments_to_data(aln, w)
  sdf$value = sdf$value + data$value
  
  count = count+1
}
hdf$value = hdf$value/count
cdf$value = cdf$value/count
sdf$value = sdf$value/count

kp <- plotKaryotype(genome="hg38")
kpDataBackground(kp, data.panel = 1)
#kpAddBaseNumbers(kp)
#kpPoints(kp, chr=data$chr, x=data$s, y=data$value)

kpHeatmap(kp, chr=hdf$chr, 
          x0=hdf$s, x1=hdf$e,
          y=pmin(hdf$value,1), r0=0, r1=0.2,
          colors = c("white", "purple"))

kpHeatmap(kp, chr=cdf$chr, 
          x0=cdf$s, x1=cdf$e,
          y=pmin(cdf$value,1), r0=0.3, r1=0.5,
          colors = c("white", "orange"))

kpHeatmap(kp, chr=sdf$chr, 
          x0=sdf$s, x1=sdf$e,
          y=pmin(sdf$value,1), r0=0.6, r1=0.8,
          colors = c("white", "blue"))

kpBars(kp, chr=data$chr, x0=data$s, x1=data$e, y1=data$value, col="#000000", border=NA)
kpRect(kp, chr=data$chr, x0=data$s, x1=data$e, y0=0,y1=1, alpha=data$value, col="#000000", border=NA)


kp <- plotKaryotype(genome="hg38")
kpDataBackground(kp, data.panel = 1)

kpBars(kp, chr=hdf$chr, 
       x0=hdf$s, x1=hdf$e,
       y1=pmin(hdf$value,1), r0=0, r1=0.2,
       col="purple", border=NA)

kpBars(kp, chr=sdf$chr, 
       x0=sdf$s, x1=sdf$e,
       y1=pmin(sdf$value,1), r0=0.3, r1=0.5,
       col="blue", border=NA)
kpBars(kp, chr=cdf$chr, 
       x0=cdf$s, x1=cdf$e,
       y1=pmin(cdf$value,1), r0=0.6, r1=0.8,
       col="orange", border=NA)

#Select some data points as "special ones"
sel.data <- rand.data[c(7, 30, 38, 52),] 
head(rand.data)

pp <- getDefaultPlotParams(plot.type = 1)

tr.i <- 1/11
tr.o <- 1/10

kp <- plotKaryotype( plot.params = pp) 

dd <- toGRanges(data.frame(chr="chr1", start=end(kp$genome[1])/50*(0:49), end=end(kp$genome[1])/50*(1:50)))
mcols(dd) <- data.frame(y=((sin(start(dd)) + rnorm(n=50, mean=0, sd=0.1))/5)+0.5)
tn <- 2
kpDataBackground(kp)
kpBars(kp, dd, y1=dd$y, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#000000", border=NA)
kpRect(kp, chr="chr1", x0=5000000, x1=45000000, y0=0.2, y1=0.8, r0=tr.o*tn, r1=tr.o*tn+tr.i, col="#EEEEEE")
kpText(kp, chr="chr1", x=25000000, y=0.5, col="red", r0=tr.o*tn, r1=tr.o*tn+tr.i, labels="kpBars", cex=0.7)