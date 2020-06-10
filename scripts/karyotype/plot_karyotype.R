source("karyotype_helper.R")

quast.directory="D:/out/quast"
#quast.directory="C:/Users/scott/Desktop/quast/quast_alignments"

supernova.target="*pseudohap.tsv"
canu.target="*contigs.tsv"
hybrid.target="*hybrid_assembly.tsv"

quast.files <- list.files(path=quast.directory, full.names=TRUE, recursive=FALSE)
NSAMPLES=32




#---------------------------------------------
#GENOME PLOT
#---------------------------------------------

genome_plot <- function(targets, nrows, colors, row.spacer=0, max.contigs=200){
  
  kp <- plotKaryotype(genome="hg38")
  kpDataBackground(kp, data.panel = 1, col="white")

  gaps = get_hg38_gaps()
  gap.data = makeGRangesFromDataFrame(get_hg38_gaps())
  kpPlotRegions(kp, col="#dddddd", data=gap.data, avoid.overlapping=F)
  
  total.rows = nrows*length(targets) + row.spacer
  y=0
  
  for (i in 1:length(targets)){
    target = targets[i]
    c0=colors[(i-1)*3 + 1]
    c1=colors[(i-1)*3 + 2]
    c2=colors[(i-1)*3 + 3]
    
    dirs <- list.files(path=quast.directory, full.names=TRUE, recursive=FALSE)
    for (dir in dirs) {
      
      print(dir)
      align.df = get_quast_alignments(dir, target)
      if (length(align.df) == 1 && is.na(align.df)){ next }
      
      hdf = top_contigs(align.df, n=max.contigs)
      
      kpHeatmap(kp, chr=hdf$chr, 
                x0=hdf$s, x1=hdf$e,
                y=hdf$value, r0=0.05 + y, r1=0.05 + y+0.75/total.rows,
                colors = c(c0,c1,c2))
      y = y + 1.0/total.rows
      
    }

  }
}


purples=c("#341752", "#D06Ec4", "#FCF0FF")
blues=c("#052B52", "#1e00d6", "#9DD2FF")
oranges=c("#764F13", "#FF8840", "#FFDDAA")

targets=c(hybrid.target, canu.target, supernova.target)
colors =c(purples, oranges, blues)
#genome_plot(targets, NSAMPLES, colors)
mc=100
rs=10
genome_plot(hybrid.target, NSAMPLES, purples, max.contigs=mc, row.spacer=rs)
genome_plot(canu.target, NSAMPLES, oranges, max.contigs=mc, row.spacer=rs)
genome_plot(supernova.target, NSAMPLES, blues, max.contigs=mc, row.spacer=rs)





kp <- plotKaryotype(genome="hg38")
kpDataBackground(kp, data.panel = 1, col="white")
expected=nrows
last=0
i=0

gaps = get_hg38_gaps()
gap.data = makeGRangesFromDataFrame(get_hg38_gaps())
kpPlotRegions(kp, col="#dddddd", data=gap.data, avoid.overlapping=F)



dirs <- list.files(path=quast.directory, full.names=TRUE, recursive=FALSE)
for (dir in dirs) {
  
  print(dir)
  align.df = get_quast_alignments(dir, target)
  if (length(align.df) == 1 && is.na(align.df)){ next }
  
  hdf = top_contigs(align.df, n=max.contigs)
  
  kpHeatmap(kp, chr=hdf$chr, 
            x0=hdf$s, x1=hdf$e,
            y=hdf$value, r0=0.05 + last, r1=0.05 + last+0.75/expected,
            colors = c(c0,c1,c2))
  last = last+1.0/expected
  i=i+1
  print(i)
  
}









kp <- plotKaryotype(plot.type = 2, chromosomes = c("chr1", "chr2", "chr3"))

### Data Panel 1 ###

#Big regions
kpRect(kp, data = big.regs.up, y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=0.8)
kpRect(kp, data = big.regs.down, y0=0, y1=1, col="#DDFFDD", border=NA, r0=0, r1=0.8)

#Data points
kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=0.8, numticks = 5, col="#666666", cex=0.5)
kpPoints(kp, data=data.points, pch=16, cex=0.5, col=dp.colors, r0=0, r1=0.8)

#Mean and sd of the data points.  
for(chr in seqlevels(kp$genome)) {
  chr.dp <- sort(keepSeqlevels(x = data.points, value = chr, pruning.mode = "coarse"))
  rmean <- rollmean(chr.dp$y, k = 6, align = "center")  
  rsd <- rollapply(data = chr.dp$y, FUN=sd, width=6)
  kpLines(kp, chr = chr, x=start(chr.dp)[3:(length(chr.dp)-3)], y=rmean, col=data.points.colors[3], r0=0, r1=0.8)
  kpPlotRibbon(kp, chr=chr, data=chr.dp[3:(length(chr.dp)-3)], y0=rmean-rsd, y1=rmean+rsd, r0=0, r1=0.8, col="#FF336633", border=NA)
}

#Markers
kpPlotMarkers(kp, data=marks, label.color = "#333333", r1=1.1, cex=0.5, label.margin = 5)

### Data Panel 2 ###

#medium regions and their coverage

kpPlotRegions(kp, data = mid.regs, r0 = 0.2, r1=1, border=NA, data.panel=2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0, col=data.points.colors[2], data.panel = 2)
kpPlotCoverage(kp, data=mid.regs, r0=0.2, r1=0.12, col=data.points.colors[1], data.panel = 2)

kpText(kp, chr=seqlevels(kp$genome), y=0.4, x=0, data.panel = 2, r0=0.2, r1=0, col="#444444", label="30x", cex=0.8, pos=2)
kpAbline(kp, h=0.4, data.panel = 2, r0=0.2, r1=0, col=data.points.colors[3])

plot of chunk Figure











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




#---------------------------------------------
#SLC26A9 PLOT
#---------------------------------------------

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

for (file in quast.files) {

  print(file)
  align.df = get_quast_alignments(file, hybrid.target)
  if (length(align.df) == 1 && is.na(align.df)){ next }
  
  align.df = align.df[align.df$chr == "chr1",]
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


#---------------------------------------------
#GENOME PLOT
#---------------------------------------------

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