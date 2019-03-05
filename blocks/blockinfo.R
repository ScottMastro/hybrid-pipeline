library(ggplot2)
library(reshape2)
#source("helper.R")
dir = "/home/scott/Dropbox/hybrid-pipeline/blocks"

df = read.csv(paste(dir,"blockinfodf.txt",sep="/"), sep="\t", header = F)

col.names <- c("canu", "nova", "csize", "nsize", "cstart", "cend", "nstart", "nend", "ccov", "cmiss",
          "coverlap", "ncov", "nmiss", "noverlap")

smallcut <- 50000
covcut <- 1.5
partial <- 0.05

colnames(df) <- col.names
df$nova <- as.character(df$nova)
df$size.ratio <- df$csize/df$nsize
df$cov.ratio <- (df$ccov + df$cmiss)/(df$ncov + df$nmiss)

canus <- as.character(unique(df$canu))
novas <- as.character(unique(df$nova))

canu<- canus[15]

filterBlocks <- function(canu){
  df.canu <- df[df$canu %in% c(canu),]
  
  df.canu$filter <- "ok"
  df.canu[df.canu$nsize < smallcut,]$filter <- "small" 
  
  #filter by coverage length
  df.canu[abs(log(df.canu$cov.ratio)) > log(covcut),]$filter <- "coverage"
  
  
  
  p <- ggplot(df.canu, aes(x=novaCount, y=(cstart+cend)/2, colour=nsize)) +
    geom_point(size=0.1, alpha=0.3)+
    geom_errorbar(aes(ymin = cstart, ymax = cend), size=1, width=0, alpha=0.5) +
    xlab('log10(nsize)') + 
    ylab('Position (0 to 1)') + 
    ggtitle(sapply(canu, paste, collapse=", ")) +
    coord_flip() +
    facet_wrap(~(filter))
  #p <- scottTheme(p)

  p

}













df$cov.ok <- abs(log(df$cov.ratio)) < log(covcut)
df$is.big <- df$ccov > smallcut & df$ncov > smallcut


canuCount<-table(unlist(df[df$cov.ok,]$canu))
novaCount<-table(unlist(df[df$cov.ok,]$nova))
df$canuCount <- sapply(df$canu, function(x) as.integer(canuCount[x]))
df$novaCount <- sapply(df$nova, function(x) as.integer(novaCount[x]))
df[is.na(df$canuCount),]$canuCount <- 0
df[is.na(df$novaCount),]$novaCount <- 0

df.filtered <- df[df$cov.ok & df$is.big,]

canu<- canus[15]
df.canu <- df[df$canu %in% c(canu),]
pltCanu <- function(canu){
  p <- ggplot(df[df$canu %in% c(canu),], aes(x=novaCount, y=(cstart+cend)/2, colour=nsize)) +
    geom_point(size=0.1, alpha=0.3)+
    geom_errorbar(aes(ymin = cstart, ymax = cend), size=1, width=0, alpha=0.5) +
    xlab('log10(nsize)') + 
    ylab('Position (0 to 1)') + 
    ggtitle(sapply(canu, paste, collapse=", ")) +
    coord_flip() +
    facet_wrap(~(cov.ok & is.big))
  #p <- scottTheme(p)
  return (p)
}
pltCanu(canu)

nova<- "16"# novas[5]
df.nova <- df[df$nova %in% c(nova),]
pltNova <- function(nova){
  p <- ggplot(df[df$nova %in% c(nova),], aes(x=canuCount, y=(nstart+nend)/2, colour=csize)) +
    geom_point(size=0.1, alpha=0.3)+
    geom_errorbar(aes(ymin = nstart, ymax = nend), size=1, width=0, alpha=0.5) +
    xlab('log10(csize)') + 
    ylab('Position (0 to 1)') + 
    ggtitle(sapply(nova, paste, collapse=", ")) +
    coord_flip() +
    facet_wrap(~(cov.ok & is.big))
  
  #p <- scottTheme(p)
  return (p)
}
pltNova(nova)



plot(log(df$csize), log(df$nsize))
plot(log(df.filtered$csize), log(df.filtered$nsize))
plot(log(df.filtered$coverlap), log(df.filtered$ccov))



plot(df.big$cov.ratio, (df.big$ccov + df.big$ncov)/2)
hist(log(df$cov.ratio), breaks = 100)
hist(log(df.big$size.ratio))
hist(df.big$csize)
hist(df.big$nsize)