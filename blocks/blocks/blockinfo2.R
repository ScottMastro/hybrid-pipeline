library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)

#source("helper.R")
dir = "/home/scott/Dropbox/hybrid-pipeline/blocks"

df = read.csv(paste(dir,"blockinfodf.txt",sep="/"), sep="\t", header = F)

col.names <- c("canu", "nova", "csize", "nsize", "cstart", "cend", "nstart", "nend", "ccov", "cmiss",
               "coverlap", "ncov", "nmiss", "noverlap")

sizeThreshold = 10000
edgeThreshold = 0.02
initialCovThreshold <- 2

colnames(df) <- col.names
df$nova <- as.character(df$nova)
df$size.ratio <- df$csize/df$nsize
df$cov.ratio <- (df$ccov + df$cmiss)/(df$ncov + df$nmiss)
df$cmiss.percent <- df$cmiss/(df$ccov + df$cmiss)
df$nmiss.percent <- df$nmiss/(df$ncov + df$nmiss)
df$ccov.percent <- (df$ccov + df$cmiss)/df$csize
df$ncov.percent <- (df$ncov + df$nmiss)/df$nsize
df = df[abs(log(df$cov.ratio)) < log(initialCovThreshold),]

canus <- as.character(unique(df$canu))
novas <- as.character(unique(df$nova))

novaSize<-table(unlist(df[df$filter,]$canu))


covThreshold <- 1.05



dff= list()
for (i in 1:length(novas)){
  
  i = round(runif(1, 1, length(novas)))
  
  
  nova = novas[i]
  df.nova <- df[df$nova %in% c(nova),]
  novaSize = df.nova$nsize[1] 
    
  if (novaSize >= sizeThreshold){
    df.fullmatch = df.nova[df.nova$nstart < edgeThreshold & df.nova$nend > (1 - edgeThreshold),]
    df.fullmatch = df.fullmatch[abs(log(df.fullmatch$cov.ratio)) < log(covThreshold),]
    
    if(nrow(df.fullmatch) == 1){
      dff[[i]] = df.fullmatch
      print("full")
    }
    else if (nrow(df.fullmatch) > 1){
      print("multiple")
    }
    
    
  }
  
 # print("skipping small contig")
}

dff = bind_rows(dff)
plot(log(dff$cov.ratio), log10(dff$cmiss), col = rgb(0, 0, log(dff$nsize)/log(max(dff$nsize))))
hist(log(dff$cov.ratio), breaks = 100)
hist(log(df.big$size.ratio))
hist(df.big$csize)
hist(df.big$nsize)




df.nova$agree.size = df.nova$ncov - df.nova$nmiss
maxAgree = max(df.nova$agree.size)
df.nova = df.nova[df.nova$agree.size > maxAgree/10,]









#class A nova contigs -> greater than 10 kb
sizeThreshold = 10000

dfa <- df[df$nsize >= sizeThreshold,]
dfb <- df[df$nsize < sizeThreshold,]


#class 1a nova contigs -> fully covered by single canu
edgeThreshold = 0.02
covcut <- 1.05
missingTreshold = 0.1


dfa$filter <- TRUE
dfa[dfa$nstart > edgeThreshold,]$filter <- FALSE
dfa[dfa$nend < (1 - edgeThreshold),]$filter <- FALSE

canuCount<-table(unlist(dfa[dfa$filter,]$canu))
novaCount<-table(unlist(dfa[dfa$filter,]$nova))
dfa$canuCount <- sapply(dfa$canu, function(x) as.integer(canuCount[x]))
dfa$novaCount <- sapply(dfa$nova, function(x) as.integer(novaCount[x]))
dfa[is.na(dfa$canuCount),]$canuCount <- 0
dfa[is.na(dfa$novaCount),]$novaCount <- 0

dfa[dfa$novaCount > 1,]$filter <- FALSE


dfa[abs(log(dfa$cov.ratio)) > log(covcut),]$filter <- FALSE
dfa[dfa$nmiss.percent > missingTreshold,]$filter <- FALSE
dfa[dfa$cmiss.percent > missingTreshold,]$filter <- FALSE

nrow(dfa[dfa$filter,])
nrow(dfa) - nrow(dfa[dfa$filter,])





dffa <- dfa[dfa$filter,]
hist(log10(dffa$cov.ratio), breaks=100)
#plot(log10(dffa$cov.ratio), log10(dffa$nsize))





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