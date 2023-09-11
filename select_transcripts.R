#!usr/bin/Rscript
args <- commandArgs(trailingOnly=T)

if(length(args) < 1){
    stop("Usage: Rscript select_transcripts.R readcounts.tab > transcripts.list")
}

file <- args[1]

geomean <- function(x,...){
    exp(mean(log(x),...))
}

df = read.table(file,sep="\t",col.names=c("sID","gID","tID","count"))
df <- df[df$count > 0,]

pseudoref <- sapply(split(df$count,df$tID),geomean)
df$ratio = df$count/pseudoref[df$tID]
normfactor <- sapply(split(df$ratio,df$sID),median)
df$norm <- df$count / normfactor[df$sID]

gmean <- sort(sapply(split(df$count,df$tID),geomean))
df$rank <- match(df$tID,names(gmean))

sCount <- length(unique(df$sID))
ValidGenes <- names(which(sapply(sapply(split(df$sID,df$gID),unique),length) == sCount))
df <- df[df$gID %in% ValidGenes,]

BestTIDs <- sapply(split(df[c("tID","rank")],df$gID),function(x){x$tID[which.max(x$rank)]})
cat(paste(BestTIDs,collapse="\n"))

