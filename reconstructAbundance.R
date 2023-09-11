##Declarations
library(Rphylopars)
library(phytools)

args <- commandArgs(trailingOnly = T)

TAbFile <- "results/readcounts.tab"
treeFile <- "metadata/tree.newick"

##Functions

geomean <- function(x,...){
    exp(mean(log(x),...))
}

LoadTAb <- function(file,tIDList=NULL){
    df <- read.table(file,sep="\t",stringsAsFactors = F,col.names=c("sID","gID","tID","count"));
    pseudoref <- sapply(split(df$count,df$tID),geomean)
    df$ratio = df$count/pseudoref[df$tID]
    normfactor <- sapply(split(df$ratio[!is.na(df$ratio)],df$sID[!is.na(df$ratio)]),median)
    df$norm <- df$count / normfactor[df$sID]
    tmp <- split(df$norm,df$tID)
    tmp2 <- split(df$sID,df$tID)
    if(is.null(tIDList)){
        tIDList <- unique(df$tID)
    }
    sIDList <- unique(df$sID)
    trait_data <- sapply(tIDList,function(nm){tmp[[nm]][match(sIDList,tmp2[[nm]])]})
    rownames(trait_data) = sIDList
    trait_data
}

anc.recon.transform <- function(trait_data,tree,transform=identity,...){
    trait_data <- transform(trait_data)
    trait_data <- trait_data[,apply(trait_data,2,function(x){sum(is.na(x) | is.infinite(x))}) == 0]
    res <- anc.recon(trait_data,tree,...)
    recon <- res
    if(length(list(...))){
        recon <- res$Yhat
        rownames(recon) <- tree$node.label
        recon <- rbind(trait_data,recon)
        res$Yhat = recon
    } else {
        rownames(res) <- tree$node.label
        res <- rbind(trait_data,res)
    }
    res
}

##Main

#stop("Functions Sourced")

TAb.mat <- LoadTAb(TAbFile)
treeObj <- read.newick(treeFile)
recon.log2 <- anc.recon.transform(TAb.mat,treeObj,log2)

recon <- round(2^(recon.log2))
message("Writing Normalized and Reconstructed values to results/Norm+ReconReadcounts.tab");
write.table(recon,"results/Norm+ReconReadcounts.tab",quote=F,sep="\t",row.names=T,col.names=T)
