#!/usr/bin/Rscript

## This script is a modified version of 
### https://raw.githubusercontent.com/ajsedgewick/epidermis-scRNA/master/src/run_zinb.R
### Which was used in the publication :
###  http://doi.org/10.1016/j.celrep.2018.09.006

args = commandArgs(trailingOnly=TRUE)
if(length(args) !=3){
    print("Usage: run_zinb.R countdata.tsv coldata.tsv output_directory")
    quit("no", status=1)
}

library(Seurat)
library(Matrix)
library(SummarizedExperiment)
library(zinbwave)
library(BiocParallel)
library(edgeR)
library(feather)
library(tibble)


countdata <- read.csv(args[1], row.names=1, header = T)  ## each column is a cell, each row is a gene
countdata <- as.matrix(countdata)
colnames(countdata) <- gsub("\\." , "-" , colnames(countdata) )
coldata <- read.table(args[2], header=T, row.names=1, check.names=F, sep=",")
outdir <- args[3]

#cpm <- cpm(countdata, log = F, prior.count = 0)
#is_quality <- rowSums(cpm  >= 500 ) >= .001*ncol(cpm) # equivalent of 5 cp10k

is_quality <- rowSums(countdata >= 5) >= 100 #5 umi in 100 cells
sum(is_quality)

genes.use <- rownames(countdata)[is_quality]
countdata <-  countdata[genes.use,]
dim(countdata)
dim(coldata)
#coldata$lognUMI <- log2(coldata$nUMI + 1)

cur.se <- SummarizedExperiment(countdata,
                               colData=coldata)

cores <- as.numeric(Sys.getenv("NSLOTS"))
register(MulticoreParam(workers=as.integer(cores)))
			       
zbc <- zinbFit(cur.se, K=20, X="~nUMI + sample + percent.mito", epsilon=1000 )
save(zbc, file=paste(outdir, "zinb_fit.RData", sep="/"))

se_norm <- zinbwave(cur.se, fitted_model=zbc, normalizedValues=TRUE,
                    residuals = TRUE,)
save(se_norm, file=paste(outdir, "zinb_wavese.RData", sep="/"))
W <- getW(zbc)
colnames(W) <- paste0("W", 1:20)
rownames(W) <- colnames(assays(se_norm)$normalizedValues)
write.table(round(W, 5), file=paste(outdir, "zinb_W.csv", sep="/"),
            sep=",", col.names=NA, row.names=T, quote=F)

write.table(genes.use, file=paste(outdir, "genes_use.txt", sep="/"), col.names=F, row.names=F, quote=F)
