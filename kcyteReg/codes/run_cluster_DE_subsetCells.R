#./run_normal_DE.R coundata.RData coldata.Rdata outdir
#OR ./run_normal_DE.R coundata.RData coldata.Rdata clusters.tsv outdir 
# where clusters.tsv is a table with cell names in first column and cluster number in the second

## This code is modified from code by AJ Sedgewick ( https://github.com/ajsedgewick/epidermis-scRNA/blob/master/src/DE_functions.R )

library(limma)
library(openxlsx)
library(clusterExperiment)

########################################################################
## FUNCTIONS 

clustDE <- function(coldata, countdata, outfile="clusterDE.xlsx"){
  #coldata$cluster <- as.factor(cluster)
  coldata$sample <- droplevels(coldata$sample)
  coldata$cluster <- droplevels(coldata$cluster)
  cur.k <- length(levels(coldata$cluster))
  
  #design <- model.matrix(~0 + cluster + sample + percent.mito + nUMI, data=coldata)
  design <- model.matrix(~0 + cluster + sample + percent.mito, data=coldata)
  #vm <- voom(countdata, design)
  #fit <- lmFit(vm, design)
  y <- new("EList")
  y$E <- edgeR::cpm(countdata, log = TRUE, prior.count = 1)
  fit <- lmFit(y, design)


  clCont <- rbind(clusterContrasts(coldata$cluster, "OneAgainstAll")$contrastMatrix, 
                  matrix(0, nrow=dim(design)[2] - cur.k, ncol=cur.k))
  rownames(clCont) <- colnames(design)
  colnames(clCont) <- colnames(design)[1:cur.k]
  
  fit <- contrasts.fit(fit, clCont)
  #fit <- eBayes(fit)
  fit <- eBayes(fit, trend = T, robust=T)

  clust.tt <- lapply(colnames(clCont), function(x) topTable(fit, coef=x, sort.by = "logFC", 
                                                            p.value=.05, n=nrow(countdata)))
  names(clust.tt) <- colnames(clCont)
  
  #xlsxfile <- paste0(prefix, cur.k, ".xlsx")
  #write.xlsx(clust.tt, file=paste(outdir, xlsxfile, sep="/"), row.names=T)
  write.xlsx(clust.tt, file=outfile, row.names=T)
  return(clust.tt)
}

######################################################################3
## Script 

args = commandArgs(trailingOnly=TRUE)

load(args[1]) #should have countdata
load(args[2]) #should have coldata 
if(length(args)==3) {
  outdir <- args[3]
} else if(length(args)>=4){
  clusts <- read.table(args[3], sep=",", row.names=1, header=F)
  if(min(clusts)==0){
    clusts <- clusts + 1
  }

  countdata <- countdata[,rownames(clusts)]
  coldata <- coldata[rownames(clusts),]
  coldata$cluster <- as.factor(clusts[,1])
  outdir <- args[4]
}

coldata$sample <- droplevels(coldata$sample)
coldata$tissue <- droplevels(coldata$tissue)

#DE of each cluster versus the rest
clust.de <- clustDE(coldata, countdata,
                    outfile=paste(outdir, "clusterDE.xlsx", sep="/"))

