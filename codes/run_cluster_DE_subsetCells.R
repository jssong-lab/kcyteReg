#./run_normal_DE.R coundata.RData coldata.Rdata outdir
#OR ./run_normal_DE.R coundata.RData coldata.Rdata clusters.tsv outdir 
# where clusters.tsv is a table with cell names in first column and cluster number in the second


source("/Users/afinneg2/projects/keratinocyteRegulators/analysis3/DE_limma/scripts/DE_functions.R")
source("/Users/afinneg2/projects/keratinocyteRegulators/analysis3/DE_limma/scripts/plot_functions.R")

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

