library(limma)
library(openxlsx)
library(clusterExperiment)

makeTissContrast <- function(tissCol){
  tissCol[tissCol > 0] <- 1/sum(tissCol>0)
  tissCol[tissCol == 0] <- -1/sum(tissCol==0)
  return(tissCol)
}

makePairTissContrast <- function(tissCol){
  tissCol[tissCol > 0] <- 1/sum(tissCol>0)
  tissCol[tissCol < 0] <- -1/sum(tissCol<0)
  return(tissCol)
}

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

tissClustDE <- function(coldata, countdata){
  
  coldata$sample <- droplevels(coldata$sample)
  #design <- model.matrix(~0 + sample + percent.mito + nUMI, data=coldata)
  design <- model.matrix(~0 + sample + percent.mito, data=coldata)
  
  fullCont <- rbind(apply(as.matrix(table(coldata$sample, coldata$tissue)), 2, makeTissContrast), 
                    matrix(0, nrow=dim(design)[2] - nlevels(coldata$sample), ncol=nlevels(coldata$tissue)))
  
  rownames(fullCont) <- colnames(design)
  colnames(fullCont) <- levels(coldata$tissue)
  
  
  #return(fullCont)
  #vm <- voom(countdata, design)
  #fit <- lmFit(vm, design)
  y <- new("EList")                                                                         
  y$E <- edgeR::cpm(countdata, log = TRUE, prior.count = 1)                              
  fit <- lmFit(y, design)

  fit <- contrasts.fit(fit, fullCont)
  #fit <- eBayes(fit)
  fit <- eBayes(fit, trend = T, robust=T)

  clust.tt <- lapply(colnames(fullCont), function(x) topTable(fit, coef=x, sort.by = "logFC", 
                                                              p.value=.05, n=nrow(countdata)))
  names(clust.tt) <- colnames(fullCont)
  
  return(clust.tt)

}

runTissClustDE <- function(coldata, countdata, outfile="clustTissueDE.xlsx"){

  #coldata$cluster <- as.factor(cluster)
  clustCells <- lapply(levels(coldata$cluster), function(x) rownames(subset(coldata, cluster==x)))
  clTisDE <- lapply(clustCells, function(x) tissClustDE(coldata[x,], countdata[,x]))
  clTisDE <- unlist(clTisDE, recursive=F)
  names(clTisDE) <- paste0(names(clTisDE), "_cl", rep(levels(coldata$cluster), each=nlevels(coldata$tissue)))
  cur.k <- length(levels(coldata$cluster))
  
  #xlsxfile <- paste0(prefix, cur.k, ".xlsx")
  #write.xlsx(clTisDE, file=paste(outdir, xlsxfile, sep="/"), row.names=T)
  write.xlsx(clTisDE, file=outfile, row.names=T)
  return(clTisDE)
}

diseaseDE <- function(coldata, countdata){
  
  coldata$sample <- droplevels(coldata$sample)
  #design <- model.matrix(~0 + sample + percent.mito + nUMI, data=coldata)
  design <- model.matrix(~0 + sample + percent.mito, data=coldata)
  
  tiss2samp <- as.data.frame.matrix(table(coldata$sample, coldata$tissue))
  pso.v.trunk <- makePairTissContrast((tiss2samp$psoriasis > 0) - (tiss2samp$trunk > 0))
  pso.v.all <- makeTissContrast(tiss2samp$psoriasis > 0)
  
  #compare pso v. trunk and pso v. all for each cluster
  sampCont <- cbind(pso.v.trunk, pso.v.all)
  fullCont <- rbind(sampCont, matrix(0, nrow=dim(design)[2] - nlevels(coldata$sample), ncol=2))
  rownames(fullCont) <- colnames(design)
  colnames(fullCont) <- c("pso.vs.trunk", "pso.vs.all")
  
  #return(fullCont)
  #vm <- voom(countdata, design)
  #fit <- lmFit(vm, design)
  y <- new("EList")                                                                         
  y$E <- edgeR::cpm(countdata, log = TRUE, prior.count = 1)                              
  fit <- lmFit(y, design)
  
  fit <- contrasts.fit(fit, fullCont)
  #fit <- eBayes(fit)
  fit <- eBayes(fit, trend = T, robust=T)

  clust.tt <- lapply(colnames(fullCont), function(x) topTable(fit, coef=x, sort.by = "logFC", 
                                                              p.value=.05, n=nrow(countdata)))
  names(clust.tt) <- colnames(fullCont)
  
  return(clust.tt)
  #xlsxfile <- paste0(prefix, cur.k, ".xlsx")
  #write.xlsx(clust.tt, file=paste(outdir, xlsxfile, sep="/"), row.names=T)
}

runDiseaseDE <- function(coldata, countdata, outfile="clustDiseaseDE.xlsx"){
  clustCells <- lapply(levels(coldata$cluster), function(x) rownames(subset(coldata, cluster==x)))
  clDisDE <- lapply(clustCells, function(x) diseaseDE(coldata[x,], countdata[,x]))
  clDisDE <- unlist(clDisDE, recursive=F)
  names(clDisDE) <- paste0(names(clDisDE), "_cl", rep(levels(coldata$cluster), each=2))

  #xlsxfile <- paste0(prefix, cur.k, ".xlsx")
  #write.xlsx(clDisDE, file=paste(outdir, xlsxfile, sep="/"), row.names=T)
  write.xlsx(clDisDE, file=outfile, row.names=T)
  return(clDisDE)
}