#!/usr/bin/env Rscript 

library("argparser")
library("S4Vectors")

cat("starting \n")
p <- arg_parser(description ="Save Robject with countdata and coldata for input to AJ DE code. Like write-countdata.colData_fromGeneList.R except that instead of 
 restricting to provided geneNames take union of these geneNames with set of genes passing expressoin threshold")
p <- add_argument(p , "--fi_rawCounts" , help = "Genes are rows cells are columns comma-separated")
p <-  add_argument(p , "--cellRows" , help = "indicate that fi_rawCounts has cells in row instead of columns", flag = T)
p <- add_argument(p , "--fi_metaData" ,help = "Rows are cells. Columns have metaData")
p <- add_argument( p , "--cellNames", help = "File listing cell names. cols of rawCounts and rows of metaData are restricted to these cells"  )
p <- add_argument(p , "--geneNames" , help = "File listing geneNames to use" )
p <- add_argument( p , "--thresh_UMI" , default = 3 , type = "numeric" , help = "When selecting genes require >= thres_UMI in >= thresh_cells")
p <- add_argument(p  , "--thresh_cell" , default = 100 , type = "numeric" , help = "When selecting genes require >= thres_UMI in >= thresh_cells")
p <- add_argument(p , "--fo_countdata", help = "output name for R object")
p <- add_argument(p , "--fo_coldata", help = "output name for R object")
argv <- parse_args(p)   
                  # exampleArgs <- c("--fi_rawCounts", "../exprMats/rawCount_foreskin.csv", "--fi_metaData" , "../exprMats/metaData-cells_foreskin.csv",
                  #                 "--cellRows" , "--geneNames" , "../sets_genes/fantomDE-allGenes_sigif-posLogFC.geneNames-scRNA.geq1pct_kcyte,../sets_genes/fantomDE-TF_sigif-posLogFC.geneNames-scRNA.geq1pct_kcyte",
                  #                 "--fo_countdata", "foreskin_countData.tmp.Rdata" , "--fo_coldata" , "foreskin_colData.tmp.Rdata")
                  # argv <- parse_args(p, exampleArgs)
                  # setwd("/Users/afinneg2/projects/keratinocyteRegulators/analysis2/DE_limma")
### Load data #########################################################################################
cat("Loading data\n")
countdata <- read.csv( argv$fi_rawCounts,  row.names = 1 , header = 1  )
if (argv$cellRows){ 
  countdata <-  t(countdata) 
}
countdata_colnames <- colnames(countdata)
countdata_rownames <- rownames(countdata)
colnames(countdata) <-  gsub( "\\." , "-" , countdata_colnames  ) 
rownames(countdata) <-   gsub( "\\." , "-" , countdata_rownames ) 

coldata <- read.csv(argv$fi_metaData ,row.names = 1 , header = 1 ) 

##### Filter cells and genes ##################################################################
#################
### Filter Cells
if (! is.na(argv$cellNames  )){
  cat("Restricing countdata cells based on provided list")
  cellNames <- scan(file = argv$cellNames, what = "character")
  stopifnot(all( cellNames %in%  colnames(countdata) ))
} else{
  cellNames <-  colnames(countdata)
}
#################
## Fitler genes - First by expression and then add geneNames regardless of expression
nCellGeqThresh <- rowSums(countdata[, cellNames] >= argv$thresh_UMI )
genes_exprThreshPass <- names(nCellGeqThresh)[ nCellGeqThresh >= argv$thresh_cell]

if (! is.na(argv$geneNames)){
  cat("Taking union of expressed genes with the provided gene list (assumed to be expressed)")
  genes_include <- unique(unlist(sapply(  strsplit( argv$geneNames, "," )[[1]] , FUN = scan , what = "character"  ) ))
  stopifnot(all( genes_include%in% row.names(countdata) ))
}
genes <- union( genes_exprThreshPass, genes_include )

### Restrict countData and colData ################################################################
countdata <-  countdata[genes ,cellNames]
cat(paste("Dimesion of countdata (after filtering) are", paste( dim(countdata) ,collapse = ",") , "(genes , cells)") )
countdata <- as.matrix( countdata )

coldata <- coldata[cellNames ,  ]
coldata  <-  DataFrame(coldata)

### Save Rdata object ##########################################################################
save(countdata , file = argv$fo_countdata)
save(coldata , file = argv$fo_coldata)



