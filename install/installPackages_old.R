
args = commandArgs(trailingOnly = TRUE )
if (length(args) >=1 ){
	if ( args[1] == "local" ) {
		cat("Installing R packages locallly\n" )
		local <-T
	} else {
		local<- F
	}
}

if (local ){
	installDir <- "./lib/R/library"
	dir.create(installDir, recursive=T )
	.libPaths(installDir)
} else {
	installDir <- .libPaths()[1]  ## the default install location
}

install.packages( c("argparser", 
		  "S4Vectors",
		  "openxlsx",
		  "Matrix",
		  "tibble"  ),
	repos='http://cran.us.r-project.org', lib=installDir )

source('https://bioconductor.org/biocLite.R')
biocLite( c("Seurat", 
                  "SummarizedExperiment",
                  "zinbwave",
                  "BiocParallel",
                  "edgeR",
		   "limma",
		   "clusterExperiment" ),  
         lib=installDir )


if ( local ) {
	


}
