
readRenviron("./Renviron")

install.packages( c("argparser", 
		  "S4Vectors",
		  "openxlsx",
		  "Matrix",
		  "tibble"  ),
	repos='http://cran.us.r-project.org' )

source('https://bioconductor.org/biocLite.R')
biocLite( c("Seurat", 
                  "SummarizedExperiment",
                  "zinbwave",
                  "BiocParallel",
                  "edgeR",
		   "limma",
		   "clusterExperiment" ) 
	)

