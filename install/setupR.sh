#!/bin/bash

## Usage ./install/setupR.sh <local>

projDir=$(pwd)

R_MAX_NUM_DLLS=500
echo "R_MAX_NUM_DLLS=500" > .Renviron
if [ "$1" == "local" ]
	then
	echo installing locally and setting R_LIBs_USER in .Renviron
	echo "R_LIBS_USER=$projDir/lib/R/library" >> .Renviron
	mkdir -p $projDir/lib/R/library
	Rscript -e "readRenviron('.Renviron'); \
install.packages( c('argparser', 'S4Vectors', 'openxlsx', 'Matrix', 'tibble'), repos='http://cran.us.r-project.org', lib=Sys.getenv('R_LIBS_USER') ); \
source('https://bioconductor.org/biocLite.R'); \
biocLite( c('Seurat', 'SummarizedExperiment','zinbwave', 'BiocParallel', 'edgeR', 'limma','clusterExperiment' ), lib=Sys.getenv('R_LIBS_USER')  )"

else
	echo "installing to default R library"
	Rscript -e "install.packages( c('argparser', 'S4Vectors', 'openxlsx', 'Matrix', 'tibble'), repos='http://cran.us.r-project.org' ); \
source('https://bioconductor.org/biocLite.R'); \
biocLite( c('Seurat', 'SummarizedExperiment','zinbwave', 'BiocParallel', 'edgeR', 'limma','clusterExperiment' ) )"
fi

## Insure all call to  R code uses definitions of .Renvrion
for d in $(find kcyteReg -maxdepth 5 -type d)
 	do 
	if [ -e $d/commands.sh  ]
		then    
		if grep -q "\.R"  $d/commands.sh     
			then
			printf "readRenviron(\"%s\")\n.libPaths(Sys.getenv(\"R_LIBS_USER\"))\nR_MAX_NUM_DLLS=as.integer(Sys.getenv(\"R_MAX_NUM_DLLS\"))\n" $projDir/.Renviron > $d/.Rprofile 
		fi
	fi
done


