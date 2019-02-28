#!/bin/bash

export nThreads=4
export nWorkers=4   ## run with 18 cores

export OPENBLAS_NUM_THREADS=$nThreads
export OMP_NUM_THREADS=$nThreads
export MKL_NUM_THREADS=$nThreads
export NSLOTS=$nWorkers
##############################
## Run zinbwave 
mkdir -p zinb
Rscript ../codes/run_zinb.R ../raw/rawCount_filtered.csv ../raw/metaData-cells_filtered.csv zinb

##############################
## Run MAGIC with zinb-wave reduced dimensional representation
#
fi_expr="../raw/rawCount_filtered.csv"      #csv file with raw counts (genes are rows, columns are cells )
fi_distData="./zinb/zinb_W.csv"               ## csv file with coordinates giving cel-cell distances (rows are cells, columns are coordinates (e.g PC coords for rows of zinb-wave matrix)  )
fo_imputed="magic_counts_t{}.csv"        ## name for output files {} is filled with corresponding value of MAGIC t parameter
fo_details="runData/impute.runData"          # information on run including runtimes
minFracExpr=0.001                            ## fraction of cells for which raw UMI > 0 needed to preform imputation
tList=4,10                                     ## comma separated list of t values for which imputation is performed
k=30
ka=10

mkdir -p runData
python ../codes/run_MAGIC.py --fi_expr "$fi_expr" --fi_distData "$fi_distData" --fo_imputed "$fo_imputed"  --fo_details "$fo_details" --geneRows --minFracExpr $minFracExpr --libNormalize --tList $tList -k "$k" --ka "$ka"

##################################
## Standardize magic_counts to 10K lib sizes
standLibSize=10000
for f_i in magic_counts_t[1-9]*.csv
	do
	fo=$f_i
	python ../codes/libNorm.py --f_i  "$f_i" --standLibSize "$standLibSize" --fo "$fo"
done

##################################
## Restrict to genes and cells of interest and serialize for fast loading

for f_i in magic_counts_t[1-9]*.csv
	do
	fo=${f_i%.csv}_cells.foreskin_genesGeq1Pct.kcyte.pkl
	python ./serializeTable.py --fi $f_i --rowNames ../raw/cellNames_foreskin.kcyte.txt --colNames ../raw/genesGeq1pct_allTissue.kcyte.txt  --fo $fo
done



