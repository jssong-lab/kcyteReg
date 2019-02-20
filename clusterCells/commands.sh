#!/bin/bash

export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4 
export MKL_NUM_THREADS=4

#######################################################
## PCA
#fi_expr="../impute/magic_counts_t10.csv"
#fi_geneNames="../raw/genesGeq5UMIgeq100Cell_rawCount_filtered.txt"
#fi_cellNames="../raw/cellNames_foreskin.kcyte.txt"
#nPC=40
#logT_pseudo=1
#fo=PCA_impute.t10_gene.geq5UMIgeq100Cell.all.csv
#
#python ../codes/run_PCA.py --fi_expr $fi_expr --nPC $nPC --logT_pseudo $logT_pseudo --fi_cellNames $fi_cellNames --fi_geneNames $fi_geneNames --fo $fo

############################################################
## KASP clustering
#f_in="./PCA_impute.t10_gene.geq5UMIgeq100Cell.all-PCcomps.csv"
#ka=10
#k=30
#nfeat=20
#nClust=8
#alpha=4
#kmeans_nJobs=1
#N_nearest=10  ## Value does not matter becuase no predictOnly file
#fo="kasp.ka${ka}.k${k}_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
#
#mkdir -p "./runData"
#python ../codes/run_specCluster.py --f_in $f_in --fo $fo --fo_runDat "runData/$fo" --nfeat $nfeat --nClust $nClust --alpha $alpha  --kmeans_nJobs "$kmeans_nJobs"  --ka "$ka" --k "$k" --N_nearest "$N_nearest" --adaptive

##########################################################################
## Reindex clusters 
f_in="./kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
l="./clustID.reindex.lookup.txt"  ## determined by examining mean stage wise expression of marker genes
fo=${f_in}.tmp

awk -v d=, -v f=22 'BEGIN{ FS=OFS=d }
 NR==FNR {map[$1] = $2; next}
{ if ($f in map) {$f = map[$f]};  print }' $l $f_in > $fo
mv $fo $f_in

############################################################################
## Visualize in ./plot_clusterCells.ipynb



