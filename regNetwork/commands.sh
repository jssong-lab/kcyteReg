#!/bin/bash
set -euo pipefail

##################################################################
## Write TFs with binding motifs enriched in BK or DK SE's
##  + flatten motifs consisiting of TF dimers

fi_motifs_arr=("../raw/motifsEnriched_SEs_BK.txt" "../raw/motifsEnriched_SEs_DK.txt" )
fo_arr=( "TFsEnriched_SEs_BK.expressed.txt"  "TFsEnriched_SEs_DK.expressed.txt" )
genes_expressed="../setsGenes/kcyte_TFs.txt"

for (( i=0; i < ${#fi_motifs_arr[@]}; i++  ))
	do
cat ${fi_motifs_arr[$i]} | python -c '
import sys
with open(sys.argv[1], "r" ) as f:
        genesAllowed = [x.strip("\n") for x in f.readlines() ]
for line in sys.stdin:
        line_list = line.strip("\n").split(",")
        if all([ x in genesAllowed  for x in line_list ] ):
                for x in line_list:
                        sys.stdout.write(x+ "\n")
' $genes_expressed | sort | uniq > ${fo_arr[$i]}
done

###################################################################
## Cluster gene and TF Modules
geneSet_arr=( "../DE/BK.vs.DK/regNetwork_targets/kcyteGenes_subtract_DEup.DK.txt" \
"../DE/BK.vs.DK/regNetwork_targets/kcyteGenes_subtract_DEup.BK.txt"    )
TF_arr=( "./TFsEnriched_SEs_BK.expressed.txt" "./TFsEnriched_SEs_DK.expressed.txt"  )
geneCorr_arr=( "../exprCorr/corr_stagebasal.stage4.pkl" \
"../exprCorr/corr_stagedifferentiated.stage4.pkl"  ) 
outDir_arr=( "correlNetwork_BK" "correlNetwork_DK" )
beta=4
linkage_regTargets="average"
linkage_TFs="average"

for (( i=0; i<${#geneSet_arr[@]}; i++ ))
	do
	mkdir -p ${outDir_arr[i]}
	python ../codes/run_clusterCoRegTFcorr.py --geneCorr ${geneCorr_arr[i]} --genes ${geneSet_arr[i]} --TFs ${TF_arr[i]} --beta $beta --fo ${outDir_arr[i]}/clustMap.pdf --linkage_regTargets $linkage_regTargets --linkage_TFs $linkage_TFs
done

#####################################################################
## Call Flat clusters - BK network
clustData="./correlNetwork_BK/clustMap.pkl"
depth_rows=4
thresh_rows=2.15
depth_cols=2
thresh_cols=0.75
outDir="correlNetwork_BK/"

python ../codes/run_getClustMapBlocks_inconsistStat.py --clustMap $clustData --depth_rows $depth_rows --thresh_rows $thresh_rows --depth_cols $depth_cols --thresh_cols $thresh_cols --outDir $outDir

#########################################################################
## Call Flat clusters - DK network
clustData="./correlNetwork_DK/clustMap.pkl"
depth_rows=16
thresh_rows=3.2
depth_cols=2
thresh_cols=0.75
outDir="correlNetwork_DK/"

python ../codes/run_getClustMapBlocks_inconsistStat.py --clustMap $clustData --depth_rows $depth_rows --thresh_rows $thresh_rows --depth_cols $depth_cols --thresh_cols $thresh_cols --outDir $outDir

##########################################################################
## Threshold on magnitude of average correlations between TF and gene pairs 
## see ./draw_regGraph.ipynb




