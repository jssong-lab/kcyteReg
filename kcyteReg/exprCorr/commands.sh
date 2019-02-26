#!/bin/bash


##################################################################################
### FILTER BY SUMMED EXPRESSION 
fi_expr="../impute/magic_counts_t4_cells.foreskin_genesGeq1Pct.kcyte.pkl"
fi_stageIDs="../clusterCells/kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
minSummedPcts="1:2,2:1,3:1,4:2,5:1,6:1,7:15,8:10"
fo_cells=cellNames.passFilter.txt
fo_summary=filterSummary.pdf

#python ../codes/run_filterCellsSummedExpr.py  --fi_expr $fi_expr --fi_stageIDs $fi_stageIDs --minSummedPcts $minSummedPcts --fo_cells $fo_cells --fo_summary $fo_summary

##############################################################################
### WRITE CORRELATION FILES -- logTpm
superstages='basal:stage1,stage2,stage3;basal.stage4:stage1,stage2,stage3,stage4;differentiated:stage5,stage6,stage7;differentiated.stage4:stage4,stage5,stage6,stage7;all.noStage8:stage1,stage2,stage3,stage4,stage5,stage6,stage7'
corrMethod="logTpm"
fo=corr.$corrMethod

python ../codes/run_calcCorr.py  --fi_expr $fi_expr --fi_stageIDs $fi_stageIDs --cells $fo_cells --corrMethod $corrMethod --superstages $superstages --fo $fo
