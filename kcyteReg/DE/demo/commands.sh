#!/bin/bash
mkdir -p ./oe
set -euo pipefail


##################################################
## Call DE : BK vs DK                #############
#################################################

#######################################
## Write cellNames_BKvsDK.txt
f_in="../../clusterCells/kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
l="../reindex_BKvsDK.txt"
fo=cellNames_BKvsDK.txt

awk -v d=, -v f=22 'BEGIN{ FS=OFS=d }
NR==FNR {map[$1] = $2; next} 
{ if ($f in map) {$f = map[$f];  print } }' $l $f_in | cut -f1,22 -d , | sort -t "," -k 2,2n -k 1,1 > $fo

#######################################
## Write Rdata files for DE
fi_rawCounts="../../raw/rawCount_filtered_foreskin.kcyte.csv"
fi_metaData="../../raw/metaData-cells_filtered.csv"
cellNames="../../raw/cellNames_foreskin.kcyte.txt"
geneNames="../../setsGenes/kcyte_genes.txt"
thresh_UMI=3
thresh_cell=20   ## this is calculated using x/100 = 21933/92889 => x \approx 23.6 
		 ##  and  x/100 = 21274/92889 => x \approx 22.9  (then round to 20).  
		 ## Where 92889={numbe cells in full data set} ,
		 ##  21933= {number of foreskin keratinocyte excluding channel stage}
		 ## and where  100 is number of expressing cells required by in prev. paper on full data 
fo_countData="DE_countData.Rdata"
fo_colData="DE_colData.Rdata"
name="write_Rdata"

Rscript ../../codes/write-countdata.colData_includeGeneList.R --fi_rawCounts $fi_rawCounts --cellRows --fi_metaData $fi_metaData --cellNames $cellNames --geneNames $geneNames --thresh_UMI $thresh_UMI --thresh_cell $thresh_cell  --fo_countdata $fo_countData --fo_coldata $fo_colData 1> ./oe/${name}.o 2> ./oe/${name}.e

#######################################
## Run DE 
fi_countData="DE_countData.Rdata"
fi_colData="DE_colData.Rdata"
fi_cluster="./../cellNames_BKvsDK.txt"
outDir="BK.vs.DK"
mkdir -p $outDir
name="run_DE_BKvPK"

Rscript ../../codes/run_cluster_DE_subsetCells.R  $fi_countData $fi_colData $fi_cluster $outDir 1> ./oe/$name.o 2> ./oe/$name.e  

####################################3
## Convert to csv
python -c '
import pandas as pd
sheetNames=["BK.vs.DK" ,  "DK.vs.BK" ]
x = pd.read_excel("./BK.vs.DK/clusterDE.xlsx", sheet_name= None)
for name, v in zip(sheetNames, x.values() ):
   v.to_csv( "./BK.vs.DK/" + name + ".csv"  )
'

##################################################
## Call DE : One vs Rest             #############
#################################################

#######################################
## Write cellNames.noStage8.txt
f_in="../../clusterCells/kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
fo="cellNames.noStage8.txt"
tail -n+2 $f_in | cut -f1,22 -d , | awk -F "," '{ if ( $2 != 8 ) { print $0 } }' | sort -t "," -k 2,2n -k 1,1  > $fo 

#######################################
## Run DE 
fi_countData="DE_countData.Rdata"
fi_colData="DE_colData.Rdata"
fi_cluster="./../cellNames.noStage8.txt"
outDir="one.vs.rest"
mkdir -p $outDir
name="run_DE_oneVrest"

Rscript ../../codes/run_cluster_DE_subsetCells.R  $fi_countData $fi_colData $fi_cluster $outDir 1> ./oe/$name.o 2> ./oe/$name.e

####################################3
## Convert to csv
python -c '
import pandas as pd
sheetNames=["stage1" ,  "stage2" , "stage3" , "stage4" , "stage5" , "stage6" , "stage7" ]
x = pd.read_excel("./one.vs.rest/clusterDE.xlsx", sheet_name= None)
for name, v in zip(sheetNames, x.values() ):
   v.to_csv( "./one.vs.rest/" + name + ".csv"  )
'

####################################################################################
### Write kcyteGenes not DE up in BK  and not DE up in DK (for network analysis) ####
####################################################################################

## Get genes DE up in each condition
lfc_thresh=0.25
pAdj_thresh=0.05

f_in="./BK.vs.DK/DK.vs.BK.csv"
f_o="./BK.vs.DK/DK_de_lfc.geq0.25.txt"
tail -n+2 $f_in | awk -F"," -v c=$lfc_thresh -v p=$pAdj_thresh '{ if ( ( $2 + 0.0 ) > c &&  ( $6 +0.0 ) < p ) { print $1 } }' > $f_o

f_in="./BK.vs.DK/BK.vs.DK.csv"
f_o="./BK.vs.DK/BK_de_lfc.geq0.25.txt"
tail -n+2 $f_in | awk -F , -v c=$lfc_thresh -v p=$pAdj_thresh '{ if ( ( $2 + 0.0 ) > c &&  ( $6 + 0.0 ) < p ) { print $1 } }' > $f_o

### Subract DE sets from kcyte genes
mkdir -p ./BK.vs.DK/regNetwork_targets/
f_o="./BK.vs.DK/regNetwork_targets/kcyteGenes_subtract_DEup.DK.txt"
comm -23 <( cat ../setsGenes/kcyte_genes.txt | sort  ) <( cat ./BK.vs.DK/DK_de_lfc.geq0.25.txt | sort ) | sort | uniq > $f_o 

f_o="./BK.vs.DK/regNetwork_targets/kcyteGenes_subtract_DEup.BK.txt"
comm -23 <( cat ../setsGenes/kcyte_genes.txt | sort  ) <( cat ./BK.vs.DK/BK_de_lfc.geq0.25.txt | sort ) | sort | uniq > $f_o  


