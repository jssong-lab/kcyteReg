#!/bin/bash

##################################################
## Call DE : BK vs DK                #############
#################################################

##########################
## Write cellNames_BKvsDK/txt

f_in="../clusterCells/kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv"
l="./reindex_BKvsDK.txt"
fo=cellNames_BKvsDK.txt

awk -v d=, -v f=22 'BEGIN{ FS=OFS=d }
NR==FNR {map[$1] = $2; next} 
{ if ($f in map) {$f = map[$f];  print } }' $l $f_in | cut -f1,22 -d , | sort -t "," -k 2,2n -k 1,1 > $fo

