#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 20:57:24 2018

@author: afinneg2
"""

from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import os
from collections import OrderedDict
import argparse

import corrFuncs as cf
import IOutils

import matplotlib.pyplot as plt
plt.rcParams["axes.labelsize"] = "x-large"
plt.rcParams["xtick.labelsize"] =  "large"
plt.rcParams["ytick.labelsize"] =  "large"
plt.rcParams["legend.fontsize"] = "x-large"

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42

############################################################################################################
## PARSE COMMAND LINE
parser = argparse.ArgumentParser(description= "Filter Cells by summed expression" )
parser.add_argument("--fi_expr" , help = "table of expression values (csv or pkl format) . Rows are cells. Columns are genes")
parser.add_argument("--fi_stageIDs" , help = "table of series (csv or pkl format) with cells as index and column clust_ID")
parser.add_argument( "--minSummedPcts" , help= "list of form <Name1>:<float0-100>, ... . \
float is the lower percentile to trim. 0 => no trim." )
parser.add_argument( "--fo_cells")
parser.add_argument("--fo_summary")
args = parser.parse_args()

#exampleArgs="--fi_expr ../exprData/impute-t4_libNorm10K_cells.foreskin_genesGeq1Pct.all.kcyte.pkl \
#--fi_stageIDs ../clusterCells/kasp.ka10.k30_impute.t4_gene.geq5UMIgeq100Cell.all.csv \
#--minSummedPcts 1:2,2:1,3:1,4:2,5:1,6:1,7:15,8:10 \
#--fo_cells cellNames.passFilter.tmp \
#--fo_summary filterSummary_tmp"
#args = parser.parse_args(exampleArgs.split())

minSummedPct_dict = IOutils.parseStrToDict(args.minSummedPcts, valueType= "float", pairSep=',')
##########################################################################################################
## LOAD DATA 
cellData = IOutils.loadCellData(OrderedDict([ ("expr", args.fi_expr) , ("pcComps" , args.fi_stageIDs) ]) )
cellData = pd.concat([cellData.loc[ : ,  ("expr" , slice(None))] ,  cellData.loc[ : ,  ("rowData" , "clust_ID")] ], 
                     axis = 1 )  ## Remove unnecessary columns

##########################################################################################################
## SPLIT CELLS BY STAGE ID, ENFORCE FILTERS, write results
#### Get type
clustID_dtype = cellData["rowData"].loc[: , "clust_ID"].dtype
if clustID_dtype.name.startswith("int"):
    minSummedPct_dict = OrderedDict([ (int(k), v) for k, v in  minSummedPct_dict.items()])
elif clustID_dtype.name.startswith("float"):
    minSummedPct_dict = OrderedDict([ (float(k), v) for k, v in  minSummedPct_dict.items()])
elif clustID_dtype.name.startswith("object"):
    pass
else:
    raise Exception( "Type of \"clust_ID\" column not recognized")
        
### Split 
exprDict = OrderedDict([
            ("stage {}".format(str(i)), cellData["expr"].loc[ cellData["rowData"].loc[: ,"clust_ID"] == i , : ].copy() ) \
                        for i in  minSummedPct_dict.keys() ])
## Filter
print("Filtering")
summedExpr_dict = cf.getSummedExpr(exprDict )
fig = cf.plotSummedExpr( summedExpr_dict,
                        minPercentile = OrderedDict([ ("stage {}".format(str(k)), v ) \
                                                             for k,v in  minSummedPct_dict.items()  ]) )
exprDFs ,  filterSummaryDF = cf.filterByTotalExpr( exprDFs = exprDict ,
                                    minPercentile_dict = OrderedDict([ ("stage {}".format(str(k)), v ) \
                                                             for k,v in  minSummedPct_dict.items()  ])  )
## Write
print("saving Results")
fig.savefig(os.path.splitext(args.fo_summary)[0] + ".pdf" , dpi = 150 , bbox_inches = "tight" )
filterSummaryDF.to_csv( os.path.splitext(args.fo_summary)[0] + ".tsv" , sep = "\t")
cellsUse = []
for exprDF in exprDFs.values():
    cellsUse.extend(list(exprDF.index ) )
cellOrder = list(cellData.index)
cellsUse = sorted( cellsUse, key = lambda x: cellOrder.index(x))
f = open(args.fo_cells , 'w' )
for x in cellsUse:
    f.write(x+ "\n")
f.close()
print("Done")


