#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 21:00:04 2018

@author: afinneg2
"""
from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import os
from collections import OrderedDict
import re
import argparse

import corrFuncs as cf
import IOutils


############################################################################################################
## PARSE COMMAND LINE
parser = argparse.ArgumentParser(description= "Calculate gene-gene correlation within stages and superstages")
parser.add_argument("--fi_expr" , help = "table of expression values (csv or pkl format). Corrlations using logTPM assumes \
that inpute expression is in units per 10000. Rows are cells. Columns are genes")
parser.add_argument("--fi_stageIDs" , help = "table of series (csv or pkl format) with cells as index and column clust_ID")
parser.add_argument( "--cells" , help= "File Listing allowed cells" )
parser.add_argument("--corrMethod" , help = "one of logTpm or logTpmMC", default = "logTpm")
parser.add_argument("--superstages")
parser.add_argument("--noSingleStageCorr" , action = "store_true")
parser.add_argument( "--fo")
args = parser.parse_args()

#exampleArgs="--fi_expr ../exprData/impute-t10_libNorm10K_cells.foreskin_genesGeq1Pct.all.kcyte.pkl \
#--fi_stageIDs ../clusterCells/kasp.ka10.k30_impute.t10_gene.geq5UMIgeq100Cell.all.csv \
#--cells ./cellNames.passFilter.tmp \
#--superstages progenitor:stage1,stage2,stage3;\
#progenitor.stage4:stage1,stage2,stage3,stage4;\
#differentiated:stage5,stage6,stage7;\
#differentiated.stage4:stage4,stage5,stage6,stage7;\
#all.noStage8:stage1,stage2,stage3,stage4,stage5,stage6,stage7 \
#--corrMethod  logTpm \
#--fo corr_tmpA"
#args = parser.parse_args(exampleArgs.split())

if args.superstages is not None:
    superstage_dict = IOutils.parseStrToDict(args.superstages , valueType= "str", pairSep=";")
    try:
        superstage_dict = OrderedDict([(k , [ int(x) for x in v.split(",")] ) for k , v in superstage_dict.items() ])
    except:
        superstage_dict = OrderedDict([(k ,[int(re.search(r'(\d+)$', x).group(1)) for x in v.split(",")]) for k , v in superstage_dict.items() ])
        
############################################################################################################
## Load data 
cellData = IOutils.loadCellData(OrderedDict([ ("expr", args.fi_expr) , ("pcComps" , args.fi_stageIDs) ]) )
cellData = pd.concat([cellData.loc[ : ,  ("expr" , slice(None))] ,  cellData.loc[ : ,  ("rowData" , "clust_ID")] ], 
                     axis = 1 )  ## Remove unnecessary columns
cells_allowed = IOutils.readListFromFile(args.cells)
print("Restricting to {:d} allowed cells".format(len(cells_allowed)))
cellData = cellData.loc[cells_allowed, :].copy()

############################################################################################################
## COMPUTE CORRELATION BY STAGE
if not args.noSingleStageCorr:
    stages = sorted( cellData["rowData"].loc[: , "clust_ID"].unique() , key = lambda x: int(x) )
    corrDict_stages = OrderedDict([])
    for stage in stages:
        print("\tCalculating correlation for stage {}".format(stage))
        exprDF = cellData["expr"].loc[ cellData["rowData"].loc[: ,"clust_ID"] == stage , : ].copy() 
        corrDict_stages[stage] = cf.pearsonCorrel_log10tpm( exprDF )    
    ##Write results
    for stage , corrDF in  corrDict_stages.items():
        fo = os.path.splitext(args.fo)[0] + "_stage{:d}.pkl".format(stage) 
        print("\tWriting {}".format(fo))
        corrDF.to_pickle( fo)
        
if args.superstages is not None:
    corrDict_superstages = OrderedDict([])
    for superstage , substageList in superstage_dict.items():
        print("\tCalculating correlation for superstage {} with substages {}".format(superstage, substageList))
        if args.corrMethod == "logTpmMC":
            exprDict = OrderedDict([ 
                        (i , cellData["expr"].loc[ cellData["rowData"].loc[: ,"clust_ID"] == i , : ].copy() )  
                                for i in substageList ] )
            corrDict_superstages[superstage] = cf.pearsonCorrel_log10tpm_mcStage_fromDict(exprDict = exprDict )
        elif args.corrMethod == "logTpm":
            exprDF = cellData["expr"].loc[ cellData["rowData"].loc[: ,"clust_ID"].map( lambda x : x in substageList) , : ].copy()
            corrDict_superstages[superstage] = cf.pearsonCorrel_log10tpm(exprDF)
        else:
            raise Exception("corrMethod {} not recognized".format(args.corrMethod))   
    ##Write results
    for superstage , corrDF in  corrDict_superstages.items():
        fo =os.path.splitext(args.fo)[0] + "_stage{}.pkl".format(superstage)
        print("\tWriting {}".format(fo))
        corrDF.to_pickle(fo)
print("Done")

