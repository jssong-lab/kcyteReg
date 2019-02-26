#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import argparse
import os


#### FUNCTION DEFINITIONS #######################################################
def fimoOut_to_countTable(fimoOut_df , pVal_thresh):
    
    fimoOut_df = fimoOut_df.loc[ fimoOut_df.loc[: , "p-value"] < pVal_thresh, :].copy()
    fimoOut_df["motif_label"] =  fimoOut_df.apply(lambda x: x["motif_id"] + " ," + x["motif_alt_id"]  , axis = 1)
    fimoOut_df.drop( labels = ["motif_id", "motif_alt_id" ] , axis = 1 , inplace = True )
    fimoOut_grouped = fimoOut_df.groupby( by = ["motif_label" , "sequence_name" ] )
    countTable = fimoOut_grouped.apply(lambda x: x.shape[0]).unstack(fill_value = 0)
    return countTable

def fimoOut_to_SElengths(fimoOut_df):
    SE_lengths = pd.concat( [ fimoOut_df.loc[: ,"sequence_name"] , 
                fimoOut_df.loc[: ,"sequence_name"].map(lambda x:  np.abs(int(x.split(":")[-1].split("-")[1] )  -  int(x.split(":")[-1].split("-")[0] )))
                        ] ,axis = 1)
    SE_lengths.columns = ["sequence_name", "length"]
    SE_lengths.drop_duplicates( subset = ["sequence_name"] , inplace =True )
    SE_lengths.set_index( keys = "sequence_name" , drop = True, inplace =True , verify_integrity=True )
    return  SE_lengths

### PARSE CMD LINE ##################################################
parser = argparse.ArgumentParser(
description= "Parse FIMO output to table of motif counts and spearman correlations among motifs\n\
Output files are\n\
o_counts  - table of motif counts in for each sequence\n\
o_counts.lengthNorm  - table of motif counts divided by sequence lengths\n\
o_cor - table of correlations among vectors of normalized occurence counts for  motifs")
parser.add_argument("-i" , help = "Fimo output (file with default name fime.tsv)")
parser.add_argument("--thresh" , help = "Threshold on pvalue when counting motifs", type = float , default = 0.0001)
parser.add_argument("--o_counts" )
parser.add_argument("--o_cor")
parser.add_argument( "--corType" , default = "pearson")
parser.add_argument("--motifIDLookup", help = "File with 1st column storing memeMotifIDs and 2nd column storing \
scRNAseq symbols. tab sep , - indicates Nan")
args = parser.parse_args()
#exampleArgs= "-i /home/groups/song/songlab2/afinneg2/projects/keratinocyteRegulators/analysis2/motifAnalysis/motifCounts_fimo/FIMO_threshE-4_NHEKP_TFfantomDE.txt \
#--thresh 0.0001 \
#--o_counts FIMO_threshE-4_NHEKP_TFfantomDE.motifCounts.txt \
#--o_cor FIMO_threshE-4_NHEKP_TFfantomDE.motifCor.txt"

### LOAD / PARSE DATA ##########################################################
fimoOut_df = pd.read_csv( args.i , sep = "\t" , header = 0 )
if args.motifIDLookup is not None:
    motifID_to_altID = pd.read_csv( args.motifIDLookup , header = 0 , sep = "\t" , index_col = 0 )
    motifID_to_altID = motifID_to_altID.iloc[:,0]
    fimoOut_df.loc[:, "motif_alt_id"] = fimoOut_df.loc[: , "motif_id"].map(motifID_to_altID )
    
countTable = fimoOut_to_countTable(fimoOut_df.copy() , pVal_thresh= args.thresh)

## normalize by SE length
SE_lengths = fimoOut_to_SElengths(fimoOut_df.copy())
SE_lengths = SE_lengths.loc[ countTable.columns.values , "length" ].copy()
countTable_normalized = countTable.divide( SE_lengths , axis = "columns" )

## Compute correlation matrix ##########################################################
if args.o_cor is not None:
    motifCor = countTable_normalized.transpose().corr(method =  args.corType.lower())

### WRITE OUTPUT FILES ##################################################################
if args.o_counts is not None:
    countTable.to_csv(args.o_counts , sep = "\t")
    countTable_normalized.to_csv( os.path.splitext(args.o_counts)[0] + ".lengthNorm" + os.path.splitext(args.o_counts)[1], 
                                 float_format='%.6e' , sep = "\t" )
if args.o_cor is not None:
    motifCor.to_csv( args.o_cor , float_format = "%.4f" , sep = "\t" )
    