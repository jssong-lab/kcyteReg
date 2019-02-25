#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import sys
import os
from collections import OrderedDict
import argparse
from scipy.stats import mannwhitneyu
from statsmodels.stats import multitest 

from IOutils import readListFromFile

import matplotlib.pyplot as plt
import seaborn as sns  

plt.rcParams["axes.labelsize"] = "x-large"
plt.rcParams["xtick.labelsize"] =  "large"
plt.rcParams["ytick.labelsize"] =  "large"
plt.rcParams["legend.fontsize"] = "x-large"

## FUNCTION DEFINITIONS #############################################################
def get_CommonLangEffectSize(x , y):
    """
    Inputs:
        x - 1d array of observations from 1st population
        y - 1d array of observations from 2nd popuation
    """
    effectSize = (x - y[: , None] > 0).sum() / (len(x)*len(y))
    return  effectSize

def get_diffOfMedian(x , y):
    """
    Inputs:
        x - 1d array of observations from 1st population
        y - 1d array of observations from 2nd popuation
    """
    rVal = np.median(x) - np.median(y)
    return rVal

def zScore_from_Ustat(x, y):
    """
    Follows notation of wikipedia.org/wiki/Mann-Whitney_U_test
    Inputs:
        x - 1d array of observations from 1st population
        y - 1d array of observations from 2nd popuation
    Returns: zScore
    """
    n1 = len(x)
    n2 = len(y)
    n = n1 + n2
    U_x = ( (x - y[: , None] ) > 0).sum()
    U_ties =   ( ( (x - y[: , None]) == 0 ).sum() )*0.5
    U_x = U_x + U_ties
    if not np.isclose( U_ties  , 0 ): 
        _ , sharedRank_counts = np.unique( np.concatenate(( x , y), axis = 0  ), return_counts = True )
        tie_correction = np.sum(sharedRank_counts**3 - sharedRank_counts) /(n*(n-1))                             
    else:
        tie_correction = 0
    
    m_u  =  n1*n2 / 2.0
    sigma_u = np.sqrt( n1*n2*(n +1 - tie_correction ) / 12.0 )
    zScore = (U_x - m_u) / sigma_u
    return zScore
    
def get_rankBiserialCorrelation(x, y):
    U_x = ( (x - y[: , None] ) > 0).sum()
    U_y = ( (y - x[: , None] ) > 0).sum()
    r = (U_x - U_y) / ( len(x)*len(y) )
    return r
def get_bestMotifRep( motifData , motifRank_df ):
    """
    Use data frame ranking motifs corresponding to same TF to select single motif
    representative.
    Inputs
        motifData - dataframe with index of the form <motifID> ,<scRNAseq symb>
        motifRank_df DataFrame with indcies that are motifIDs and columns symbol_scRNAseq rank
    Returns:
        motifData_filtered 
    """
    col_order = ["symbols_scRNAseq" , "motifHeader"]
    cols_drop = ["reject" ]
    
    ## split the <motifID> ,<scRNAseq symb> parts of motifData index 
    motifData["symbols_scRNAseq"] =   list(map( lambda x: x.split(" ,")[-1].strip() ,  list(motifData.index) ))
    motifData.index = list(map( lambda x: x.split(" ,")[0].strip() ,  list(motifData.index) ))
    ## select best representative
    motifData.index.name = "motifHeader"
    motifData[["rank" , "origin"]] = motifRank_df.loc[motifData.index , ["rank" , "origin"]  ]
    motifData_by_symbol =  motifData.reset_index().groupby(by = "symbols_scRNAseq")
    motifData_filtered = motifData_by_symbol.apply( lambda x : x.loc[ x["rank"] == x["rank"].min(),  : ]  )  
    motifData_filtered =  motifData_filtered.reset_index(drop = True)
    #motifCounts_df_filtered.drop(labels = ["rank" , "origin" , "symbols_scRNAseq", "motifHeader"], axis = 1, inplace = True)
    cols_allowed =  col_order + [ x for x in  motifData_filtered.columns if x  not in set( cols_drop).union(col_order )] 
    motifData_filtered = motifData_filtered.loc[: ,cols_allowed  ].copy()
    
    return motifData_filtered


### GLOBALS ##############################################################################
effectMeasure_dict={"comonLangEffectSize" :  get_CommonLangEffectSize ,
                    "differenceOfMedian":  get_diffOfMedian, 
                    "zScore" : zScore_from_Ustat  , 
                    "rankBiserialCorrelation" : get_rankBiserialCorrelation }

### PARSE COMMAND-LINE ###############################################################
parser = argparse.ArgumentParser(description= "Test differential enrichment of Feature counts (e.g. motif counts) \
between two sets of observations (e.g Genomic intervals associated with two different biological conditions)" )
parser.add_argument("--f1" , help = "First table of feature counts. Rows are features (e.g motifs), columns are observations (e.g genomic intervals)" )
parser.add_argument("--f2" , help = "Second table of feature counts. Rows are features, columns (e.g motifs) are observations (e.g. genomic intervals)")
parser.add_argument("--f1_allowedObs" , help = "list of names of allowed observations")
parser.add_argument("--f2_allowedObs" , help = "list of names of allowed observations")
parser.add_argument("--motifRanks")
parser.add_argument( "--effectMeasure" , default = "differentOfMedian" , help = "comma separated list with allowed values:\n" + ", ".join(effectMeasure_dict.keys()) )
parser.add_argument( "--alpha" , default = 0.01 , type = float)
parser.add_argument("--o_allMotif" ,  help = "output filename. Positive effectSizes indicate enrichment for features in f2" )
parser.add_argument("--o_bestMotif" , help = "output file name. If multiple motifs for same TF are enriched only write the motif with highest KL divergence")
parser.add_argument("--plotEffectMeasure" , action = "store_true")

#exampleInput = "--f1 ../motifCounts_fimo/FIMO_TFsfantom-Klein-misc_bkground.combinedStages/FIMO_threshE-4_NHEKP_TFsfantom-Klein-misc.threshE-4.counts.lengthNorm.tsv \
#--f2 ../motifCounts_fimo/FIMO_TFsfantom-Klein-misc_bkground.combinedStages/FIMO_threshE-4_NHEKD_TFsfantom-Klein-misc.threshE-4.counts.lengthNorm.tsv  \
#--f1_allowedObs ./NHEKP-K27Ac-SE.unique.names \
#--f2_allowedObs ./NHEKD-K27Ac-SE.unique.names \
#--o_allMotif ./NHEKP.unique_Vs_NHEKD.unique_TFsfantom-misc.threshE-4.mannWhit.tsv \
#--effectMeasure comonLangEffectSize,differenceOfMedian,zScore,rankBiserialCorrelation \
#--alpha 0.001 \
#--plotEffectMeasure"
args = parser.parse_args()

effectMeasure_list = args.effectMeasure.split(",")
assert( all( [x in effectMeasure_dict.keys() for x in effectMeasure_list] ) )

### LOAD DATA ############################################################################
df1 = pd.read_csv(args.f1, index_col = 0 , header = 0 , sep = "\t" )
df2 =  pd.read_csv(args.f2, index_col = 0 , header = 0 , sep = "\t" )

if args.f1_allowedObs is not None:
    df1_allowedObs = readListFromFile(args.f1_allowedObs)
    df1_allowedObs_filtered = list(set( df1_allowedObs).intersection( set(df1.columns ) ))
    if len(df1_allowedObs_filtered) < len(df1_allowedObs):
        print( "{} of {} observations provided in f1_allowedObs occur in f1".format(len(df1_allowedObs_filtered) ,  len(df1_allowedObs) ) )
    df1 = df1.loc[: ,   df1_allowedObs_filtered ].copy()

if args.f2_allowedObs is not None:
    df2_allowedObs = readListFromFile(args.f2_allowedObs)
    df2_allowedObs_filtered = list(set( df2_allowedObs).intersection( set(df2.columns ) ))
    if len(df2_allowedObs_filtered) < len(df2_allowedObs):
        print( "{} of {} observations provided in f2_allowedObs occur in f2".format(len(df2_allowedObs_filtered) ,  len(df2_allowedObs) ) )
    df2 = df2.loc[: ,   df2_allowedObs_filtered ].copy()
    
##### DIFFERENTIAL ENRICHMENT TEST ########################################################
test_order = list(df1.index) 
mannWhit_dict = OrderedDict([])
for feature in test_order:
    x =  df2.loc[feature, :].values 
    y = df1.loc[feature, :].values
    U , p =  mannwhitneyu(x = x , y= y , alternative = "two-sided" )
    feature_results = [ U , p]
    for effectMeasure in effectMeasure_list:
        effectMeasure_result = effectMeasure_dict[ effectMeasure ]( x=x ,y=y  )
        feature_results.append(effectMeasure_result)
    mannWhit_dict[feature] = feature_results     
    
    mannWhit_df = pd.DataFrame.from_dict( data = mannWhit_dict, orient= "index" )
mannWhit_df.columns = ["U_mannWhitney" , "p",  ] + effectMeasure_list

## multiple test correction : Benjaminni 
reject , p_adj, _, _= multitest.multipletests( mannWhit_df.loc[:, "p"].values.copy() ,alpha = args.alpha , method = "fdr_bh" , 
                                                  is_sorted = False , returnsorted = False)
mannWhit_df["p_adj"] = p_adj 
mannWhit_df["reject"] = reject

#### FILTER TO BEST MOTIF REPRESENTATAVE ###############

if args.o_bestMotif is not None:
    ### Load the data frame of motif ranks
    motifRank_df = pd.read_csv( args.motifRanks  , index_col = 0 , header = 0 , sep = "\t" )
    ### Restrict to significant motifs 
    mannWhit_df_signif = mannWhit_df.loc[mannWhit_df.loc[: , "reject"] , : ].copy()
    mannWhit_df_bestMotif = get_bestMotifRep( motifData =mannWhit_df_signif.copy(), motifRank_df =   motifRank_df )
    mannWhit_df_bestMotif.sort_values( by = "p" , axis = 0, inplace = True )
    mannWhit_df_bestMotif.to_csv( args.o_bestMotif , sep = "\t", index = False)

### Sort by P-value and write
    
mannWhit_df.sort_values( by = "p" , axis = 0, inplace = True )
mannWhit_df.to_csv(args.o_allMotif , sep = "\t")

#### Make Plots #################################################################################
if args.plotEffectMeasure:
    ## Scatter plot P_adj vs others
    ncols = 2
    nrows = -1*((-len(effectMeasure_list))//ncols )
    fig , axes = plt.subplots(nrows = nrows , ncols = ncols, figsize = (7.5 , 4*nrows))
    axes = np.ravel(axes)
    for ax , effectMeasure in zip(axes , effectMeasure_list):
        ax.scatter( x = -np.log10(mannWhit_df.loc[: , "p_adj"].values),
                      y = mannWhit_df.loc[: , effectMeasure].values )
        ax.set_xlabel("p_adj")
        ax.set_ylabel(effectMeasure)
        if effectMeasure == "differenceOfMedian":
            ax.set_ylim(-0.004 , 0.004)
    fig.tight_layout()
    fig.savefig( os.path.splitext( args.o_allMotif )[0]  + "pAdj_v_effectSize.pdf" , format = "pdf" )

