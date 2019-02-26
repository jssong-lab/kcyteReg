#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:37:46 2018

@author: afinneg2
"""
from __future__ import print_function
import numpy as np
import pandas as pd
from collections import OrderedDict

import matplotlib.pyplot as plt 
import seaborn as sns  


##########################################################################################
## FILTER CELLS BY SUMMED EXPRESSION 
def getSummedExpr( exprDFs ):

    summedExprDict = OrderedDict([ ])
    for clustName , df  in exprDFs.items():
         summedExprDict[clustName] =  df.sum(axis = 1)
    return summedExprDict

def filterByTotalExpr( exprDFs , minPercentile_dict):
    """
   
    """
    summedExprDict =  getSummedExpr( exprDFs )
    filterSummaryDF = pd.DataFrame(index = summedExprDict.keys(),
                                   columns = ["minSummedPerentile", "minSummedExpr", "nFiltered"] ,
                                  dtype = float)

    for key, summedExpr in summedExprDict.items():
        minPercentile =  minPercentile_dict[key]
        if minPercentile is not None:
            minSummedExpr = np.percentile( summedExpr.values ,  minPercentile)
        else:
            minSummedExpr = summedExpr.min() -1
        cells_pass = summedExpr.loc[summedExpr > minSummedExpr].index
        exprDF = exprDFs[key]
        exprDFs[key] = exprDF.loc[ cells_pass , : ].copy()
        filterSummaryDF.loc[key,  "minSummedPerentile" ] = minPercentile
        filterSummaryDF.loc[key,  "minSummedExpr"] = minSummedExpr
        filterSummaryDF.loc[key, "nFiltered"] = np.count_nonzero((summedExpr <= minSummedExpr).values)

    return exprDFs ,  filterSummaryDF

def plotSummedExpr(summedExpr_dict, figwidth = 7.5 , axHieght = 4, bins = 100,
                   minPercentile = None):

    if minPercentile  is not None:
        if isinstance(minPercentile  , float) or isinstance(minPercentile  , int):
            minPercentile_dict = { key : minPercentile for key in summedExpr_dict.keys() }
        else:
            minPercentile_dict = minPercentile

    nrows = len(summedExpr_dict)
    fig , axes = plt.subplots(nrows= nrows , ncols = 2 , figsize = (figwidth, axHieght*nrows)  )

    for rowIdx , key in enumerate(summedExpr_dict.keys()):
        ## get position of min percentile  
        vline = np.percentile( summedExpr_dict[key].values, minPercentile_dict[key] )

        ax = axes[rowIdx , 0]
        ax.hist( summedExpr_dict[key].values , bins = bins )
        ax.set_xlabel("Summed Expression")
        ax.set_ylabel("Cell Count")
        ax.set_title(key)
        ax.axvline(vline , color = "r" )

        ax = axes[rowIdx , 1]
        ax.hist( summedExpr_dict[key].values , bins = bins, cumulative= True , normed = True  )
        ax.set_xlabel("Summed Expression")
        ax.set_ylabel("Cumulative Cell Fraction")
        ax.set_title("Cumulative " + key + "\nminPercentile {:.2f} minSummedExpr {:.0f}".format(minPercentile_dict[key] ,
                                                                                              vline) )
        ax.axvline(vline , color = "r" )
        ax.axhline(minPercentile_dict[key]/100 , color = "r" )

    fig.tight_layout()
    return fig
##############################################################################################
## CALCULATE CORRELATIONS
    
def pearsonCorrel_log10tpm(exprDF):
    """ 
    exprDF - rows are cells, columns are genes
    """
    log10p1_arr = np.log10(exprDF.values*100 + 1 ) 
    cor_arr = np.corrcoef(log10p1_arr , rowvar = False )
    
    cor_df = pd.DataFrame(cor_arr , index = exprDF.columns , columns = exprDF.columns)
    return cor_df 

def pearsonCorrel_log10tpm_mcStage(exprDF, clustIDs ):
    """ 
    exprDF - rows are cells, columns are genes
    clustIDs - series, index is cells
    """
    exprDF_log10p1 = np.log10(exprDF*100.0 + 1.0 )
    exprDF_log10p1_mc = exprDF_log10p1.groupby(by = clustIDs  ).apply(lambda x : x.subtract( x.mean(axis = 0)  , axis = 1) )
    cor_arr = np.corrcoef(exprDF_log10p1_mc.values  , rowvar = False )
    cor_df = pd.DataFrame(cor_arr , index = exprDF_log10p1_mc.columns , columns = exprDF_log10p1_mc.columns)

    return  cor_df

def pearsonCorrel_log10tpm_mcStage_fromDict(exprDict):
    """ 
    exprDF - rows are cells, columns are genes
    clustIDs - series, index is cells
    """
    exprDF_log10p1_dict = OrderedDict([
                            (stageName , np.log10(df*100.0 + 1.0 )) for stageName , df in exprDict.items()
                                        ] )
    exprDF_log10p1_dict_mc = OrderedDict([
                            (stageName, df.subtract(df.mean(axis = 0), axis = 1)) for stageName, df in exprDF_log10p1_dict.items()
                                             ] )
    exprDF_log10p1_mc = pd.concat( [df for df in exprDF_log10p1_dict_mc.values()] , axis = 0 , verify_integrity=True)
    cor_arr = np.corrcoef(exprDF_log10p1_mc.values , rowvar = False )
    cor_df = pd.DataFrame(cor_arr , index = exprDF_log10p1_mc.columns , columns = exprDF_log10p1_mc.columns)

    return cor_df

#############################################################################################
## CORRELATION DISTRIBUTIONS 
    
def getMarginal_corr(corrDF, bins = 100, kwargs_hist = {}):
    
    plotData = corrDF.values.copy()[ np.triu_indices(n = corrDF.shape[0] , k =1 )   ]
    counts, bin_edges  = np.histogram( plotData , bins = bins , **kwargs_hist)
    return counts, bin_edges

def plotHistData(counts_dict , binEdges, col_wrap = 4 , FacetGrid_kws = {"aspect": 5 , "size": 3}):
    """
    counts_dict - dictionary with arrays of counts as keys. Or dictoray of dictionaries. 1st level keys are clusters. 2nd level
                dictionary values are arrays of counts
    """
    xvals = (binEdges[:-1]  + binEdges[1:] )/2
    toConcat= []
    if isinstance(counts_dict[list(counts_dict.keys())[0] ] , dict):
        ##Parse dictionary of dictionaries into a long form data frame for plotting
        for key, subDict in counts_dict.items():
            df = pd.DataFrame( subDict  )
            df.loc[: , "corr_val"] =xvals
            df_melted = pd.melt(df ,id_vars = ["corr_val"] ,
                                var_name='cor_type', value_name='count')
            df_melted.loc[:, "cluster"] =  key
            toConcat.append(df_melted)
        plotDF = pd.concat(toConcat , axis = 0)
        g = sns.FacetGrid(data = plotDF , col = "cluster" ,  col_wrap =col_wrap, hue = "cor_type", **FacetGrid_kws)
        g = g.map( plt.plot , "corr_val" , "count",  )
    else:
        df = pd.DataFrame(counts_dict )
        df.loc[: , "corr_val"] =xvals
        plotDF = pd.melt(df ,id_vars = ["corr_val"] , var_name='cor_type', value_name='count')
        g = sns.FacetGrid(data = plotDF , hue = "cor_type", **FacetGrid_kws)
        g = g.map( plt.plot , "corr_val" ,"count")
    g.add_legend()
    return g,  plotDF
    
