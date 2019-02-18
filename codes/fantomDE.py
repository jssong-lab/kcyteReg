#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import numpy as np
import pandas as pd
from scipy.stats import  ttest_ind ,  mannwhitneyu
from collections import OrderedDict
from  statsmodels.stats.multitest import multipletests


########################################################################################################
#### Utilities
def calc_foldChange( group1_df , group2_df ):
    """
    calculate log2( <expression_group1> / <expression_group2>)  where <> indicates mean
    Inputs:
         group1_df  - (nObservations , nResponse variables) pandas.DataFrame
         group2_df -  (nObservations , nResponse variables) pandas.DataFrame
    """
   
    assert group1_df.shape[1] == group2_df.shape[1]
    logFC =  group1_df.mean(axis = 0).map( np.log2 ) - group2_df.mean(axis = 0).map( np.log2 )
    return logFC

def calc_foldChange_median( group1_df , group2_df ):
    """
    calculate log2( <expression_group1> / <expression_group2>)  where <> indicates median
    Inputs:
         group1_df  - (nObservations , nResponse variables) pandas.DataFrame
         group2_df -  (nObservations , nResponse variables) pandas.DataFrame
    """
   
    assert group1_df.shape[1] == group2_df.shape[1]
    logFC =  group1_df.median(axis = 0).map( np.log2 ) - group2_df.median(axis = 0).map( np.log2 )
    return logFC
def run_ttestInd( group1_df , group2_df , progressUpdate =200 , ttest_ind_kwargs = {},
                 pAdj_method =  'bonferroni' ,  alpha = 0.05):
    """
    Run scipy.stats.ttest_ind test between corresponding columns of group1_df and group2_df.
    The returned t_stat is for group1 cloumn vs group2 column
    Inputs:
         group1_df  - (nObservations , nResponse variables) pandas.DataFrame
         group2_df -  (nObservations , nResponse variables) pandas.DataFrame
    returns:
        results_df - (nResponse variables , 2) data frame where columns store t_stat and p_vals
    """
    assert group1_df.shape[1] == group2_df.shape[1]
    nTests = group1_df.shape[1]
    results_df = pd.DataFrame(-1*np.ones(shape =(nTests, 2), dtype= float ) ,
                                index = group1_df.columns ,
                                columns = ["t_stat" , "p_val" ] )
    testIdx = 0
    for responseVar in group1_df.columns:
        t_stat , p_val =  ttest_ind(  group1_df.loc[: , responseVar].values ,
                                    group2_df.loc[: , responseVar].values ,
                                   **ttest_ind_kwargs )
        results_df.loc[responseVar, :] = t_stat , p_val
        testIdx+=1
        if testIdx % progressUpdate ==0:
            print("completed {} of {} tests".format(testIdx , nTests))
    if pAdj_method:
        _ ,p_adj, _ ,_ = multipletests(pvals = results_df.loc[: ,"p_val" ].values ,alpha =   alpha  ,method =  pAdj_method)
        results_df['p_adj'] = p_adj
    return results_df

def run_mannWhitneyU( group1_df , group2_df , progressUpdate =200 ,
                     alternative = "greater" ,
                     pAdj_method =  'bonferroni',
                     alpha = 0.05):
    """
    Run scipy.stats.mannwhitneyu  test between corresponding columns of group1_df and group2_df.
    The returned u_stat is for group1 cloumn vs group2 column
    Inputs:
        group1_df  - (nObservations , nResponse variables) pandas.DataFrame
        group2_df -  (nObservations , nResponse variables) pandas.DataFrame
        alternative - the alternative kwarg to scipy.stats.mannwhitneyu. "greater" means
                     test with alternative hypothesis that group1 > group2
        alpha - alpha kwarg to multipletests (see AFutils.stats.multitest)
    returns:
        results_df - (nResponse variables , 2) data frame where columns store t_stat and p_vals
    """
    assert group1_df.shape[1] == group2_df.shape[1]
    nTests = group1_df.shape[1]
    results_df = pd.DataFrame(-1*np.ones(shape =(nTests, 2), dtype= float ) ,
                                index = group1_df.columns ,
                                columns = ["u_stat" , "p_val" ] )
    testIdx = 0
    for responseVar in group1_df.columns:
        u_stat , p_val =  mannwhitneyu(  group1_df.loc[: , responseVar].values ,
                                           group2_df.loc[: , responseVar].values ,
                                           alternative = alternative )
        results_df.loc[responseVar, :] = u_stat , p_val
        testIdx+=1
        if testIdx % progressUpdate ==0:
            print("completed {} of {} tests".format(testIdx , nTests))
    if pAdj_method:
        _ ,p_adj, _ ,_ = multipletests(pvals = results_df.loc[: ,"p_val" ].values ,alpha = alpha, method =  pAdj_method)
        results_df['p_adj'] = p_adj
    return results_df

def sort_absCol(df, by , ascending, absMask ,inplace= True   ):
    """
    sort data frame by column value, allowing for sorting by abs(colValue)
    Inputs:
        df - dataFrame
        by - list of column names
        ascending - list of bools 
        absMask - bools indicating whether abs is applied
    Returns
        df
    """
    absCols = [ col for col, maskElem in zip( by , absMask ) if maskElem  ]
    for colName in absCols:
        df[colName + "_abs"] = df.loc[: , colName].map( np.abs )
    sortCols = [  col if not maskElem else col + "_abs" for col, maskElem in zip( by , absMask )  ]
    df.sort_values(by = sortCols , ascending = ascending , inplace= inplace)
    df.drop(labels = [colName + "_abs" for colName in absCols] , axis = 1 , inplace = True)
    return df

def sort_signifFirst(df, threshCol , thresh , by , ascending, absMask ):
    """
    
    """
    signifMask = df.loc[: , threshCol] < thresh
    df_signif = df.loc[signifMask , :].copy()
    df_nonSignif = df.loc[ ~ signifMask , :].copy()

    df_signif = sort_absCol(df_signif, by= by ,
                            ascending = ascending,
                            absMask= absMask ,
                            inplace= True   )
    df_nonSignif = sort_absCol( df_nonSignif ,
                                by = [threshCol] + by,
                                ascending = [True] + ascending,
                                absMask = [False] + absMask,
                                inplace= True)
    df_signifFirst = pd.concat( [df_signif  ,  df_nonSignif ],
                               axis = 0)
    return df_signifFirst

def calc_avgExpr(data_df, groupLabels_series ):
    """
    Inputs:
        data_df - (nObservations , nRepsonseVars) dataFrame
        groupLabels - (nObservations,) pandas.Series of group indices for each observation
    Returns:
        results_df - (nRepsoneVars , nGroups) each entry is mean of responcse variable for 
                    observations in the group
    """
    data_df_grouped = data_df.groupby(groupLabels_series)
    groups = np.sort(list(data_df_grouped.groups.keys() ) )
    resultsDict = OrderedDict([])
    for groupID in groups:
        avgExpr =  data_df_grouped.get_group(groupID).mean(axis = 0)
        resultsDict["{}_avgExpr".format(groupID)] = avgExpr
    results_df = pd.DataFrame(resultsDict )
    return results_df

########################################################################################################
#### Workflow functions of calling DE genes between two groups
   
def callDE_genes( group1_df , group2_df ,
                    group1Name = "group1" , group2Name = "group2",
                    signifThresh = 0.05,
                    testName = "ttestInd",
                    testKwargs = {'progressUpdate': 1000,
                                "pAdj_method" : "fdr_bh"} ,
                    minGroup1_expr_mean = None , minGroup2_expr_mean = None,
                    minGroup1_expr_median = None , minGroup2_expr_median = None,
                    minLog2FC_mean = None , minLog2FC_median = None,
                    sortByLogFC_median = False):
    """
    A wrapper performing the follows steps for calling DE genes in group1_df relative
    to group1:
        1 remove genes with 
            mean(group1_expr) < minGroup1_expr  or 
            mean(group2_expr) < minGroup2_expr
        2. remove genes with:
            log2FC <  minLog2FC_magnitude
        3. Perform DE test on each of the remaining genes
    Inputs:
        group1_df : (nObservations , nResponseVars) pandas dataframe
        group2_df : (nObservations , nResponseVars) pandas dataframe
        signifThresh :  threshold for sorting genes by p_adj
        ...
    Returns: 
        results_df - rows are nResponseVars columns are: 
                stat, p_val ,p_adj , meanExpr_group1 , meanExprGroup2, logFC
    """
    ## calacualte average/median expression within groups and logFC between groups
    group1_avgExpr = group1_df.mean(axis = 0)
    group1_avgExpr.name = "avgExpr_{}".format(group1Name)
    group1_medExpr = group1_df.median(axis = 0)
    group1_medExpr.name = "medianExpr_{}".format(group1Name)
    group2_avgExpr = group2_df.mean(axis = 0)
    group2_avgExpr.name = "avgExpr_{}".format(group2Name )
    group2_medExpr = group2_df.median(axis = 0)
    group2_medExpr.name = "medianExpr_{}".format(group2Name )
    foldChangeResult = calc_foldChange(group1_df, group2_df )
    foldChangeResult.name = "logFC"
    foldChangeResult_median = calc_foldChange_median( group1_df , group2_df )
    foldChangeResult_median.name = "logFC_median"
    ## mask to select genes passing filtering criteria
    mask  = pd.Series( np.ones_like(group1_df.columns.values, dtype = bool), index = group1_df.columns.values)
    if minGroup1_expr_mean is not None:
        print("filtering on group1_expr_mean with minimum {}".format(minGroup1_expr_mean)  )
        mask =  mask  & ( group1_avgExpr > minGroup1_expr_mean)
    if minGroup2_expr_mean is not None:
        print("filtering on group2_expr_mean with minimum {}".format(minGroup2_expr_mean)  )
        mask =  mask  & ( group2_avgExpr > minGroup2_expr_mean )
    if minGroup1_expr_median is not None:
        print("filtering on group1_expr_median with minium {}".format(minGroup1_expr_median))
        mask = mask  & ( group1_medExpr > minGroup1_expr_median )
    if minGroup2_expr_median is not None:
        print("filtering on group2_expr_median with minium {}".format(minGroup2_expr_median))
        mask = mask & (  group2_medExpr > minGroup2_expr_median)
    if minLog2FC_mean  is not None:
        print("filtering on  minLog2FC_mean with minimum {}".format( minLog2FC_mean ))
        mask = mask & (foldChangeResult >   minLog2FC_mean )
    if minLog2FC_median  is not None:
        print("filtering on minLog2FC_median with minimum {}".format( minLog2FC_median ) )
        mask = mask & ( foldChangeResult_median >  minLog2FC_median   )
    print("Excluding {:d} genes from DE test".format(np.count_nonzero(~mask )))
    print("Running DE test on {:d} genes".format(np.count_nonzero(mask)) )

    ## run DE test
    allowedTests = {'ttestInd' : run_ttestInd ,
                   'mannWhitneyU' : run_mannWhitneyU}
    if 'alpha' not in  testKwargs.keys():
        testKwargs['alpha'] = signifThresh
    testFunc = allowedTests[testName]
    print("comparing groups with test {}".format(testName))
    test_results =  testFunc( group1_df.loc[:, mask].copy() ,
                            group2_df.loc[: , mask].copy() ,
                             **testKwargs )
    results_df = pd.concat( [test_results,
                             group1_avgExpr, group2_avgExpr , foldChangeResult ,
                             group1_medExpr , group2_medExpr ,foldChangeResult_median ] ,
                            axis =1 , verify_integrity=True  )
    if sortByLogFC_median:
        results_df = sort_signifFirst(results_df , threshCol = "p_adj" , thresh = signifThresh,
                                         by =['logFC_median'] ,
                                         ascending = [False] ,
                                         absMask = [True] )
    else:
        results_df = sort_signifFirst(results_df , threshCol = "p_adj" , thresh = signifThresh,
                                         by =['logFC'] ,
                                         ascending = [False] ,
                                         absMask = [True] )

    return results_df
