#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 16:17:22 2018

@author: afinneg2
"""

from __future__ import  division
import numpy as np
import pandas as pd
import sys
from collections import OrderedDict
import argparse
from  magic.MAGIC import impute_fast ,  compute_markov
import time

#### Function definitions
def magic_progressive_t(data_df, dist_data_df , t_list, 
                        k = 30 , ka = 10 , epsilon =1 , rescale= 99,
                       fo_imputed ="" ,assessImputFunc =None,
                       runTimes  =False):
    """
    data - rows are observations (cells) columns are variables (genes)
    dist_data. - data used to compute distances between observations (cells) 
                that define the markov transition matrix. Rows are cells
                columns are coordinate in a reduced dimensional representation
                (ie if dimension reduction by PCA with 20 pc components 
                dist_data.shape = (nbCells , 20 ))
    assessImputFunc  - function taking:  imputedData ,
                        and returning imputationQualityMetric
    """
    data = data_df.values
    dist_data = dist_data_df.values
    if runTimes:
        runTimes_dict = OrderedDict([])
        startcalc = time.time()
    W = compute_markov(dist_data , knn =k , epsilon = epsilon , distance_metric = 'euclidean' , knn_autotune =ka)
    if runTimes:
        endcalc =  time.time()
        runTimes_dict["time:compute_markov (s)"] = "{:.2f}".format( endcalc - startcalc )
    W_t = None
    tprev = None
    qualityMetrics = []
    for t in t_list:
        print("working on t={}".format(t))
        if runTimes:
            startcalc=time.time()
        ## do imputation
        data_new , W_t = impute_fast(data, L = W, t = t, L_t = W_t, 
                                     tprev = tprev, rescale_percent = rescale )
        if runTimes:
            endcalc = time.time()
            runTimes_dict["time:imputeFast_t={:d} (s)".format(t)] = "{:.2f}".format(endcalc - startcalc)
        tprev = t
        ## write result
        if fo_imputed:
            pd.DataFrame( data_new , index = data_df.index.values , 
                         columns = data_df.columns.values  ).to_csv(fo_imputed.format(t),
                                                                   float_format = "%.4f")
        ## assess imputation quality
        if assessImputFunc:
            qualityMetric =  assessImputFunc(data_new)
            qualityMetrics.append((t , qualityMetric ))
    if runTimes:
        return qualityMetrics , runTimes_dict
    else:
        return qualityMetrics
    
def libNorm_MAGIC(rawCounts, returnLibSizes = False):
    """
    Perfom MAGIC normalization of each cell's library size. Normalization is
    Dnew_i,j = [D_{i,j} / sum_j(D_{i,j})] * median_i(sum_j(D_{i,j})
    Inputs:
        rawCounts - dataFrame rows are cells , cols are genes
    Reurns:
        rawCounts  - normalized dataFrame
    """
    libSize = rawCounts.sum(axis=1)
    medianLibSize = np.median(libSize.values)
    rawCounts =  rawCounts.divide(libSize , axis = 0) * medianLibSize 
    if returnLibSizes:
        return rawCounts , libSize , medianLibSize 
    else:
        return rawCounts

if __name__=="__main__":

    ###########################################################################
    parser = argparse.ArgumentParser(description='Do magic imputation for range of t vals')
    ### file names
    parser.add_argument('--fi_expr', help='csv file of libSizeNormalized count data. \
    Rows index cells, cols index genes', required=True)
    parser.add_argument('--fi_distData' , help = "csv file of cell coordinates used to compute cell cell distance.  \
     Rows index cells,  cols index coordinates. If provided npc, PCAlogExpr, fo_details are ignored.  \
     If not provided coordinates used to compute distances will be ontained from fi via PCA decomponsition."
     , default = '' )
    parser.add_argument('--fo_imputed' , help="eg. baseName_t{}.csv. {} filled in with appropriate t. \
    If unspecified don't save imputed counts, just compute quality metric", required=True)
    parser.add_argument('--fo_details' , help= 'file for output of details on run parameters and timing)', required  =True)
    ### filtering  / normalization of raw data
    parser.add_argument("--geneRows" , action = 'store_true' , help = "indicate if fi_expr has genes in rows. \
    In this case transpose is taken." )
    parser.add_argument('--minFracExpr' , type = float,  default = 0.0 , 
                        help = "genes expressed in fewer that this fraction of cells are not imputed")
    parser.add_argument('--libNormalize' , action = 'store_true' , help ="Perform magic library size normalization on fi_expr")
    ## Magic Parameters
    parser.add_argument('--tList', help='comma separated list of t values', required = True)
    parser.add_argument('-k' , type = int , default = 30)
    parser.add_argument('--ka' , type = int , default = 10)
    parser.add_argument( '--epsilon', type = int , default = 1)
    parser.add_argument('--rescale' , type =int , default = 99)
    
    args = parser.parse_args()
    tList = [int(t) for t in args.tList.split(",")]
    
    ####### Load expression data  #############################################
    print('Loading cell coordinates for computing distances')
    dist_data_df = pd.read_csv(args.fi_distData, index_col = 0, header = 0)
    print('Loading expression Data')
    data_df = pd.read_csv(args.fi_expr , index_col = 0, header = 0)
    if args.geneRows:
        print("tranposing fi_expr so that genes are in columns")
        data_df = data_df.transpose()
    ### check the cellNames,labeling columns of data_df and dist_data_df agree,
    assert np.all( dist_data_df.index.values == data_df.index.values)
    print("shape of loaded expression data is {}".format(data_df.shape))
    ### filter genes and lib normalize
    if args.minFracExpr > 0 :
        print("Calculating fraction of cells expressing each gene")
        fracExpressed = (data_df > 0).sum(axis = 0) /  data_df.shape[0]
        geneMask = fracExpressed  >= args.minFracExpr 
    if args.libNormalize:
        print("performing MAGIC lib size normalization")
        startcalc = time.time()
        data_df , _ , medianLibSize= libNorm_MAGIC(data_df, returnLibSizes = True)
        endcalc = time.time()
        libNormTime = endcalc - startcalc
    else:
        libNormTime = "NA"
        medianLibSize = "NA"
    if args.minFracExpr > 0:
        print("Removing genes with fraction of cell expression less than {:.4f}".format(args.minFracExpr))
        #data_df= data_df.loc[: , geneMask].copy()  ## change to below method in attempt to reduce mem
        data_df.drop(labels = geneMask.index[~geneMask].values, axis=1   ,inplace = True )
        print("{:d} genes remain".format(data_df.shape[1]))
    
    runDetails = OrderedDict([("fi_expr", args.fi_expr) , 
                             ("fi_distData" , args.fi_distData), 
                             ("k" , args.k) ,
                             ("ka" , args.ka) ,
                             ("epsilon" , args.epsilon) , 
                             ("rescale" , args.rescale),
                             ("nCell" , data_df.shape[0] ),
                             ("nGene" , data_df.shape[1]),
                             ("minFracExpr" , args.minFracExpr) , 
                             ("medianLibSize" , medianLibSize),
                             ("time:libNorm" , libNormTime)] )
    
    ## Run magic
    _, runTimes_dict= magic_progressive_t(data_df, dist_data_df , t_list = tList,
                                            k = args.k , ka= args.ka , epsilon = args.epsilon,
                                            rescale = args.rescale , fo_imputed = args.fo_imputed,
                                            runTimes  =True )
    
    ### write runtimes to file
    runDetails.update( runTimes_dict ) 
    runDetails_str = "\t".join([ "{}:{}".format(key, value) for key, value in runDetails.items() ])
    f =open(args.fo_details , 'a')
    f.write(runDetails_str + "\n")
    f.close()

