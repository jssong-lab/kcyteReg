#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 09:58:41 2018

@author: afinneg2
"""

from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import os
from collections import OrderedDict
import argparse
import pickle

from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import  squareform

import IOutils

import matplotlib.pyplot as plt 
import seaborn as sns
plt.rcParams["axes.labelsize"] = "x-large"
plt.rcParams["xtick.labelsize"] =  "x-large"
plt.rcParams["ytick.labelsize"] =  "x-large"
plt.rcParams["legend.fontsize"] = "x-large"

###########################################################################################################################
## FUNCTION DEFINITONS
def clusterMap_genesCoReg_TFsCorr( simMat , TFs ,  linkageMethod_TFs = "complete",
                                distMetric_regTargets =  "euclidean"  , linkageMethod_regTargets = "complete", 
                                distMetric_regTargets_kws = {} , 
                                kws_clustGrid = {"z_score" : 1 , "cmap" : "RdBu_r" , 
                                                "center" : 0 ,  "vmax" : 6}  ,
                                figsize = (20,20) ):
    """
    Inputs:
        simMat - dataFrame  rows and columns are genes values store signed gene similarity based on expression
        linkageMethod_SE= "complete"
        distMetric_TFsim = "euclidean"
        linkageMethod_TFs = "complete"
    Returns:
        clustGrid , Y_genes, Y_TFs - seaborn object, linkage matrix for SE , linkage matrix for TFs
    """
    ##############
    # CLUSTER TFs based on gene correlation patterns
    dists = ((1 - simMat.loc[TFs, TFs ].copy().values) + (1 - simMat.loc[TFs, TFs ].copy().values).transpose() ) /2.0  ## insure symmetric
    dists[np.arange( dists.shape[0]), np.arange( dists.shape[0]) ] = 0.0  ## insure 0 self-distance
    TF_dists =  squareform( dists )
    Y_TFs = sch.linkage( TF_dists  , method= linkageMethod_TFs )
    ##############
    ## Cluster Genes based on correlation patterns with TFs
    regTargets = [ x for x in  simMat.index if x not in TFs  ] 
    simMat_cluster  = simMat.loc[ regTargets, TFs  ].copy()  ## rows are regulatory targets cols are TFs
    regTargets_dist = pdist( simMat_cluster.values.copy(), distMetric_regTargets, **distMetric_regTargets_kws  )
    Y_regTargets = sch.linkage( regTargets_dist,  method = linkageMethod_regTargets )
    ################################
    ## MAKE PLOT 
    clustGrid =  sns.clustermap( data = simMat_cluster , figsize =  figsize, row_linkage=Y_regTargets , 
                                col_linkage= Y_TFs , **kws_clustGrid)
    return clustGrid,  Y_regTargets , Y_TFs

###########################################################################################################################
## PARSE COMMAND LINE  
parser = argparse.ArgumentParser(description= "cluster TFs by correlation in expression. Cluster genes by coReg stat" )
parser.add_argument("--geneCorr" , help =  "dataframe of gene correlations" )
parser.add_argument("--genes", default = "" , help =  "file listing gene names(1 per line). If comma sep list take union")
parser.add_argument("--TFs" , help = "File listing TFs"  )
parser.add_argument("--beta" , default = 4, type = int)
#parser.add_argument("--combineStages" , help= "e.g. <superMat_name1>:<Name1>,<Name2>;<sumperMat_name2>:<Name3>")
parser.add_argument("--fo" , help = "name output figure. Extension indiciates figure type. Data files saved with same name but different extension" )
#parser.add_argument("--stagesPlot")
parser.add_argument("--vmin" , type = float)
parser.add_argument("--vmax" , type = float)
parser.add_argument("--linkage_regTargets" , default = "complete")
parser.add_argument("--linkage_TFs" ,  default = "complete")
args = parser.parse_args()

###########################################################################################################################
## LOAD DATA  
corrDF = IOutils.loadDF(args.geneCorr)
genes = list(set().union(*[ IOutils.readListFromFile(x) for x in args.genes.split(",")] ))
genes_missing = [x for x in genes if not (x in corrDF.index)  ]
if len(genes_missing) > 0:
    print("WARNING the following genes are not in correlation data {}".format( genes_missing) )
    genes = [x for x in genes if x in corrDF.index  ]
TFs = list(set().union(*[ IOutils.readListFromFile(x) for x in args.TFs.split(",")] ))
TFs_missing = [x for x in TFs if not (x in corrDF.index)  ]
if len(TFs_missing) > 0:
    print("WARNING the following TFs are not in correlation data {}".format( TFs_missing ) )
    TFs = [x for x in TFs if x in corrDF.index ]

############################################################################################################################# 
## COMPUTE SIGNED ADJACENCY  
genes_and_TFs = sorted(set(genes ).union(set( TFs )) )
corrDF = corrDF.loc[ genes_and_TFs , genes_and_TFs ].copy()
signedAdj_df =  np.sign(corrDF.values) * (corrDF.abs()**(args.beta))

##############################################################################################################################
### Generate Cluster Maps
clustGrid, Y_regTargets, Y_TFs =  clusterMap_genesCoReg_TFsCorr(simMat = signedAdj_df.loc[genes_and_TFs, TFs].copy() ,
                                                            TFs = TFs , 
                                                            linkageMethod_TFs = args.linkage_TFs,
                                                            distMetric_regTargets =  "euclidean" ,
                                                            linkageMethod_regTargets = args.linkage_regTargets, 
                                                            distMetric_regTargets_kws = {} , 
                                                            kws_clustGrid = {"z_score" : None , "cmap" : "RdBu_r" , 
                                                                    "center" : 0 , "vmin" : args.vmin ,"vmax" : args.vmax} )
###############################################################################################################################
## Save figure and data 
clustGrid.savefig( args.fo , format = os.path.splitext(args.fo)[-1][1:])
f = open( os.path.splitext(args.fo)[0] + ".pkl"  , 'wb' )
pickle.dump(obj = {"clustGridData":  clustGrid.data.copy() ,
                   "dendro_regTargets" : clustGrid.dendrogram_row.dendrogram  ,
                   "linkageMat_regTargets":  Y_regTargets ,
                    "dendro_TFs" :  clustGrid.dendrogram_col.dendrogram ,
                   "linkageMat_TFs" : Y_TFs }, 
                file=  f)
f.close()

regTargets =[clustGrid.data.index[i] for i in clustGrid.dendrogram_row.dendrogram["leaves"]  ]
f = open( os.path.splitext(args.fo)[0] + "_targetsDendroOrder.txt"  , 'w'  )  
for x in regTargets:
    f.write(x+ "\n")
f.close()

TFs_colOrder =  [clustGrid.data.columns[i] for i in clustGrid.dendrogram_col.dendrogram["leaves"]]
f = open( os.path.splitext(args.fo)[0] + "_TFsDendroOrder.txt"  , 'w'  ) 
for x in TFs_colOrder:
    f.write(x + "\n")
f.close()

    





