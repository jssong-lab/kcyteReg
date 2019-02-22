#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 17:11:39 2018

@author: afinneg2
"""

from __future__ import division
from __future__ import print_function

import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import sys
import os
import argparse
import pickle


import matplotlib.pyplot as plt
import seaborn as sns  
plt.rcParams["axes.labelsize"] = "x-large"
plt.rcParams["xtick.labelsize"] =  "medium"
plt.rcParams["ytick.labelsize"] =  "medium"
plt.rcParams["legend.fontsize"] = "x-large"

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42

###############################################################################################
## GLOBALS
cmap = plt.cm.tab10
colorNames = ["blue" , "orange" , "green" , "red" , "purple" , "brown" , "pink" , "gray" , "yellow" , "cyan"]

###############################################################################################
## FUNCTION DEFINITION
def fcluster_clustMap(  clustMap_data , linkMat_rows, t_rows, depth_rows,
                        linkMat_cols , t_cols, depth_cols, cmap= cmap, colorNames = colorNames,
                     figsize = (10.5,10.5)):
    
    ############################################################
    ## CALL FLAT ROW CLUSTERS
    fclust_rows = sch.fcluster( linkMat_rows, t =t_rows, depth = depth_rows )
    rowClusters_df = pd.DataFrame(fclust_rows[ :, None], index =  clustMap_data.index , columns = ["inconsistClust_IDs" ] )
    rowClusters_df["Gene module"] = rowClusters_df.loc[: ,"inconsistClust_IDs" ].map( lambda x: cmap( (x-1) % len(cmap.colors) )  )
    rowClusters_df["Gene module names"] = rowClusters_df.loc[:, "inconsistClust_IDs" ].map(
                                                        lambda x: colorNames[(x-1) % len(cmap.colors) ] + " {:d}".format((x-1)//len(cmap.colors) + 1) )

    ############################################################
    ## CALL FLAT COL CLUSTERS
    fclust_cols = sch.fcluster(linkMat_cols , t = t_cols, depth = depth_cols )
    colClusters_df = pd.DataFrame(fclust_cols[ :, None], index =  clustMap_data.columns , columns = ["inconsistClust_IDs" ] )
    colClusters_df["TF module"] = colClusters_df.loc[: ,"inconsistClust_IDs" ].map( lambda x: cmap( (x-1) % len(cmap.colors) )  )
    colClusters_df["TF module names"] = colClusters_df.loc[:, "inconsistClust_IDs" ].map(
                                                        lambda x: colorNames[(x-1) % len(cmap.colors) ] + " {:d}".format((x-1)//len(cmap.colors) + 1) )

    #########################################################################
    ## CONSTRUCT CLUSTER MAP
    g = sns.clustermap(  clustMap_data, row_linkage = linkMat_rows, col_linkage= linkMat_cols, cmap =  plt.cm.bwr, center = 0 ,
                      figsize =  figsize , row_colors = rowClusters_df.loc[: ,  ["Gene module"]].copy() ,
                      col_colors = colClusters_df.loc[: ,  ["TF module"]].copy())
    return g ,  rowClusters_df ,  colClusters_df

if __name__ == "__main__":
    ###################################################################################################
    ## PARSE CMD LINE
    parser = argparse.ArgumentParser( description= "Color rows and columns of gene/TF cluster Grid by modules (flat clusters). \
Flat clusters are created with inconsistency statistic")
    parser.add_argument("--clustMap" , help =  'pickled dictionary with keys "clustGridData" ,  "linkageMat_regTargets" "linkageMat_TFs" ')
    parser.add_argument( "--depth_rows", default = 4, type  = int )
    parser.add_argument("--thresh_rows" , default = 2.15 , type  = float  )
    parser.add_argument("--depth_cols" , default = 4 , type  = int  )
    parser.add_argument( "--thresh_cols" ,default = 0.85 , type = float )
    parser.add_argument("--outDir" , required =True )
    args = parser.parse_args()
    
#    exampleArgs = "--clustMap ./NHEKP_lfcThresh0.25_cluster.t10_geq1pctAllTissue/regLink.average_TFlink.average_beta4.pkl \
#    --depth_rows 4 --thresh_rows 2.15 --depth_cols 2 --thresh_cols 0.75 --outDir tmp"
#    args = parser.parse_args(exampleArgs.split())
    
    #######################################
    ## LOAD DATA
    f = open(args.clustMap , 'rb')
    data = pickle.load(f)
    f.close()    
    cgData = data['clustGridData']
    linkMat_rt = data['linkageMat_regTargets']
    linkMat_TF = data['linkageMat_TFs']
    ##########################################################
    ## CALL FLAT CLUSTERS
    g, rowClusters_df, colClusters_df = fcluster_clustMap(  clustMap_data = cgData.copy() ,
                                                    linkMat_rows = linkMat_rt, t_rows = args.thresh_rows, depth_rows= args.depth_rows ,
                                                    linkMat_cols=  linkMat_TF , t_cols = args.thresh_cols, depth_cols = args.depth_cols , 
                                                    cmap= cmap, colorNames = colorNames,
                                                    figsize = (10.5,10.5) )
    #########################################################################
    ## WRITE RESULTS
    os.makedirs(args.outDir , exist_ok=True )
    g.savefig(os.path.join(args.outDir , os.path.splitext(os.path.basename(args.clustMap))[0]+ "_modules.pdf" ) , format ="pdf")
    ## Write gene modules
    geneModule_names = sorted(rowClusters_df.loc[: , "Gene module names"].unique() )
    for name in geneModule_names:
        genes = list(rowClusters_df.loc[rowClusters_df.loc[: , "Gene module names"] == name, : ].index )
        f =open(os.path.join(args.outDir , "genes_" + name.replace( " ", "_" ) +  ".txt" ) , 'w')
        for gene in genes:
            f.write(gene + "\n")
        f.close()
        
    ## Write TF modules
    TFmodule_names = sorted(colClusters_df.loc[: , "TF module names"].unique() )
    for name in TFmodule_names:
        TFs = list(colClusters_df.loc[colClusters_df.loc[: , "TF module names"] == name, : ].index )
        f =open(os.path.join(args.outDir ,  "TFs_" + name.replace( " ", "_" ) + ".txt" ), 'w'  )
        for TF in TFs:
            f.write(TF + "\n")
        f.close()


    
    
