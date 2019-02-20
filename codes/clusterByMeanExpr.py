#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 12:22:02 2018

@author: afinneg2
"""
from __future__ import  print_function
import numpy as np
import pandas as pd
from cycler import cycler

## H cluster funcs 
from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import  squareform

import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams["axes.labelsize"] = "x-large"
plt.rcParams["xtick.labelsize"] =  "large"
plt.rcParams["ytick.labelsize"] =  "large"
plt.rcParams["legend.fontsize"] = "x-large"

########################################################################################################################
## UTILITY FUNCTIONS 
def plotExprVsStage(df, ax = None, figsize = (7.5 , 5), colors = None , styles = None, bbox_to_anchor = (1.0 ,1.0) ,
                   legend_ncol= 2): 
    """ 
    FUNCTION NOT YET REVISED
    df - rows are genes, columns are stages
    """
    if colors is None:
        colors = plt.cm.tab10
    if styles is None:
        styles = ["-"]
    
    if ax is None:
        fig = plt.figure(figsize = figsize)
        ax = fig.add_subplot(1,1,1)
    
    ax.set_prop_cycle( cycler("color" , colors.colors*len(styles) ) + cycler("linestyle", list(np.repeat(styles, len(colors.colors) )) ) ) 
    xvals = np.arange( len(df.columns.values) )
    for gene , exprSeries in df.iterrows():
        ax.plot(xvals , exprSeries.values, label = gene  )
    ax.set_xticks(xvals)
    ax.set_xticklabels(df.columns.values, rotation = 45) 
    
    ax.legend(loc = "upper left", bbox_to_anchor = bbox_to_anchor , ncol = legend_ncol)
    
    return ax, colors
    
def standardize_DF(df, rows = True):
    """
     FUNCTION NOT YET REVISED
    rows- if true : standardize rows else standardize columns
    """
    if rows:
        means = df.mean(axis = 1)
        stds = df.std(axis = 1)
        result =   df.subtract(means , axis = 0 ).divide(stds , axis =0  )
    else:
        means = df.mean(axis = 0)
        stds = df.std(axis = 0)
        result = result =   df.subtract(means , axis = 1 ).divide(stds , axis =1 )
    return result

def meanCenter_DF(df , rows = True):
    """
     FUNCTION NOT YET REVISED
    """
    if rows:
        means = df.mean(axis = 1)
        result = df.subtract(means , axis = 0 )
    else:
        means = df.mean(axis = 0)
        result = result =   df.subtract(means , axis = 1 )
    return result

def logFC_oneVsRest(df , rows = True):
    """
     FUNCTION NOT YET REVISED
    """
    def logFC_oneVsRes_series(series):
        series_sum  = series.sum()
        seriesLen = series.shape[0]
        series_result = series.apply(  lambda x : np.log2(x *(seriesLen  - 1) / ( series_sum - x ) ) )
        return  series_result
    if rows:
        result = df.apply( logFC_oneVsRes_series , axis =1)
    else:
        result = df.apply( logFC_oneVsRes_series , axis =0)
    return result

def logFC_oneVsMax(df , rows = True):
    """
    FUNCTION NOT YET REVISED
    """
    def logFC_oneVsMax_series(series):
        series_max  = series.max()
        series_result = series.apply(  lambda x : np.log2(x /  series_max  ) )
        return  series_result
    if rows:
        result = df.apply( logFC_oneVsMax_series , axis =1)
    else:
        result = df.apply( logFC_oneVsMax_series , axis =0)
    return result


#######################################################################################################################
##### PLOTTING

def clusterHmap(data, linkageMat, data_secondary= None, 
                cmap = plt.cm.bwr, cmap_secondary = None , 
                cbarLabel = "" , cbarLabel_secondary = "" ,   xtickLabelFontsize = 12 ,
                ytickLabelFontsize = 12 , grid_lw=2 ,
               cbarLabelSize = 12, cbarTickLabelSize = 12 , xlabelHeight = 0.15 ,
               ylabelWidth = 0.15 , matrixWidth = 0.7 ,color_threshold = 0.7 ,matrixHeight_secondary = 0.1,
               figsize=( 15,15), colorMinMaxPct = None, dendroColorPalette  = None ):
    """
    Create a figure of the form:
        |-------------------|
        |   dendrogram      |
        ----------------------
        |     Heat map      |
    Where heat map columns correspond to plotData rows and are ordered according to the 
    clustering in the dendrogram
    Inputs:
        data - dataFrame. Rows are Random variables (e.g. genes) columns are observations (e.g. differentiation stage)
        linkageMat - output of sch.linkage where clustering is done on the rows or data
        data_secondary - None or 1 column dataFrame with index matching that of linkageMat -
        
    Returns:
        fig
    """
    class MidpointNormalize(colors.Normalize):
        """
        FROM : http://chris35wills.github.io/matplotlib_diverging_colorbar/
        Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
    
        e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
        """
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            colors.Normalize.__init__(self, vmin, vmax, clip)
    
        def __call__(self, value, clip=None):
            # I'm ignoring masked values and all kinds of edge cases to make a
            # simple example...
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
    ##########################################################################################################
    ##  SET UP AXIS POSITIONS (DENDROGRAM , HEATMAP , AND HEATMAP SECONDARY) 
    fig = plt.figure(figsize=figsize )
    ## set up axes posiitions. Coordinates are left , bottow,  width , height
    position_AxDendro = [ylabelWidth, xlabelHeight + matrixWidth + 0.01,  matrixWidth, 1 - ( xlabelHeight + matrixWidth + 0.01 ) ]
    if data_secondary is None:
        position_AxHeat_main  = [ylabelWidth, xlabelHeight, matrixWidth, matrixWidth ]
        position_axCbar_main = [ylabelWidth+ matrixWidth + 0.02, xlabelHeight, 0.02, matrixWidth]
    else:
        position_AxHeat_main  =[ylabelWidth, xlabelHeight + matrixHeight_secondary , matrixWidth, matrixWidth - matrixHeight_secondary ]
        position_axCbar_main = [ylabelWidth+ matrixWidth + 0.02, xlabelHeight  +  matrixHeight_secondary, 0.02, matrixWidth- matrixHeight_secondary ]
        position_AxHeat_secondary = [ ylabelWidth, xlabelHeight  , matrixWidth, matrixHeight_secondary ]
        position_axCbar_secondary = [ylabelWidth+ matrixWidth + 0.02, xlabelHeight  - 0.001, 0.02, matrixHeight_secondary - 0.001]
    #########################################################################################################
    ## PLOT Dendrogram 
    axDendro = fig.add_axes( position_AxDendro ) 
    if dendroColorPalette is not None:
        hierarchy.set_link_color_palette( dendroColorPalette )
    Z = sch.dendrogram( Z = linkageMat, color_threshold = color_threshold, ax = axDendro )

    axDendro.set_xticks([])
    axDendro.set_yticks([])
    #########################################################################################################
    ### PLOT PRIMARY HEATMAP
    axmatrix = fig.add_axes(position_AxHeat_main )
    heatmap_vals = data.copy().values.transpose()
    RVnames = np.array(list(data.index))
    observNames = np.array(list(data.columns))
    idx = Z['leaves']  ## an array whre Z['leaves'][i] stores the leaf at ith position when traversing
                         ## dendrogram leaves from left ot right
    heatmap_vals = heatmap_vals[:, idx ]   ## Reorder heatmap columns by clustering
    if colorMinMaxPct is not None:
        vmin = np.percentile(np.ravel( heatmap_vals), q  = colorMinMaxPct[0] )
        vmax = np.percentile(np.ravel( heatmap_vals), q  = colorMinMaxPct[1] )
        im = axmatrix.matshow(heatmap_vals, aspect='auto', origin='upper', cmap= cmap , vmin = vmin , vmax = vmax,
                                  norm=MidpointNormalize(midpoint=0,vmin=vmin, vmax=vmax ))
    else:
        im = axmatrix.matshow(heatmap_vals, aspect='auto', origin='upper', cmap= cmap,
                                  norm=MidpointNormalize(midpoint=0,vmin=None, vmax=None ) )
    ## yticklabels
    axmatrix.set_yticks( np.arange(len(observNames))  )
    axmatrix.set_yticklabels(observNames, fontsize =ytickLabelFontsize)
    ## draw grid
    if  grid_lw >0 :
        axmatrix.grid(color = 'w' , lw = grid_lw, which = 'minor')
    ## add color bar  
    axcolor = fig.add_axes( position_axCbar_main )
    cbar  = plt.colorbar(im, cax= axcolor )
    cbar.ax.tick_params(labelsize = cbarTickLabelSize )
    cbar.set_label(label = cbarLabel , fontsize = cbarLabelSize )
    #########################################################################################################
    ### PLOT SECONDARY HEATMAP
    if data_secondary is not None:
        cmap_secondary.set_bad("grey")
        axmatrix_secondary = fig.add_axes(  position_AxHeat_secondary ) 
        heatmap_vals_secondary =  data_secondary.loc[ data.index.values[  idx  ], : ].values.transpose()
        im_secondary =  axmatrix_secondary.matshow( heatmap_vals_secondary, aspect='auto', origin='upper', cmap= cmap_secondary ,
                                                     norm=MidpointNormalize(midpoint=0,vmin=None, vmax=None )   )
    
        axmatrix_secondary.set_yticks( np.arange(len(data_secondary.columns))  )
        axmatrix_secondary.set_yticklabels( list(data_secondary.columns) ,  fontsize =ytickLabelFontsize )
        
        ## add cbar_secondary
        axcolor_secondary = fig.add_axes( position_axCbar_secondary ) 
        cbar_secondary  = plt.colorbar( im_secondary, cax=axcolor_secondary   )
        cbar_secondary.ax.tick_params(labelsize = cbarTickLabelSize )
        cbar_secondary.set_label(label = cbarLabel_secondary , fontsize = cbarLabelSize )
    ## xticklabels
    if data_secondary is None:
        axmatrix.xaxis.set_ticks_position("bottom")
        axmatrix.set_xticks(np.arange(len(RVnames)) )
        axmatrix.set_xticklabels(RVnames[idx], rotation = 90 , fontsize =xtickLabelFontsize )
    else:
        axmatrix.xaxis.set_ticks_position("bottom")
        axmatrix.xaxis.set_ticks([])
        axmatrix_secondary.xaxis.set_ticks_position("bottom")
        axmatrix_secondary.set_xticks(np.arange(len(RVnames)) )
        axmatrix_secondary.set_xticklabels(RVnames[idx], rotation = 90, fontsize =xtickLabelFontsize )
        
    return fig ,Z
    
def correlCluster_lineplot(data ,linkageMethod = "complete", fclust_kwargs = {"t": 0.7 , "criterion" : "distance"},
                           exprTransform = "z-score" ,
                           ncols = 2, figWidth = 20 , axHeight = 4,
                           cmap = plt.cm.tab10 ,styles = ["-", "--"] , wspace = 0.7 , hspace = 0.5 ,
                           ylabel = ""):
    """
    Function copied from analysis2/scripts/clusterByMeanExpr_funcs.py  - has not yet been revised as part of code cleanUp
    """
    ## cluster data ###
    cosineDists = 1 - data.transpose().corr() ## 1 - pearson correlation of rows of data
    obsNames = data.index.values
    featNames = data.columns.values

    Y = sch.linkage(squareform(cosineDists.values), method= linkageMethod)## output is hierarchical
                                        ## clustering reported as a linkage matrix. Linkage matrix is
                                        ## (nbObservations -1  ,4 ) array where rows index iterations of cluster
                                        ## joining and cols 0 ,1 store joined cluster idx,  col2 indictes distance
                                        ## between joined clusters and col 3 is number of origonal observations in
                                        ## the new cluster
                                        ## Noten that because input is a matrix of distance  
                                        ##. sch.linkage API requires that distances are 
                                        ## represented as 1d array using scipy.spatial.distance.squareform on 
                                        ## a symmetric 2d distance matri
    fclust = sch.fcluster(Z = Y , **fclust_kwargs )
    fclustIdxs =  np.unique( fclust )

    ## allow different transformations of the 1d array of expression values before plotting
    if exprTransform == "z-score":
        data_plot = standardize_DF(data, rows = True)
    elif exprTransform == "logFC_oneVRest"   :
        data_plot = logFC_oneVsRest(data , rows = True)
    elif exprTransform == "logFC_oneVsMax" :
        data_plot = logFC_oneVsMax(data , rows = True)
    elif exprTransform =="meanCenter":
          data_plot  = data.apply(lambda x:  x - x.mean() , axis = 1)
    elif exprTransform in ["None" , "none" , ""] :
        data_plot = data
    else:
       raise ValueError( "exprTransform value {} not recognized".format(exprTransform) )

    nrows = int(np.ceil( len(fclustIdxs) / ncols ))
    fig , axes = plt.subplots(figsize = (figWidth , axHeight* nrows) , nrows = nrows , ncols= ncols)
    axes = np.ravel(axes)
    genesByClust  =[]
    for fclustIdx, ax in zip(fclustIdxs, axes ):
        genes = np.copy(obsNames[fclust== fclustIdx])
        genesByClust.append(genes)
        plotExprVsStage(data_plot.loc[ genes , : ].copy() , ax = ax,
                        colors =  cmap , styles = styles, bbox_to_anchor = (1.0 ,1.0))
        if ax.is_first_col():
            ax.set_ylabel(ylabel)
        if not (ax.is_last_row() or (fclustIdx ==fclustIdxs[-1] ) ):
            ax.tick_params(labelbottom = "off")
    fig.subplots_adjust(hspace = hspace , wspace = wspace)


    return fig , genesByClust

