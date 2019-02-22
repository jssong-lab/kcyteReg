#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import sys
import os
from collections import OrderedDict

import networkx as nx
import pygraphviz as pgv
import re

import IOutils

import matplotlib.pyplot as plt 
from matplotlib import patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns  

plt.rcParams["axes.labelsize"] = "xx-large"
plt.rcParams["xtick.labelsize"] = "xx-large"
plt.rcParams["ytick.labelsize"] = "xx-large"
plt.rcParams["legend.fontsize"] = "xx-large"

##############################################################################################################
### FUNCTION DEFINITIONS

#####################################################################
### LOAD DATA
def loadGeneModules(fnames, fnameToModName = lambda x: x):
    """
    fnames - list of strings
    fnameToModName - function: str -> str
    """
    moduleDict =OrderedDict([ (fnameToModName(x), IOutils.readListFromFile(x)) for x in fnames ] )
    return moduleDict

def get_meanVal_modulePairs( scoreDF , rowModules , colModules ):
    """
    Inputs
        rowModules - dictionary {"ModuleName" : <List of Row Elements>}
        colModules - dictionary {"ModuleName": <List of Col Elements>}
    Returns:
        meanScoreDF - rows are keys or rowModules , cols are keys of colModules
    """
    rowElemToModule = pd.Series( {elem : k for k , v in  rowModules.items() for elem in v } )
    colElemToModule = pd.Series( {elem : k for k , v in  colModules.items() for elem in v } )

    scoreDF_byRowMod = scoreDF.groupby( by =  rowElemToModule , axis = 0  )
    meanScoreDF = scoreDF_byRowMod.apply( lambda x: x.groupby(by =  colElemToModule, axis = 1 ).apply( lambda y : y.mean().mean() ) )
    return  meanScoreDF 

###########################################################################
### MAKE REGULATORY FIGURE
def getNodeLayout( nodes_l, nodes_r, edgeDict, x_pad = 0.05 , y_pad = 0.05, arrow_pad = 0.15,
                      ax_pad = 0.01, ax_maxHeight = 0.5 ):
    """
    Inputs:
        nodes_l - list of nodeNames
        nodes_r - list of nodeNames
        edgeDict - OrderedDict([  (( NodeName1 , NodeName2) , edgeParamDict ), ... ])
    Returns: 
        nodes1_params - dictonary of node paramseters: 
                        OrderedDict([ (<NodeName>,{"pos_fig": [0.0 , 0.3 , 0.4, 0.4] ,
                                        "nAttach" : 2 , "attachSide" :"right" ), ... ])
        nodes2_parasms - dictonary of node paramseters: 
                        OrderedDict([ (<NodeName>,{"pos_fig": [0.0 , 0.3 , 0.4, 0.4] ,
                                        "nAttach" : 2 , "attachSide" :"left" ), ... ])
    """
    def parseDotFormat( s, nodes1 , nodes2 ):
        """
        PARSING pydotvis OUTPUT
        Inputs:
            s - string description of networkx AGraph object for 
                bipartite graph
            nodes1 - list nodes composition one part of bipartite graph
            nodes2 - list nodes composition other part of bipartite graph
        Returns:
            nodes1_ordered  ,nodes2_ordered
        """
        nodes = nodes1 + nodes2
        pos_dict = {}
        lines = s.split("\n")
        pair = []
        for line in lines:
            if line.find('--') > 0:
                break
            indicator = [ line.strip().startswith(x) for x in nodes  ]
            if any(indicator):
                node = [ x for x ,y in zip(nodes, indicator )  if y ]
                assert( len(node) ==1) , "Muiltiple Matches"
                node = node[0]
                pair.append(node)
            pos_match = re.match(r'\s*pos="(\d+\.*\d*),(\d+\.*\d*)",$' , line )
            if pos_match is not None:
                pos = (float(pos_match.group(1) ), float( pos_match.group(2)) )
                pair.append(pos)
                assert(len(pair) ==2) , "Names and pos values not paired"
                pos_dict[pair[0]] = pair[1]
                pair = []
        assert( np.all( np.diff(sorted([pos_dict[node][0] for node in nodes1 ])) > 0 )) , "overlapping nodes in nodes1"
        nodes1_ordered = sorted( nodes1 , key = lambda x:  pos_dict[x][0]  )
        assert( np.all( np.diff(sorted([pos_dict[node][0] for node in nodes2 ])) > 0 )) , "overlapping nodes in nodes2"
        nodes2_ordered = sorted( nodes2 , key = lambda x:  pos_dict[x][0]  )    
        return nodes1_ordered  ,nodes2_ordered
    #######################
    ## GET pygraphviz NODE ORDER
    B = nx.Graph()
    B.add_nodes_from(nodes_l , bipartite=0)
    B.add_nodes_from(nodes_r, bipartite=1)
    B.add_edges_from( list(edgeDict.keys()))
    A = nx.nx_agraph.to_agraph(B)
    one = A.add_subgraph(nodes_l ,rank='same')
    two = A.add_subgraph(nodes_r ,rank='same')
    A.layout(prog = "dot")
    nodes_l_ordered, nodes_r_ordered = parseDotFormat( A.to_string(), nodes_l, nodes_r )
    
    #########
    ## SETUP nodes1_params
    ax_height = min( (1.0 -2*y_pad + ax_pad)/float(len(nodes_l_ordered)) - ax_pad ,  ax_maxHeight )
    ax_width = (1.0 - 2*x_pad - arrow_pad)/2.0
    yll_arr = np.arange( y_pad, 1.0, ax_height + ax_pad)
    xll_arr = np.ones(len(nodes_l_ordered), dtype= float)*x_pad
    nodes_l_kws = OrderedDict([])
    for i , node in enumerate( nodes_l_ordered ):
        nAttach = sum( [1 for nodePair in edgeDict.keys() if  nodePair[0] == node ] )
        nodes_l_kws[node] =  {"pos_fig": [xll_arr[i] , yll_arr[i], ax_width, ax_height] , "nAttach" : nAttach, "attachSide" : "right"}
    #########
    ## SETUP nodes2_params
    ax_height = min( (1.0 -2*y_pad + ax_pad)/float(len(nodes_r_ordered )) - ax_pad ,  ax_maxHeight )
    ax_width = (1.0 - 2*x_pad - arrow_pad)/2.0
    yll_arr = np.arange( y_pad, 1.0, ax_height + ax_pad)
    xll_arr = np.ones(len(nodes_r_ordered), dtype= float)*(ax_width + x_pad + arrow_pad )
    nodes_r_kws = OrderedDict([])
    for i , node in enumerate( nodes_r_ordered ):
        nAttach = sum( [1 for nodePair in edgeDict.keys() if nodePair[1] == node ] )
        nodes_r_kws[node] =  {"pos_fig": [xll_arr[i] , yll_arr[i], ax_width, ax_height] , "nAttach" : nAttach, "attachSide" : "left"} 
    return nodes_l_kws, nodes_r_kws
    
    
def makeAxNetwork(nodes_set1_kws , nodes_set2_kws, edges, figsize= (12,12),
                  FancyArrowPatch_kws = dict(arrowstyle = "simple",alpha = 0.7),  facecolor = "whitesmoke",
                  strength_to_mut =  lambda x:x  ):
    """
    nodes_set1_kws - OrderedDict([("node_s1_1", {"pos_fig": [0.0 , 0.3 , 0.4, 0.4] , "nAttach" : 2 , "attachSide" : "right"} ) , ... ])
                       arrows start from these nodes
    nodes_set2_kws - OrderedDict([("node_s2_1", {"pos_fig": [0.6 , 0.0 , 0.4, 0.4] , "nAttach" : 1 , "attachSide" : "left"} ), ...])
                          arrows end on these nodes
    edges = OrderedDict([ (("node_s1_1","node_s2_1"), {"strength" : 40 , "color": "g" , **FancyArrowPatch_kws}), .. ]) 
    
    """
    def makeNode(pos_fig, nAttach , attachSide = "top", fig= None, name = ""):
        """
        Uses figure , axis transformation ideas from:
        https://www.cilyan.org/blog/2016/01/23/matplotlib-draw-between-subplots/
        Inputs:
            pos_fig - position in figure coordinates
        returns - nodeDict - {"ax" : ax , "attachPts" : <list>, "name" :  <name> }
        """
        if fig is None:
            fig = plt.gcf()
        ax = fig.add_axes(pos_fig)
        ########################
        ## setup attach points
        figtr = fig.transFigure.inverted()
        if attachSide.lower() == "top":
            attachPts = [(x , 1.0) for x in np.linspace(0,1 ,nAttach+2)[1:-1] ]
        elif  attachSide.lower() == "bottom":
            attachPts = [(x , 0.0) for x in np.linspace(0,1 ,nAttach+2)[1:-1] ]
        elif attachSide.lower() == "left":
            attachPts = [(0.0, y) for y in np.linspace(0,1 ,nAttach+2)[1:-1] ]
        elif attachSide.lower() == "right":
            attachPts = [(1.0, y) for y in np.linspace(0,1 ,nAttach+2)[1:-1] ]
        else:
            raise Exception("attachSide {} not recognized".format( attachSide ) )
        attachPts_figCoords = list(map( lambda x: figtr.transform(ax.transAxes.transform(x)),  attachPts ) )
        ###########################
        ## aesthetics
        ax.set_xticks([])
        ax.set_yticks([])
        nodeDict= {"ax": ax, "attachPts" : attachPts_figCoords , "name": name }
        return nodeDict        
    def makeDummyAx(pos_fig, visible = True ):
        ax = fig.add_axes(pos_fig)
        ax.set_visible( visible )
        ax.set_xticks([])
        ax.set_yticks([])    
        return ax
    def edgeParams_to_arrowParams(edgeParams, strength_to_mut = lambda x:x ):
        """
        Inputs
            edgeParams - dictionary
        Returns
            arrowParams - dict of kws to mpl.patches.FancyArrowPatch
        """
        if "strength" in edgeParams.keys():
            strength = edgeParams.pop("strength")
            edgeParams["mutation_scale"] = strength_to_mut( strength )
        return edgeParams
    fig = plt.figure(figsize = figsize, facecolor = facecolor, edgecolor = "k", frameon=True)
    makeDummyAx([0,0, 0.001 , 0.001] ) ## hack
    makeDummyAx([1.0,1.0, 0.001 , 0.001] )
    ############
    ## DRAW NODES
    nodes_set1 = OrderedDict([])
    for name , nodeKws in nodes_set1_kws.items():
        nodes_set1[name] = makeNode(fig= fig, **nodeKws)

    nodes_set2 = OrderedDict([])
    for name , nodeKws in  nodes_set2_kws.items():
        nodes_set2[name] =  makeNode(fig= fig, **nodeKws)    
    ############
    ## DrawEdges
    ## Order edges
    nodes_set1_ordered = list( nodes_set1.keys())
    nodes_set2_ordered = list( nodes_set2.keys())
    edges_ordered_keys = sorted(list(edges.keys()) , 
                           key = lambda x: (nodes_set1_ordered.index(x[0]),nodes_set2_ordered.index(x[1])) )
    edges_ordered = OrderedDict([ (k,  edges[k]) for k in  edges_ordered_keys ])
    for nodePair, edge_kws in edges_ordered.items():
        node1 , node2 = nodePair
        arrowParams = edgeParams_to_arrowParams(edge_kws , strength_to_mut =  strength_to_mut )
        arrowParams.update(FancyArrowPatch_kws)
        arrow = patches.FancyArrowPatch( nodes_set1[node1]["attachPts"].pop(0) ,
                                        nodes_set2[node2]["attachPts"].pop(0) ,
                                       transform = fig.transFigure , 
                                       **arrowParams)
        fig.patches.append(arrow)
    return fig , nodes_set1 , nodes_set2
    
###############################################################
## Decorate nodes
def plot_moduleExprByStage( plotData, genes , stages = [1,2,3,4,5,6,7] , ax = None,
                              kws_lineplot= {"ci":"sd" }, figsize = (7.5,4), xlabel = "",
                          ylabel = "", title = "" , title_kws = {'fontsize': "large"}):
    """
    Inputs
        plotData - DataFrame , columns are genes, index is stages , values are log10( imputed_tpm + 1 )
        genes - list of genes in module (column subset)
        stages - index subset
    Returns 
        ax
    """
    if ax is None:
        _ ,ax = plt.subplots(nrows = 1 , ncols = 1, figsize = figsize)
    plotData = plotData.loc[stages , genes].copy()
    plotData_zscore = plotData.subtract( plotData.mean(axis =0) , axis = 1).divide( plotData.std(axis = 0), axis = 1  )
    plotData_zscore.index.set_names("stage ID" , inplace = True)
    plotData_zscore = plotData_zscore.reset_index()
    plotData_zscore= plotData_zscore.melt(id_vars = ["stage ID"] ,var_name = "gene" , value_name= "Expression z-score" )
    ax = sns.lineplot(x="stage ID", y="Expression z-score", data=plotData_zscore,ax = ax, marker = "o" ,  **kws_lineplot)
    ## Make sure all stages are labeled
    ax.set_xticks( sorted(  plotData_zscore.loc[: , "stage ID"].unique()) )
    ax.set_xticklabels( sorted(  plotData_zscore.loc[: , "stage ID"].unique()) )
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title( title, **title_kws  )
    return ax

def add_moduleText(ax, geneList, name ,wordWrap=4 , moduleName_lookup = None , 
                  text_kws = { "x" :0.5,  "y": 0.5 , "ha": "center" , "va" : "center", "fontsize" : 14 }):
    if  moduleName_lookup is None:
         moduleName_lookup = {name: name}
    
    if len(geneList) > wordWrap:
        nRows = -((-len(geneList))//wordWrap)
        TFstr = "\n".join( [ ", ".join( geneList[i*wordWrap : (i+1)*wordWrap] ) for i in range(nRows)] )
    else:
        TFstr = ", ".join(geneList)
    text = "{:s}\n{}".format(moduleName_lookup[name], TFstr  )  
    ax.text( s = text, transform=ax.transAxes, **text_kws)
    return ax

###############################################################################
## Utilities
def cdf_empirical(obs ):
    """
    Inputs:
        obs - 1d array
    Returns:
        x, cumProb
    """
    obs_sorted = np.array( sorted(list(obs)) )
    cumProb = np.arange(1 , len(obs_sorted )+1, 1) / float(len(obs_sorted ) )
    return obs_sorted,  cumProb 

def make_pairAvgRegFig( meanScoreDF,
                       thresh_regMagnitude = 0.12, figsize = (7.5, 6),
                      xticklabel_size = 'medium' , yticklabel_size = 'medium'):

    fig , axes = plt.subplots( nrows = 2 , ncols = 2 , figsize = figsize )
    xticklabelsize_orig =  plt.rcParams["xtick.labelsize"]
    yticklabelsize_orig = plt.rcParams["ytick.labelsize"]
    plt.rcParams["xtick.labelsize"] =  xticklabel_size 
    plt.rcParams["ytick.labelsize"] = yticklabel_size 
    subplot2grid_shape = (2,5)
    
    ## module pair avg reg
    ax = plt.subplot2grid( shape =subplot2grid_shape , loc = (0,0), colspan=3 )
    ax = sns.heatmap( meanScoreDF, cmap = plt.cm.bwr , center = 0, ax =ax, xticklabels = 1, yticklabels=1)
    ax.set_xlabel("TF modules")
    ax.set_ylabel("Gene modules")
    
    ## Inverse CDF
    ax =plt.subplot2grid( shape =subplot2grid_shape , loc = (0,3),colspan=2 )
    x,y = cdf_empirical(np.ravel(np.abs(meanScoreDF.values)) )
    ax.plot(y*100,x )
    ax.axhline( thresh_regMagnitude, c = 'r' )
    ax.set_xlabel("Percentile")
    ax.set_ylabel("abs( mean regulation score )")
    
    ## module pair avg reg masked
    ax =plt.subplot2grid( shape =subplot2grid_shape , loc = (1,1),colspan=3 )
    ax.set_facecolor("grey")
    ax = sns.heatmap( meanScoreDF, cmap = plt.cm.bwr, center=0, ax=ax, mask= meanScoreDF.abs() < thresh_regMagnitude,
                        xticklabels = 1, yticklabels=1)
    ax.set_xlabel("TF modules")
    ax.set_ylabel("Gene modules")    

    plt.rcParams["xtick.labelsize"] = xticklabelsize_orig
    plt.rcParams["ytick.labelsize"] = yticklabelsize_orig
    fig.tight_layout(h_pad = 2 , w_pad = 1 )
    return fig
    
