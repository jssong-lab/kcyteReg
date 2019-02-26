#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 13:12:14 2018

@author: afinneg2
"""

from __future__ import  division
import numpy as np 
import pandas as pd

import argparse
import os
import re
import itertools
from textwrap import wrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["axes.labelsize"] = "large"
plt.rcParams["axes.titlesize"] = "large"
plt.rcParams["xtick.labelsize"] =  "large"
plt.rcParams["ytick.labelsize"] =  "large"
plt.rcParams["legend.fontsize"] = "large"
mpl.rcParams["pdf.fonttype"] = 42

##### FUNCTION DEFINITIONS #####################################################
######################
### PARSING GO FILES 
  
def add_p(fnameTerms_orig, fnameCol = "fname" , TermCol = "GOterm", pVal_col = "Benjamini" ):
    
    def get_pVal( term , fname, pVal_col = "Benjamini"  ):
        troubleChrs = {'(' :  '\(', ')' : '\)'}
        for chr_find , chr_replace in troubleChrs.items(): term = term.replace(  chr_find , chr_replace  )
    
        GO_df= pd.read_csv( fname , sep = "\t", header = 0 )
        mask =  GO_df.loc[ : , "Term" ].map( lambda x: True if (re.match( term+ ".*" , x) is not None) else False )
        if mask.sum() > 1:
            print("More than 1 match found for term {} in file {} using first hit".format( term , os.path.basename(fname)  ))
        if mask.sum() > 0:
            p = GO_df.loc[mask , pVal_col ].values[0]
        else:
            p = np.nan
        return p
    
    p_series = pd.Series(index= fnameTerms_orig.index, dtype = float )
    for idx, row in fnameTerms_orig.iterrows():
         p_series.loc[idx] =  get_pVal(row.loc[TermCol], row.loc[fnameCol],  pVal_col = pVal_col )
    fnameTerms_orig["p"] = p_series.loc[ fnameTerms_orig.index  ]
    return  fnameTerms_orig
    
##### FUNCTION DEFINITIONS #####################################################
######################
### PLOTTING 
def make_barPlot( plotDF,y = "GOterm_short" , x = "-log10p_orig" ,addHatches = True, figsize = (7.5,3), 
                 n_ax =1, title = None, asterisk_col = None, asterisk_thresholds = [] ,
                 horizontal = False, xlabel = "", ylabel = "",  asterisk_pad = 0.1,
                 xticklabel_rot = 0 , yticklabel_rot = 0):
    """
    Inputs
        plotDF - data frame. Each row is a term
        y - (str) column name of plotDF
        x - (str) column name of plotDF
        asterisk_col - (str) column name of plotDF, Significance is indicated by comparing
                        values in this column to a list of thresholds
        asterisk_tresholds - a list of threholds on values in asterisk_col
        horizontal  - (bool) if True make horizontal bar plot
    """
    def get_barPosition( patch  , loc = "topCenter", pad_top = 0.1, pad_right = 0.0 ):
        if loc == "topCenter":
            center = patch.get_x() + patch.get_width()/2.0 + pad_right
            top =patch.get_y()+  patch.get_height() + pad_top
            x_out , y_out =  center , top
        elif loc == "rightCenter":
            center = patch.get_y()+  patch.get_height()/2.0 + pad_top
            right = patch.get_x() + patch.get_width() + pad_right
            x_out , y_out = right , center
        else:
            raise NotImplementedError
        return (x_out , y_out)
    def addHatchColumn( color_cluster_df , hatchOrder =['/', '//' , '+' , '-' , 'x'] ):
        """
        color_cluster_df - has columns ClusterName
        """
        clustersUnique = sorted(color_cluster_df.loc[: , "ClusterName"].unique())
        hatches = itertools.cycle(hatchOrder)
        hatchDict = { clust: hatch for clust , hatch in zip(clustersUnique , hatches ) }
        color_cluster_df["hatch"] =  color_cluster_df.loc[: ,  "ClusterName"].map( hatchDict )
        return  color_cluster_df
    ################################################
    ## Stetup plotDF
    plotDF = plotDF.loc[: , [z for z in [x, y,"ClusterName" , "color",  asterisk_col] if z is not None  ]].copy()
    plotDF["label_tmp"] = list(range(plotDF.shape[0]))
    if addHatches :
        plotDF =   plotDF.groupby( by = "color" ).apply( lambda x:  addHatchColumn(color_cluster_df = x ) )
    ####################################################
    ## Make plot
    barsPerAx = -((-plotDF.shape[0]) // n_ax )
    fig , axes = plt.subplots(figsize = figsize , nrows = 1,ncols = n_ax )
    axes = np.ravel(axes)
    for i in range(0, n_ax):
       ax = axes[i]
       sns.barplot( y= "label_tmp", x= x , data = plotDF.iloc[i*barsPerAx: (i+1)*barsPerAx, :].copy() ,
                   orient = "h" , ax  = ax)
       ax.set_yticks( list(range(  len(plotDF["GOterm_short"].values[i*barsPerAx: (i+1)*barsPerAx]) ) ))
       ytickLabels = list( plotDF.loc[: , "GOterm_short" ].values[i*barsPerAx : (i+1)*barsPerAx] )
       ax.set_yticklabels( [ '\n'.join(wrap(l, 30)) for l in ytickLabels ] )
    ####################################################
    ## color bars 
    for bar , row in zip( ax.patches , plotDF.iterrows() ):
        row_idx , row = row
        bar.set_facecolor(row["color"])
        if addHatches :
            bar.set_hatch(row["hatch"])
        bar.set_label(row["ClusterName"])
    ################################################
    ## Add Legend. (for loop removes repeat labels)
    handles, labels = ax.get_legend_handles_labels()
    handle_list, label_list = [], []
    for handle, label in zip(handles, labels):
        if label not in label_list:
            handle_list.append(handle)
            label_list.append(label)
    ax.legend(handle_list, label_list, loc = "upper left" , bbox_to_anchor = (1.0 , 1.0)  )
    ##############################################################
    ## add asterisks
    if asterisk_thresholds:
        asterisk_thresholds = sorted( asterisk_thresholds , reverse = True)
        markers = ["*"*i for i in range(1, len( asterisk_thresholds) + 1 )]
        for bar , row in zip( ax.patches , plotDF.iterrows() ):
            row_idx , row = row
            asterisk_thresholds_satisified =  [ idx for idx ,  thresh in enumerate(asterisk_thresholds) if row[asterisk_col]  < thresh]
            if  asterisk_thresholds_satisified:
                marker = markers[max( asterisk_thresholds_satisified )]
                marker_x , marker_y =  get_barPosition( bar  , loc = "rightCenter" if  horizontal else "topCenter",
                                                        pad_top = asterisk_pad if not horizontal else 0.04,
                                                        pad_right = asterisk_pad if horizontal else 0 )
                ax.text( marker_x , marker_y , marker, fontsize = 20  , horizontalalignment = 'center' , verticalalignment = 'center')
    ##########################################################
    ## Tick and axis labels:
    if horizontal:
        ax.set_yticklabels( list( plotDF.loc[:, y] ) , rotation =yticklabel_rot )
    else:
        ax.set_xticklabels( list( plotDF.loc[:, x]  ) , rotation = xticklabel_rot )
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        if not horizontal:
            ax.set_xlabel(x)
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        if horizontal:
            ax.set_ylabel(ylabel)
    return fig 


####### SCRIPT #############################################################################
if __name__ == "__main__":
    ### READ CMD LINE #############################################
    parser = argparse.ArgumentParser( description = "Make bar plot of GO pVals from multiple DAVID output files" )
    parser.add_argument( "-i" , help = "file indicating GO files and terms to plot. Each row is a term, \
    Columns are tab separated and are <Term full>\t<Tem abbreviation>\t<Fname>\t<ClusterName>\t<plotKwargs (comma and  colon separated)>")
    parser.add_argument("--basedir" , help = "path to GO dir")
    parser.add_argument( "-o" ,  help = "output file" )
    parser.add_argument("--addHatches" , help = "Do we use hatches do distingush different clusters assiged same color", action = "store_true")
    parser.add_argument('--title' )
    parser.add_argument( '--figsize' , default = '7.5,4' )
    parser.add_argument('--pVal_col' ,default = "Benjamini")
    parser.add_argument('--n_cols', default = 1, type = int, help = "Number of columns of horizonal barplots")
    parser.add_argument('--asteriskThresh', help = "float or comma separated list of floats")
    args = parser.parse_args()
    
    # exampleArgs= ("-i GO_termsInterest_example.txt --basedir  /Users/afinneg2/projects/keratinocyteRegulators/analysis2/motifAnalysis/coReg_coOccur/coRegGenes.fantom.logFCthresh0.25.average_corrTFs.motifsEnrichedNHEKP/\clustGrid_progenitor_beta4_colorBranches.fine/GO_bkgnd.geq1Pct  -o test.pdf   --addHatches  --title Enriched_terms_for_progenitor_gene_clusters")
    # args = parser.parse_args(exampleArgs.split() ) 
    
    ### LOAD DATA #################################################################################################
    termsInterest_df =pd.read_csv(args.i , sep = "\t" , names = ["GOterm" , "GOterm_short" , "fname", "ClusterName" , "plotKwargs"] )
    print( termsInterest_df)
    termsInterest_df.loc[: ,"fname"] = termsInterest_df.loc[: , "fname"].map( lambda x: os.path.join(args.basedir , x.strip())  )
    ## Parse kwargs column
    print(termsInterest_df)
    plotKwargs = list(set().union( *list(termsInterest_df.loc[:, "plotKwargs"].map( lambda x:  [ y.split(":")[0] for y in x.split(",") ] ).values) )  )
    for kwarg in plotKwargs:
        termsInterest_df[kwarg] =  termsInterest_df.loc[: , "plotKwargs"].map( lambda x: [y.split(":")[1] for y in x.split(",") if  y.split(":")[0]  == kwarg ][0]  )

    ### Get pVals from GOfiles ####################################################################
    termsInterest_df = add_p(fnameTerms_orig = termsInterest_df, fnameCol = "fname", TermCol = "GOterm", pVal_col = args.pVal_col )
    termsInterest_df["-log10p"] = termsInterest_df.loc[: , "p"].map( lambda x: -1*np.log10(x) )    
    ### PLOT #########################################################################################
    if args.title is not None:
        title = args.title.strip('"')
    else:
        title = None
    
    asteriskThresh = [ float(x) for x in args.asteriskThresh.split(",") ] if args.asteriskThresh is not None else None
    print(asteriskThresh )
    fig = make_barPlot( plotDF = termsInterest_df.copy() , y= "GOterm_short" , x = "-log10p" , 
                       addHatches = args.addHatches, title =title ,horizontal = True , 
                       figsize =  tuple([float(x) for x in args.figsize.split(",") ]), n_ax = args.n_cols,
                       asterisk_thresholds = asteriskThresh,  asterisk_col=  "p" , asterisk_pad = 0.05,
                       xlabel = r'$-log_{10}(p)$', ylabel = "GO term")

    fig.tight_layout()
    fig.savefig( args.o ,bbox_inches ='tight' ,format = os.path.splitext(args.o)[-1][1:]   ,dpi = 300 )
    #### Write table of plotted data ###################################################################    
    termsInterest_df.loc[: , ["GOterm" , "GOterm_short" , "ClusterName" , "p", "-log10p" ] ].to_csv(
                                os.path.splitext(args.o)[0] + ".tsv" , sep = "\t" , index = False )
    
    
